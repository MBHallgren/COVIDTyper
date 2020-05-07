#!/usr/bin/env python3

# Copyright (c) 2020, Malte Bjørn Hallgren Technical University of Denmark
# All rights reserved.
#

#Import Libraries

import sys
import os
import argparse
import time
import COVIDTyperFunctions as mtf


parser = argparse.ArgumentParser(description='.')
parser.add_argument('-i_path_illumina', action="store", type=str, dest='i_path_illumina', default="", help='The path to the directory containing ONLY the input illumina files. Should be used when analyzing >5 read-files at a time.')
parser.add_argument('-i_path_nanopore', action="store", type=str, dest='i_path_nanopore', default="", help='The path to the directory containing ONLY the input nanopore files. Should be used when analyzing >5 read-files at a time.')
parser.add_argument('-i_mfa', action="store", type=str, dest='i_mfa', default="", help='Multi fasta assembly input. Only takes one single multi fasta file as input. Please give complete path to the file, not the directory.')
parser.add_argument('-i_path_assemblies', action="store", type=str, dest='i_path_assemblies', default="", help='The path to the directory containing the assembly files')
parser.add_argument("-pe", action="store_true", dest="paired_end", default = False, help="If paipred ends are used give input as True (-pe True). If Paired-ends are used, it is important that the files are written in the correct order, such as: sample1_1.fasta sample1_2.fasta sample2_1.fasta sample2_1.fasta")
parser.add_argument("-bc", action="store", type=float, default = 0.7, dest="bc", help="Base calling parameter for nanopore KMA mapping. Default is 0.7")
parser.add_argument("-db", action="store", type=str, default = "", dest="ref_kma_database", help="Comeplete path for the ref_kma_database for KMA mapping")
parser.add_argument("-thread", action="store", default = 1, dest="multi_threading", help="Set this parameter to x-number of threads that you would like to use during KMA-mapping.")
parser.add_argument("-ref", action="store", type=str, default = "", dest="reference", help="KMA will by default determine the best template against the given database. However, if you want to align your query sequences against a reference of your own choice, use this function. If this is left blank, KMA will determine the optimal reference.")
parser.add_argument('-version', action='version', version='COVIDTyper 1.0.0', help = "current version of COVIDTyper")
parser.add_argument("-exepath", action="store", dest="exepath", default = "", help="Complete path to the COVIDTyper repo that you cloned, in which the executables are located")
parser.add_argument("-o", action="store", dest="output_name", help="Name that you would like the output directory to be called.")
args = parser.parse_args()

def researchPipeline(i_path_illumina, i_path_nanopore, i_mfa, paired_end, bc,
                     ref_kma_database, multi_threading, reference, output_name, exepath, assemblies):
    if assemblies != "":
        assembly_flag = True
        i_path_illumina = assemblies
    else:
        assembly_flag = False

    #if used on server and output path is provided:
    if output_name[0] == "/":
        target_dir = output_name + "cgeout/"
        cmd = "mkdir " + target_dir
        os.system(cmd)
        cmd = "chmod 775 " + target_dir
        os.system(cmd)
        logfilename = target_dir + "logfile"
        logfile = open(logfilename, 'w')
    else:
        current_path = os.getcwd()
        target_dir = current_path + "/" + output_name + "/"
        cmd = "mkdir " + output_name
        os.system(cmd)
        logfilename = target_dir + "logfile_" + output_name
        logfile = open(logfilename, 'w')


    kma_database_path = ref_kma_database
    cmd = "mkdir " + target_dir + "DataFiles"
    os.system(cmd)

    kma_path = exepath + "kma/kma"


    # Print messages
    #Very old logfile, need update
    startTime = time.time()
    print("# Running COVIDTyper 1.0.0 with following input conditions:", file=logfile)
    mtf.logfileConditionsResearch(logfile, bc, ref_kma_database, multi_threading, reference, output_name, paired_end)
    if paired_end == True:
        print("# -pe", file=logfile)
    if bc != 0:
        print("# -bc: " + str(bc), file=logfile)
    if ref_kma_database != "":
        print("# -db: " + ref_kma_database, file=logfile)
    if multi_threading != 1:
        print("# -thread: " + str(multi_threading), file=logfile)
    if reference != "":
        print("# -ref: " + reference, file=logfile)
    print ("loading input")
    print ("MFA: {}".format(i_mfa))

    if i_mfa != "":
        if i_path_illumina != "": #Autodetection on CGE webserver
            i_mfa = i_path_illumina + os.listdir(i_path_illumina)[0]
            print (i_mfa, file = logfile)


        assembly_flag = True
        print ("Multifasta file was given as input.", file=logfile)
        cmd = "mkdir {}singlefastas".format(target_dir)
        os.system(cmd)
        cmd = "chmod 775 " + target_dir + "singlefastas/"
        os.system(cmd)
        try:
            infile = open(i_mfa, 'r')
            sequence = ""
            for line in infile:
                line = line.rstrip()
                if line[0] == ">" and sequence != "":
                    print (sequence, file=outfile)
                    outfile.close()
                    sequence = ""
                    filename = line.split()[0][1:]
                    outfile = open("{}singlefastas/{}".format(target_dir, filename), 'w')
                    print (line, file = outfile)
                elif line[0] == ">" and sequence == "":
                    filename = line.split()[0][1:]
                    outfile = open("{}singlefastas/{}".format(target_dir, filename), 'w')
                    print(line, file=outfile)
                else:
                    sequence += line
            print(sequence, file=outfile)
            outfile.close()
        except IOError as err:
            print ("IOERROR", file=logfile)
            print("Cant open file:", str(err));
            sys.exit(1)
        i_path_illumina = "{}singlefastas/".format(target_dir)

    illumina_files = mtf.load_illumina(i_path_illumina)
    nanopore_files = mtf.load_nanopore(i_path_nanopore)
    illumina_files = mtf.generate_complete_path_illumina_files(illumina_files, i_path_illumina)
    nanopore_files = mtf.generate_complete_path_nanopore_files(nanopore_files, i_path_nanopore)
    total_filenames = mtf.combine_input_files(illumina_files, nanopore_files)

    if 3 > len(nanopore_files) + len(illumina_files):
        sys.exit("You did not supply 2 or more input files. Please run the program again with correct input")


    best_template, templatename = mtf.findTemplateResearch(total_filenames, target_dir, kma_database_path, logfile, reference, kma_path)

    print ("performing KMA mapping")
    if paired_end == True:
        print ("Paired end illumina input was given")
        mtf.illuminaMappingPE(illumina_files, best_template, target_dir, kma_database_path, logfile, multi_threading, reference, kma_path)
    else:
        print ("Only forward illumina reads was given")
        mtf.illuminaMappingForward(illumina_files, best_template, target_dir, kma_database_path, logfile, multi_threading, reference, kma_path)

    mtf.nanoporeMapping(nanopore_files, best_template, target_dir, kma_database_path, logfile, multi_threading, bc, reference, kma_path)

    print ("calculating distance matrix")

    ccphylo_path = exepath + "ccphylo/ccphylo"

    if assembly_flag == True:
        cmd = "{} dist -i {}*.fsa -o {}{} -r \"{}\" -f 9 -mc 1 -nm 0 -nv {}nucleotideVarriance.gz &>> {}distance_matrix_logfile".format(ccphylo_path, target_dir, target_dir, "distmatrix.phy", templatename, target_dir, target_dir)
        os.system(cmd)
    else:
        cmd = "{} dist -i {}*.fsa -o {}{} -r \"{}\" -f 1 -mc 1 -nm 0 -nv {}nucleotideVarriance.gz &>> {}distance_matrix_logfile".format(ccphylo_path, target_dir, target_dir, "distmatrix.phy", templatename, target_dir, target_dir)
        os.system(cmd)


    infile = open("{}distance_matrix_logfile".format(target_dir),'r')
    for line in infile:
        line = line.rstrip()
        print (line, file=logfile)
    infile.close()

    cmd = "rm {}distance_matrix_logfile".format(target_dir)
    os.system(cmd)



    cmd = "{} tree -i {}{} -o {}outtree.newick".format(ccphylo_path, target_dir, "distmatrix.phy", target_dir)
    os.system(cmd)

    mtf.cleanUp(target_dir, illumina_files, nanopore_files, paired_end, reference)
    endTime = time.time()
    dTime = endTime - startTime
    print("COVIDTyper total runtime: " + str(dTime) + " seconds", file=logfile)
    logfile.close()

    cmd = "cat {}DataFiles/*.vcf.gz > {}combined.vcf.gz".format(target_dir, target_dir)
    print (cmd)
    os.system(cmd)
    print ("COVIDTyper has completed")


    mtf.varriansfileRenamer(total_filenames)

def main():
    inputCheck = args.i_path_illumina + args.i_path_nanopore + args.i_mfa + args.i_path_assemblies
    if inputCheck == "":
        sys.exit("No input was given.")

    researchPipeline(args.i_path_illumina, args.i_path_nanopore, args.i_mfa, args.paired_end, args.bc, args.ref_kma_database, args.multi_threading, args.reference, args.output_name, args.exepath, args.i_path_assemblies)
if __name__== "__main__":
    main()


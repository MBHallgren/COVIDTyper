import  os
cmd = "chmod a+x CoronaTyper"
os.system(cmd)
cmd = "chmod a+x CoronaTyperFunctions.py"
os.system(cmd)
#Kma
cmd = "git clone https://bitbucket.org/genomicepidemiology/kma.git"
os.system(cmd)
cmd = "cd kma && make"
os.system(cmd)
cmd = "mv kma/kma .
#ccphylo
cmd = "git clone https://bitbucket.org/genomicepidemiology/ccphylo.git"
os.system(cmd)
cmd = "cd ccphylo && make"
os.system(cmd)
cmd = "mv ccphylo/ccphylo .

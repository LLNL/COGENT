#! /usr/bin/env python
#To run this program type "pythonProgram"
#to make emacs understand python, type <esc>python-mode
import math 
import subprocess
import os
import time
import sys

#parse command line 
if len(sys.argv) != 1:
  print "totalRun reads inputs from input.totalRun"
  sys.exit() 

devNull = open('/dev/null','w')

#read nCell for various file names
base = open("input.totalRun","r")
inFile = base.readlines()

for line in inFile:
  thisLine = line.split()
  if len(thisLine) > 0:
    if thisLine[0] == "log2MaxRes_BrackbillSalzman":
      thisBLine = thisLine
    if thisLine[0] == "nCell":
      thisELine = thisLine  
    if thisLine[0] == "evenSpacingCase":
      thisOnlyELine = thisLine 

BrackbillSalzmanRes = 2**(int(thisBLine[2]))
evenSpacingRes      = int(thisELine[2])
evenSpacingCase     = (thisOnlyELine[2])
 
# convert to boolean
if evenSpacingCase == "true":
  evenSpacingCaseOnly = True
elif evenSpacingCase == "false":
  evenSpacingCaseOnly = False
else:
  print "Must set :evenSpacingCase = true or evenSpacingCase = false. Can't use evenSpacingCase = %s"%(evenSpacingCase)
  sys.exit()

base.close()

#compile initNewBd
os.chdir("newBoundaryNodesIC")
print " Compile new boundary node code"
print "    Compiling: make all DEBUG=FALSE MPI=FALSE OPT=TRUE -j32 DIM=2"
success = subprocess.call(["make","all","DEBUG=FALSE","MPI=FALSE","OPT=TRUE","-j32","DIM=2"],stdout = devNull)
if success != 0:
  print "failed to compile initNewBdSolution"
  sys.exit(success)
os.chdir("..")
  
#compile addMap
os.chdir("utilities")
print " Compile add maps code"
print "    Compiling: make all DEBUG=FALSE MPI=FALSE OPT=TRUE -j32 DIM=2"
success = subprocess.call(["make","all","DEBUG=FALSE","MPI=FALSE","OPT=TRUE","-j32","DIM=2"],stdout = devNull)
if success != 0:
  print "failed to compile addMap"
  sys.exit(success)
os.chdir("..")

#compile Brackbill-Salzman
os.chdir("blocks")  
print " Compile Brackbill-Salzman code"
print "    Compiling: make all DEBUG=FALSE MPI=FALSE OPT=TRUE -j32 DIM=2"
success = subprocess.call(["make","all","DEBUG=FALSE","MPI=FALSE","OPT=TRUE","-j32","DIM=2"],stdout = devNull)
if success != 0:
  print "failed to compile bs.cpp"
os.chdir("..")
 
#compile ghostCellExtrapBS
os.chdir("ghostCellExtrapBS")
print " Compile ghost cell code"
print "    Compiling: make all DEBUG=FALSE MPI=FALSE OPT=TRUE -j32 DIM=2   "
success = subprocess.call(["make","all","DEBUG=FALSE","MPI=FALSE","OPT=TRUE","-j32","DIM=2"],stdout = devNull)
if success != 0:
  print "failed to compile ghostCellExtrapBS"
  sys.exit(success)
os.chdir("..")

#make Results in grid_generator to store data from each block and each part of the algorithm
topLevel = "Results"
if os.path.exists(topLevel):
  timeStamp = '.%02d.%02d'%(time.localtime().tm_mon,time.localtime().tm_mday)
  topLevel += timeStamp
    
if os.path.exists(topLevel):
  hourMin = '.%02d.%02d'%((time.localtime().tm_hour)%12,time.localtime().tm_min)
  topLevel += hourMin

topLevel += "/"

print ""

subprocess.call(["mkdir", "%s"%(topLevel,) ],stdout = devNull)

#change to Results directory
os.chdir(topLevel)
subprocess.call(["mkdir", "utilities"              ],stdout = devNull)
subprocess.call(["mkdir", "utilities/outputMapping"],stdout = devNull)

subprocess.call(["mkdir", "newBoundaryNodesIC"               ],stdout = devNull)
subprocess.call(["mkdir", "newBoundaryNodesIC/outputMapping" ],stdout = devNull)
subprocess.call(["mkdir", "newBoundaryNodesIC/diagnostic"    ],stdout = devNull)

subprocess.call(["mkdir", "ghostCellExtrapBS"              ],stdout = devNull)
subprocess.call(["mkdir", "ghostCellExtrapBS/mappingFiles" ],stdout = devNull)
subprocess.call(["mkdir", "ghostCellExtrapBS/diagnostic"   ],stdout = devNull)
subprocess.call(["mkdir", "ghostCellExtrapBS/outputMapping"],stdout = devNull)
for block in range(0,10):
  print " Creating directory for block%d"          %(block,)

  subprocess.call(["mkdir", "block%d"              %(block,)],stdout = devNull)
  subprocess.call(["mkdir", "block%d/diagnostic"   %(block,)],stdout = devNull)
  subprocess.call(["mkdir", "block%d/checkpoint"   %(block,)],stdout = devNull)
  subprocess.call(["mkdir", "block%d/outputMapping"%(block,)],stdout = devNull)
  subprocess.call(["mkdir", "block%d/mappingFiles" %(block,)],stdout = devNull)

print ""
#run newBoundaryNodesIC

#change to Results/newBoundaryNodesIC
os.chdir("newBoundaryNodesIC")
success = subprocess.call(["ln","-s", "../../newBoundaryNodesIC/newBdInit2d.Linux.64.g++.gfortran.OPT.ex","."],stdout = devNull)
if success != 0:
  print "failed to link executable newBdInit"
  sys.exit(success)
success = subprocess.call(["cp", "../../input.totalRun","."],stdout = devNull)
if success != 0:
  print "failed to copy ../../input.totalRun"
  sys.exit(success)
newBd = "newBdInit2d.Linux.64.g++.gfortran.OPT.ex"
out  = open ("out",'w')

success = subprocess.call([newBd,"input.totalRun"],stdout = out)
if success != 0:
  print "failed to run newBd executable"
  sys.exit(success)
#change to Results/
os.chdir("..")

# symbolic link to executable, initial map,and DCT file in each Result directory
#run Brackbill-Salzman
bs = "bs2d.Linux.64.g++.gfortran.OPT.ex"
even     = "newBoundaryNodesIC/outputMapping/even.curr" 
aligned  = "newBoundaryNodesIC/outputMapping/aligned.curr" 

for block in range(0,10):
  #change to Results/block%d
  os.chdir("block%d"%(block,))
  success = subprocess.call(["ln","-s", "../../blocks/%s"%(bs,),"."],stdout = devNull)
  if success != 0:
    print "failed to link Brackbill-Salzman executable"
    sys.exit(success)
  
  #change to   Results/block%d/mappingFiles 
  os.chdir("mappingFiles")
  success = subprocess.call(["ln","-s", "../../%s.%03d.2d.txt"%(even   ,evenSpacingRes      + 1)  ,"even.curr.%03d.2d.txt"   %(evenSpacingRes      + 1,)],stdout = devNull)
  if success != 0:
    print "failed to link to evenMap in newBoundaryNodesIC/outputMapping"
    sys.exit(success)
  success = subprocess.call(["ln","-s", "../../%s.%03d.2d.txt"%(aligned,BrackbillSalzmanRes + 1)  ,"aligned.curr.%03d.2d.txt"%(BrackbillSalzmanRes + 1,)],stdout = devNull)
  if success != 0:
    print "failed to alignedMap in newBoundaryNodesIC/outputMapping"
    sys.exit(success)
  success = subprocess.call(["ln","-s", "../../../mappingFiles/DCT_coefficients_3-22-16"        ,"."                                                ],stdout = devNull)
  if success != 0:
    print "failed to link to coefficient file ../../../mappingFiles/DCT_coefficients_3-22-16"
    sys.exit(success)
 #change to Results/block%d
  os.chdir("..")

  #change to Results/
  os.chdir("..")

#make directory to put the final map
finalMap = "outputMapping"
success = subprocess.call(["mkdir", finalMap],stdout = devNull)
finalMap += "/"

minRes = "log2MinRes = %d" %(BrackbillSalzmanRes,)
maxRes = "log2MaxRes = %d" %(BrackbillSalzmanRes,)

# remove done file from each block
for block in range(0,10):
  os.chdir("block%d"%(block,))
  success = subprocess.call(["rm","-f","done"],stdout = devNull)
  os.chdir("..")
print ""

for block in range(0,10):
  print "start  Brackbill-Salzman on block%d"%(block,)
  os.chdir("block%d"%(block,))
  
  base = open("../../input.totalRun","r")
  local = open("input","w")
  success = subprocess.call(["sed","-e","s/BLOCK/%d/g"%(block,)],stdin = base,stdout = local)

  local.write(minRes + "\n")
  local.write(maxRes + "\n")
  base.close() 
  local.close()
  
  # run brackbill-salzman simultaneously in each directory
  out  = open ("out",'w')
  subprocess.Popen([bs,"input"],stdout = out)
  os.chdir("..")

while True:
  time.sleep(5)
  finish = True
  for block in range(0,10):
    if not os.path.exists("block%d/done"%(block,)):
      finish = False
      print  "still working on block%d"%(block,)
    else:
     #print  "finish Brackbill Salzman on block%d"%(block,)        
     pass

  print ""
  if finish == True:
    break
      
os.chdir("utilities") 
print "start  combining block maps"
addMap = "addMap2d.Linux.64.g++.gfortran.OPT.ex"
success = subprocess.call(["ln","-s", "../../utilities/%s"%(addMap,),"."],stdout = devNull)
if success != 0:
  print "failed to link to ../../utililities/addMapp...ex"
  sys.exit(success)
success = subprocess.call(["cp", "../../input.totalRun","."],stdout = devNull)
if success != 0:
  print "failed to copyto ../../input.totalRun"
  sys.exit(success)
os.chdir("outputMapping")
if evenSpacingCaseOnly == False:
  print "evenSpacingCaseOnly = False"
  subprocess.call(["ln","-s", "../../../utilities/outputMapping/even.curr.%03d.2d.txt"%(evenSpacingRes + 1,),"."],stdout = devNull)
  if success != 0:
    print "failed to link to ../../../utilities/outputMapping/even.curr"
    sys.exit(success)
os.chdir("..")

out  = open ("out",'w')
subprocess.call([addMap,"input.totalRun"],stdout = out)
if success != 0:
  print "failed to run executable addMap"
  sys.exit(success)
print "finish combining block maps"
print ""

# if only doing the (preliminary step of) evenSpacingCaseOnly = true, copy the result to utilities/mapping/ and exit
if evenSpacingCaseOnly == True:
   success = subprocess.call(["cp","outputMapping/even.curr.%03d.2d.txt"%(evenSpacingRes + 1,),"../../utilities/outputMapping/."],stdout = devNull) 
   if success != 0:
    print "failed to copy outputMapping/even...txt into ../../utilities/outputMapping"
    sys.exit(success)
   print "Writing the map based on evenly spaced boundary nodes to utilities/outputMapping"
   sys.exit()

#change to Results
os.chdir("..")

#print os.getcwd()

os.chdir("ghostCellExtrapBS")
ghostCellExtrap = "ghostCell.BS.10.block2d.Linux.64.g++.gfortran.OPT.ex"
print "start ghost cell extrapolation"
success = subprocess.call(["ln","-s", "../../ghostCellExtrapBS/%s"%(ghostCellExtrap,),"."],stdout = devNull)
if success != 0:
    print "failed to link to executable ../../ghostCellExtrapBS/ghostCell...ex"
    sys.exit(success)
success = subprocess.call(["cp", "../../input.totalRun" ,"."],stdout = devNull)
if success != 0:
    print "failed to copy ../../input.totalRun"
    sys.exit(success)

out  = open ("out",'w')

subprocess.call([ghostCellExtrap,"input.totalRun"],stdout = out)
print "Writing map with ghost cells to %s%s"%(topLevel,finalMap)
success = subprocess.call(["cp","outputMapping/extended_map.%03d.2d.txt"%(BrackbillSalzmanRes + 1,), "../%s/"%(finalMap,)],stdout = devNull)
if success != 0:
    print "failed to copy outputMapping/extended_map...txt"
    sys.exit(success)
print "finsh ghost cell extrapolation"
print ""
os.chdir("..")

#change to grid_generator
os.chdir("..")
#print os.getcwd()

sys.exit()





  






# Copyright (C) <2025>  <The Ohio State University>       

# This program is free software: you can redistribute it and/or modify                              
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or    
# (at your option) any later version.                                                                                       
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of           
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
# GNU General Public License for more details.                                                                             

# You should have received a copy of the GNU General Public License 
# along with this program.  If not, see <https://www.gnu.org/licenses/>;.

import os
import subprocess
from subprocess import Popen,PIPE
import sys

inputChrom = sys.argv[1]

infileConstantsName = "KdandBestSites/Kd"+inputChrom
infileConstants = open(outfileConstantsName, "r")

outfileConstantsName = "KdandDeltaG/KdDG"+inputChrom
outfileConstants = open(outfileConstantsName, "w")
#openfiles


outstring = "seq \t start \t end \t Kd \t BestinSeq \t BestinMotif \t zeroProbCount \t withdeletion \t startD \t endD \t constantD \t BestinSeqD \t BestinMotifD \t indelLoc \t ref \t alt \t deletion \t delta_en \t delta_enD \t ddG \t bindingSite \n"
outfileConstants.write(outstring)
#write header

t = "\t"
s = " "

# check whether this is necessary
if inputChrom == "1":
    infileConstants.readline()

for line in infileConstants:
    fields = line.split()
    seq = fields[0]
    BestinSeq = fields[4]
    BestinMotif = fields[5]
    withdeletion = fields[7]
    BestinSeqD = fields[11]
    BestinMotifD = fields[12]

    if float(BestinSeq.split(":")[1]) - 0.01 <= float(BestinMotif.split(":")[1]) and float(
            BestinSeqD.split(":")[1]) - 0.01 <= float(BestinMotifD.split(":")[1]):
        bindingSite = "inMotif"
    else:
        bindingSite = "outOfMotif"

    proteinStart = int(BestinMotif.split(":")[0])
    proteinEnd = proteinStart + 6
    proteinStartD = int(BestinMotifD.split(":")[0])
    proteinEndD = proteinStartD + 6

    Command1 = "./findEnergy" + s + str(proteinStart) + s + str(proteinEnd) + s + seq
    myPipe1 = Popen(Command1, shell=True, stdout=PIPE)
    out1 = myPipe1.communicate()[0].decode("utf-8")

    Command2 = "./findEnergy" + s + str(proteinStartD) + s + str(proteinEndD) + s + withdeletion
    myPipe2 = Popen(Command2, shell=True, stdout=PIPE)
    out2 = myPipe2.communicate()[0].decode("utf-8")

    ddG = str(float(out2) - float(out1))
    string = "\t".join(fields) + t + out1 + t + out2 + t + str(ddG) + bindingSite + "\n"
    outfileConstants.write(string)

outfileConstants.close()
infileConstants.close()
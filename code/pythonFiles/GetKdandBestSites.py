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

hurPath = " /home/owusu-ansah.11/Desktop/HumanRNA/inputdata/HuR.rnacompete "
inputChrom = sys.argv[1]
#global Variables

infileSequencesName = "sequences150/Sequences"+inputChrom
infileSequences = open(infileSequencesName, "r")

outfileConstantsName = "KdandBestSameSites/Kd"+inputChrom
outfileConstants = open(outfileConstantsName, "w")
#openfiles

outstring = "seq \t start \t end \t Kd \t BestinSeq \t BestinMotif \t zeroProbCount \t withdeletion \t startD \t endD \t constantD \t BestinSeqD \t BestinMotifD \t indelLoc \t ref \t alt \t deletion \n"
outfileConstants.write(outstring)
#write header

t = "\t"

if inputChrom == "1":
    infileSequences.readline()

for line in infileSequences:
    fields = line.split()
    seq = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    withdeletion = fields[3]
    startD = int(fields[4])
    endD = int(fields[5])
    indel = int(fields[6])
    ref = fields[13]
    alt = fields[14]
    deletionSize = len(seq) - len(withdeletion)
    footprint=7

    Command1a = "./findConstant" + hurPath + seq
    myPipe1a = Popen(Command1a, shell=True, stdout=PIPE)
    out1a = myPipe1a.communicate()[0].decode("utf-8")

    Command1b = "./findConstant" + " -c " + out1a + hurPath + seq
    myPipe1b = Popen(Command1b, shell=True, stdout=PIPE)
    out1b = myPipe1b.communicate()[0].decode("utf-8")

    ###############################################################

    Command2a = "./findConstant" + hurPath + withdeletion
    myPipe2a = Popen(Command2a, shell=True, stdout=PIPE)
    out2a = myPipe2a.communicate()[0].decode("utf-8")

    Command2b = "./findConstant" + " -c " + out2a + hurPath + withdeletion
    myPipe2b = Popen(Command2b, shell=True, stdout=PIPE)
    out2b = myPipe2b.communicate()[0].decode("utf-8")

    ###############################################################

    probList1a = [float(i) for i in out1b.split()]
    probList2a = [float(i) for i in out2b.split()]

    #added recently

    summedProbs = []

    for i in range(len(probList2a) - (footprint - 1)):
        if i <= (indel - 1) - (footprint - 1):
            summedProbs.append(probList1a[i] + probList2a[i])
        else:
            summedProbs.append(probList1a[i + (footprint - 1) + deletionSize] + probList2a[i + (footprint - 1)])

    maxinSummedProbs = summedProbs.index(max(summedProbs)) + 1

    if maxinSummedProbs - 1 <= (indel - 1) - (footprint - 1):
        maxinSeq1a = maxinSummedProbs
        maxinSeq2a = maxinSummedProbs
    else:
        maxinSeq1a = maxinSummedProbs + (footprint - 1) + deletionSize
        maxinSeq2a = maxinSummedProbs + (footprint - 1)


    maxInSeqPair1a = str(str(maxinSeq1a) + ":" + str(probList1a[maxinSeq1a - 1]))
    maxInSeqPair2a = str(str(maxinSeq2a) + ":" + str(probList2a[maxinSeq2a - 1]))


    ######################################################################################################
    if start <= indel + deletionSize <= end or start <= indel - footprint + 1 <= end:
        probList1b = probList1a[start - 1:indel - footprint + 1] + probList1a[indel + deletionSize:end]
        probList2b = probList2a[startD - 1:indel - footprint + 1] + probList2a[indel:endD]

    else:
        probList1b = probList1a[int(start) - 1: int(end)]
        probList2b = probList2a[int(startD) - 1: int(endD)]

    actualProbList1b = probList1a[int(start) - 1: int(end)]
    actualProbList2b = probList2a[int(startD) - 1: int(endD)]
    ############################################################################################

    summedProbs2 = []

    for i in range(len(probList1b)):
        summedProbs2.append(probList1b[i] + probList2b[i])

    maxInMotif = summedProbs2.index(max(summedProbs2))

    maxInMotif1b = int(start) + actualProbList1b.index(probList1b[maxInMotif])
    maxInMotifPair1b = str(str(maxInMotif1b) + ":" + str(probList1a[maxInMotif1b - 1]))

    maxInMotif2b = int(startD) + actualProbList2b.index(probList2b[maxInMotif])
    maxInMotifPair2b = str(str(maxInMotif2b) + ":" + str(probList2a[maxInMotif2b - 1]))

    ###############################################################################################

    zerocount = str(probList1a.count(0.0))

    inseq = seq[maxinSeq1a - 1: maxinSeq1a + 5]
    inseq_d = seq[maxinSeq2a - 1: maxinSeq2a + 5]
    motif = seq[maxInMotif1b-1: maxInMotif1b+5]
    motif_d = seq[maxInMotif2b - 1: maxInMotif2b + 5]

    outstring = inseq + t + inseq_d + t + motif + t + motif_d + "\n"

    #outstring = seq + t + str(start) + t + str(end) + t + out1a + t + maxInSeqPair1a + t + maxInMotifPair1b + t + zerocount + t + withdeletion + t + str(startD) + t + str(endD) + t + out2a + t + maxInSeqPair2a + t + maxInMotifPair2b + t + str(indel) + t + ref + t + alt + t + str(deletionSize) + "\n"
    outfileConstants.write(outstring)
    # output data

outfileConstants.close()
infileSequences.close()
#close files

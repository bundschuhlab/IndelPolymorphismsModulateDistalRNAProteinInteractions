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
import copy
import sys

# change directory
# os.chdir('/mnt/c/Users/ok/Desktop/currentStudies/indels/data')

# global variable
window = 50

# argument
chromInput = sys.argv[1]

# open files

infilenameM = "ranges/Ranges" + chromInput
infileM = open(infilenameM, 'r')

infilenamefai = "/home/owusu-ansah.11/Desktop/HumanRNA/inputdata/transcripts.fa.fai"
infilefai = open(infilenamefai, 'r')

outfilenameM = "sequences/Sequences" + chromInput
outfileM = open(outfilenameM, 'w')

tempName = "temporary/Temp" + chromInput
temp = open(tempName, 'w+')


# function to complement rna
def complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna])


# write headers in outfile
t = "\t"

if chromInput == "1":
    outstringheaders = "string" + t + "Start" + t + "End" + t + "withdeltion" + t + "StartD" + t + "EndD" + t + "indelLoc" + t + "strand" + t + "size" + t + "chrom" + t + "startinChrom" + t + "endinChrom" + t + "indelLocinChrom" + t + "ref" + t + "alt" + t + "foldstart" + t + "foldend"
    outfileM.write(outstringheaders + "\n")

# skip header line in infile

next(infileM)

# store info about binding site and indel

for line in infileM:
    fields = line.split()
    chrom = fields[0]
    pos = int(fields[1])
    start = int(fields[6])
    end = int(fields[7])
    seq = ""
    status = fields[4]
    size = int(fields[5])
    ref = fields[2]
    allele = fields[3]
    strand = fields[8]

    for letter in allele:
        iscorrect = (letter == "C" or letter == "G" or letter == "A" or letter == "T")
        if iscorrect is False:
            break

    if iscorrect is False: continue

    for letter in ref:
        iscorrect = (letter == "C" or letter == "G" or letter == "A" or letter == "T")
        if iscorrect is False:
            break

    if iscorrect is False: continue

    if size > 10: continue

    #############################################################################################

    # collect transcript sequence using samtools

    for line in infilefai:
        name = line.split()[0]
        fieldsFai = name.split("|")
        startFai = int(fieldsFai[2])
        endFai = int(fieldsFai[3])
        chromFai = fieldsFai[4]

        if fieldsFai[1] == "1":
            strandFai = "+"
        else:
            strandFai = "-"

        # find transcript in the same chromosome as binding site and indel in the index file
        # binding site and indel location should be within transcript
        # there should be at least 1 window of nucleotides before and after indel. size is the size of insertion/deletion

        if chrom == chromFai and strandFai == strand and startFai < start < end < endFai and startFai < (
                pos - window - size) < pos < (pos + window + size) < endFai:
            break

    infilefai.seek(0)

    # if it has reached end of  file without finding any suitable transcripts
    # skip to the next indel
    if not chrom == chromFai or not strandFai == strand or not startFai < start < end < endFai or not startFai < (
            pos - window - size) < pos < (pos + window + size) < endFai:
        continue

    # find desired sequence in the transcripts file and store in temporary file
    command = '/home/owusu-ansah.11/Desktop/HumanRNA/apps/samtools/bin/samtools faidx /home/owusu-ansah.11/Desktop/HumanRNA/inputdata/transcripts.fa ' + "'" + name + "'"
    retcode = subprocess.run(command, shell=True, stdout=temp)

    # extract sequence from temporary file and clear the file

    temp.seek(0)

    for line in temp:
        if line[0] != ">":
            seq += line[:-1]

    temp.seek(0)

    temp.truncate()

    #################################################################
    realAllele = allele
    realRef = ref

    if strand == "-":
        # the transcript on which it is located will probably be reverse stranded too
        # to avoid starting from the end, we reverse the sequence and treat it like a positive stranded transcript until the end
        # take the complement of the reference and allele sequences in the indels because they are in the reference genome
        # we do not take the reverse of these sequences because they are already compatible with the transcripts
        reverseSeq = seq[::-1]
        refComplement = complement(ref)
        alleleComplement = complement(allele)

        seq = reverseSeq
        ref = refComplement
        allele = alleleComplement
        realAllele = allele[::-1]
        realRef = ref[::-1]

    ###############################################################

    # add 1 because we want the first nucleotide to be a 1

    relativeStart = start - startFai + 1
    relativeEnd = end - startFai + 1
    relativePos = pos - startFai + 1

    # turn the transcript sequence into a list and record its length
    seqlist = [char for char in seq]
    scriptLength = len(seqlist)

    ################################################################

    # if we have a deletion
    if status == "D":

        # cocantenate the items of the ref sequence that are going to be replaced
        # relPos-1 because python is 0 based
        seqlist[(relativePos - 1):(relativePos - 1) + len(ref)] = [
            ''.join(seqlist[(relativePos - 1):(relativePos - 1) + len(ref)])]

        # experiment
        startfold = startFai + (relativePos) - window + len(ref) - 1
        endfold = startFai + (relativePos) + window - 1 + (len(ref) - 1)

        # len(ref) added to lower limit because ref nucleotides indexed together
        # 1 added to upper limit so that length is 2*window
        foldSeq = seqlist[(relativePos - 1) - window + len(ref): (relativePos - 1) + window + 1]

        # add newPos to lower limit of previous line, you get 1 based relative position so this must be correct
        newPos = window - len(ref) + 1

        # allele seq is a copy of original sequence with the ref sequences substituted for the allele
        alleleSeq = copy.deepcopy(foldSeq)
        alleleSeq[newPos - 1] = allele

        # turning lists into strings
        foldstring = "".join(foldSeq)
        withdeletion = "".join(alleleSeq)

        # relativeStart - (1-based start of fold sequence -1) -1 because we want the position of the nucleotide before the start nucleotide
        newStart = relativeStart - (relativePos - window + len(ref) - 1)
        newEnd = relativeEnd - (relativePos - window + len(ref) - 1)

        # what are the start and end positions of the bindings sites in sequences to be folded
        # if both the start and end of the binding site come after the indel

        if relativeStart > relativePos+size:
            newStartD = newStart - size
            newEndD = newEnd - size

            #start is in indel

        elif relativePos <= relativeStart <= relativePos+size:
            newStartD = newStart - (relativeStart-relativePos)
            newEndD = newEnd - size

            # if both the start and end of binding site come before the indel

        elif relativeEnd < relativePos:
            newStartD = newStart
            newEndD = newEnd

            # if relative start comes before indel and relative end comes after indel

        elif relativeStart < relativePos and relativeEnd > relativePos+size:
            newStartD = newStart
            newEndD = newEnd - size

            # end is in indel

        elif relativePos <= relativeEnd <= relativePos+size:
            newStartD = newStart
            newEndD = newEnd - (relativeEnd-relativePos)


    ################################################################################

    # if we have an insertion
    if status == "I":

        # experiment
        startfold = startFai + relativePos - window + len(ref) - 1
        endfold = startFai + (relativePos) + window - size - 1 + (len(ref) - 1)

        # this portion is similar to that of the deletion
        seqlist[relativePos - 1:(relativePos + len(ref) - 1)] = [
            ''.join(seqlist[relativePos - 1:(relativePos + len(ref) - 1)])]
        foldSeq = seqlist[(relativePos - 1) - window + len(ref): (relativePos - 1) + window + 1 - size]
        newPos = (window - len(ref)) + 1

        alleleSeq = copy.deepcopy(foldSeq)
        alleleSeq[newPos - 1] = allele

        # this time the fold string is the allele (with insertions) and the deletion string is the original
        foldstring = "".join(alleleSeq)
        withdeletion = "".join(foldSeq)

        newStartD = relativeStart - (relativePos - window + len(ref) - 1)
        newEndD = relativeEnd - (relativePos - window + len(ref) - 1)

        # if both the start and end of the binding site come after the indel
        if relativeStart > relativePos + size:
            newStart = newStartD + size
            newEnd = newEndD + size

        # start is in indel
        elif relativePos < relativeStart <= relativePos + size:
            newStart = newStartD + size
            newEnd = newEndD + size

        # if both the start and end of the binding site come before the indel
        elif relativeEnd < relativePos:
            newStart = newStartD
            newEnd = newEndD

        # if start comes before indel and end comes after indel
        elif relativeStart < relativePos and relativeEnd > relativePos+size:
            newStart = newStartD
            newEnd = newEndD + size
        ##
        elif relativePos < relativeEnd <= relativePos + size:
            newStart = newStartD
            newEnd = newEndD + size

        # if start of the binding region is equal to the location of the indel
        elif relativeStart == relativePos:
            newStart = newStartD
            newEnd = newEndD + size

        # end of binding region is equal to the location of the indel
        elif relativeEnd == relativePos:
            newStart = newStartD
            newEnd = newEndD

    ##########################################################################

    if strand == "-":
        # reverse fold string and deletion string so that they represent sequences on the transcript
        # the start binding positions become the end binding positions and the ends become the start
        reverseFoldString = foldstring[::-1]
        reverseWithDeletion = withdeletion[::-1]

        newStartReflected = len(foldstring) - newStart + 1
        newEndReflected = len(foldstring) - newEnd + 1
        newStartDReflected = len(withdeletion) - newStartD + 1
        newEndDReflected = len(withdeletion) - newEndD + 1

        if status == "D":
            newPosReflected = (len(foldstring) - newPos + 1) - (len(ref) - 1)
        else:
            newPosReflected = (len(withdeletion) - newPos + 1) - (len(ref) - 1)

        foldstring = reverseFoldString
        withdeletion = reverseWithDeletion
        newStart = newEndReflected
        newEnd = newStartReflected
        newStartD = newEndDReflected
        newEndD = newStartDReflected
        newPos = newPosReflected

    ##########################################################################

    # if the binding site and indel location are in the fold sequence, print out info about the fold sequence and deletion

    if 0 < newStart < newEnd < 2 * window and 0 < newStartD < newEndD < 2 * window - size and 0 < newPos < 2 * window:

        outstring = foldstring + t + str(newStart) + t + str(newEnd) + t + withdeletion + t + str(newStartD) + t + str(
            newEndD) + t + str(newPos) + t + strandFai + t + str(size) + t + chrom + t + str(start) + t + str(
            end) + t + str(
            pos) + t + realRef + t + realAllele + t + str(startfold) + t + str(endfold) + "\n"

        outfileM.write(outstring)

    else:
        continue

infileM.close()
infilefai.close()
outfileM.close()
temp.close()
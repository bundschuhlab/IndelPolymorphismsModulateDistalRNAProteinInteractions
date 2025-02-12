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
import sys

#global variable
window = 200
chromInput = sys.argv[1]

#opening files
infilenameM = "/home/owusu-ansah.11/Desktop/HumanRNA/inputdata/mukherjee/Mukherjee"+chromInput
infileM = open(infilenameM, 'r')

infilenameDel = "/home/owusu-ansah.11/Desktop/HumanRNA/inputdata/indels/Indels" +chromInput
infileDel= open(infilenameDel, 'r')

outfile = "ranges/Ranges"+chromInput

#outfiles = ["ranges/Ranges" +str(i+1) for i in range(22) ] + ["ranges/RangesX", "ranges/RangesY"]

for x in outfiles:
    exec("%s = open(x, 'w')"%("out"+x[7:]))


#writng header
headerString = "chrom del_pos ref alt del/Insert size binding_start binding_end strand \n"

for i in range(22):
    exec("outRanges%s.write(headerString)" % (i + 1))

outRangesX.write(headerString)
outRangesY.write(headerString)


# iterate through indel file store indel info
for lineDel in infileDel:
    if lineDel[0] != "#":
        fieldsDel = lineDel.split()
        chromDel = fieldsDel[0]
        pos = int(fieldsDel[1])
        ref = fieldsDel[3]
        alt = fieldsDel[4].split(",")[0]

        # if ref sequence is greater than allele sequence call it a deletion else call it an insertion
        if len(ref) > len(alt):
            stat = "D"
            size = len(ref) - len(alt)
        else:
            stat = "I"
            size = len(alt) - len(ref)

        # iterate through HuR file (Mukherjee) and store HuR info
        for lineM in infileM:
            fieldsM = lineM.split()
            chromM = fieldsM[0]
            start = int(fieldsM[1])
            end = int(fieldsM[2])
            strand = fieldsM[5]

            # max distance of binding site to indel
            disToDel = max(abs(pos - start), abs(pos - end))

            # if the maximum distance to the indel is less than the window, store HuR info and indel info in the output file
            if chromM == chromDel and disToDel <= window:
                string = "\t".join([chromM, str(pos), ref, alt, stat, str(size), str(start), str(end), strand, "\n"])

                exec("outRanges%s.write(string)" % chromM)

        infileM.seek(0)

#closing files

infileM.close()
infileDel.close()

for i in range(22):
    exec("outRanges%s.close()" % (i + 1))

outRangesX.close()
outRangesY.close()
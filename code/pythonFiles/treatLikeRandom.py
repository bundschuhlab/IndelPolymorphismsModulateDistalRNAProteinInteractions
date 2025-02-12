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

from subprocess import Popen, PIPE
import sys

inputChrom = sys.argv[1]

infileName = "sequences150/Sequences"+inputChrom
infile = open(infileName, "r")

outfileName =  "humanLikeRandom/data"+inputChrom
outfile = open(outfileName, "w")

header = "#\t"
for i in range(101, 145):
    header += str(i) + "\t"
header += "dels"

if inputChrom == "1":
    outfile.write(header)
    infile.readline()

counter = 1

for line in infile:
    fields = line.split()
    seq = fields[0]

    command = "./indel_human seq"
    myPipe = Popen(command, shell=True, stdout=PIPE)
    output = myPipe.communicate()[0].decode("utf-8")

    outfile.write(output)

    if counter == 250:
        break

    counter+=1




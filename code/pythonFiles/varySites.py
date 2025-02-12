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

outfileName =  "varySites/data"+inputChrom
outfile = open(outfileName, "w")

header = "#\t"
for i in range(101, 145):
    header += str(i) + "\t"
header += "dels \n"

if inputChrom == "1":
    outfile.write(header)
    infile.readline()

counter1 = 1
counter = 1

for line in infile:
    fields = line.split()
    seq = fields[0]
    indel_loc = fields[6]
    size = fields[8]

    if size==1 and counter1 > 2000:
        continue


    command = "./varySites " + seq + " " + indel_loc + " " + size
    myPipe = Popen(command, shell=True, stdout=PIPE)
    output = myPipe.communicate()[0].decode("utf-8")

    fields = output.split()
    AllDDGs = fields[:-1]
    deletion = fields[-1]

    if len(AllDDGs) <= 44:
        continue

    DDGs = "\t".join(AllDDGs[:44])+"\t"+deletion

    outfile.write(DDGs + "\n")

    if counter == 10000:
        break

    counter+=1
    if size == 1:
        counter1 += 1




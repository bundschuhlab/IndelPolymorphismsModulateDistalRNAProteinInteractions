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


#opening files
inPath = "/home/owusu-ansah.11/Desktop/HumanRNA/inputdata/"

infilenameM = inPath + "MukherjeeNew.bed"
infileM = open(infilenameM, 'r')

infilenameDel = inPath + 'indels.vcf'
infileDel= open(infilenameDel, 'r')

outPath = "/home/owusu-ansah.11/Desktop/HumanRNA/inputdata/"

outfilesIndels = [outPath + "indels/Indels" +str(i+1) for i in range(22) ] + [outPath+"indels/IndelsX", outPath+"indels/IndelsY", outPath+"indels/IndelsM"]
outIndels = []
outfilesMuk = [outPath + "mukherjee/Mukherjee" +str(i+1) for i in range(22) ] + [outPath+"mukherjee/MukherjeeX", outPath+"mukherjee/MukherjeeY", outPath+"mukherjee/MukherjeeM" ]
outMukherjee = []

#writing headers

headerIndels = infileDel.readline()

for x in outfilesIndels:
    name = "out" + x[(len(outPath) + 7):]
    outIndels.append(name)
    exec("%s = open(x, 'w')" % (name))
    exec("%s.write(str(%s))" % (name, headerIndels))

for x in outfilesMuk:
    name = "out" + x[(len(outPath) + 10):]
    outMukherjee.append(name)
    exec("%s = open(x, 'w')" % (name))
    
#writing files

for line in infileM:
    chrom = line.split()[0]
    linestring = "\t".join(line.split())
    exec("outMukherjee%s.write(linestring)" % chrom)

for line in infileDel:
    chrom = line.split()[0]
    linestring = "\t".join(line.split())
    exec("outIndels%s.write(linestring)" % chrom)
    
#closing files
for i in outIndels:
    exec("%s.close()" % i)

for i in outMukherjee:
    exec("%s.close()" % i)

infileDel.close()
infileM.close()
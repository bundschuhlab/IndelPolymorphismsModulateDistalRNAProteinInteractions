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

from subprocess import call

# file names: DEL_proteinSize_numIterations_omitRegion_seqLength_randomSequenceSeed.dat;

file_names = ["testData/DEL_10_2_50_150_"+str(i+1)+".dat" for i in range(2)]

with open("testData/DEL_10_4_50_150.dat", 'w') as outfile:
    for file_name in file_names:
        with open(file_name) as infile:
            for line in infile:
                outfile.write(line)

for file_name in file_names:
    command = "rm " + file_name
    retcode = call(command, shell=True)
    print(retcode)
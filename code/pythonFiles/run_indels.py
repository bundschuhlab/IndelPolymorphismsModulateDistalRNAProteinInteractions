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

#"usage: ./DEL.x <protein_size> <num_iterations> <seq_length> <random_sequence_seed> <header>" << endl << "seq_length has specifications and header is either 1 or 0: read README ";

command = "nohup ./DEL.x" + " 7 2 150 1 1 &"
retcode = call(command, shell=True)
print(retcode)

for i in range(1):
    command = "nohup ./DEL.x" + " 7 2 150 "  + str(i+2) + " 0 &"
    retcode = call(command, shell=True)
    print(retcode)
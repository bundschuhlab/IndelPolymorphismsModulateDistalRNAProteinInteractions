// Copyright (C) <2025>  <The Ohio State University>       
// 
// This program is free software: you can redistribute it and/or modify                              
// it under the terms of the GNU General Public License as published by 
// the Free Software Foundation, either version 3 of the License, or    
// (at your option) any later version.                                                                                       
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of           
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
// GNU General Public License for more details.                                                                             
// 
// You should have received a copy of the GNU General Public License 
// along with this program.  If not, see <https://www.gnu.org/licenses/>;.

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <sstream>
extern "C"
{
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <string.h>
#include  <data_structures.h>
#include  <params.h>
#include  <utils.h>
#include  <eval.h>
#include  <fold.h>
#include  <part_func.h>
}


using namespace std;

//compile with: g++ -g -fopenmp -I /mnt/c/Users/ok/Desktop/currentStudies/indels/ViennaRNA2/src/ViennaRNA -L /mnt/c/Users/ok/Desktop/currentStudies/indels/ViennaRNA2/lib -o ../executables/findEnergy reduced_indels.cpp -l RNA

void usage(void) {
    cout << "usage: ./findEnergy <protein_start> <protein_end> <sequence>" << endl;
    exit(0);
}


int seq_length; 
int protein_start; 
int protein_end; 
int protein_length;


int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        usage();
    }

    //set parameters
    protein_length = 7;
    

    
    seq_length = strlen(argv[3]);
    protein_start = atof(argv[1]);
    protein_end = atof(argv[2]);

    //cout << " Protein_Length: " << protein_length << ", seq_length: " << seq_length << ", protein_start:" << protein_start << ", protein_end:" << protein_end << endl ;

    //char seq[seq_length + 1] = &argv[4];
    
    char seq[seq_length + 1];
    char* seq_ptr = seq;
    seq_ptr = argv[3];
    seq_ptr[seq_length] = 0;
    
    

    //cout << seq;

    //char* seq_ptr;
    //seq_ptr = seq;
   // seq_ptr[seq_length] = 0;

    //cout << "seqLength: "<< strlen(seq_ptr) << "seq: " << seq_ptr << endl;

    //initialize parameters for folding
    char* mfe_structure = (char*)vrna_alloc(sizeof(char) * (strlen(seq_ptr) + 1));
    char* prob_string = (char*)vrna_alloc(sizeof(char) * (strlen(seq_ptr) + 1));
    vrna_fold_compound_t* vc = vrna_fold_compound(seq_ptr, NULL, VRNA_OPTION_DEFAULT);

    char* prob_string_protein = (char*)vrna_alloc(sizeof(char) * (strlen(seq_ptr) + 1));

    double mfe = (double)vrna_mfe(vc, mfe_structure);
    vrna_exp_params_rescale(vc, &mfe);
    double en = vrna_pf(vc, prob_string);


    //cout << " MFE structure: " << mfe_structure;
    //cout << ", Minimum free energy: " << mfe << ", Free energy: " << en;
    //cout << ", Probability string: " << prob_string << endl;


    vrna_hc_init(vc);

    for (int i = protein_start; i < protein_end + 1; i++)
    {
        vrna_hc_add_up(vc, i, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS);
        //cout << i << endl;
    }

    double en_protein = vrna_pf(vc, prob_string_protein);

    cout << en - en_protein;

    //frees allocated memory
    free(prob_string);
    free(mfe_structure);
    vrna_fold_compound_free(vc);
    free(prob_string_protein);

    
     
    return 0;
}


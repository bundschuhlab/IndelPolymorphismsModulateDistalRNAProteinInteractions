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

//compile with: g++ -g -fopenmp -I /home/owusu-ansah.11/Desktop/randomRNA/ViennaRNA/include/ViennaRNA -I /home/owusu-ansah.11/Desktop/randomRNA/ViennaRNA/include -L /home/owusu-ansah.11/Desktop/randomRNA/ViennaRNA/lib -o varySites varySites.cpp -l RNA
// compile with g++ -g -fopenmp -I /mnt/c/Users/ok/Desktop/currentStudies/indels/ViennaRNA2/src/ViennaRNA -L /mnt/c/Users/ok/Desktop/currentStudies/indels/ViennaRNA2/lib -o varySites varySites.cpp -l RNA


void usage(void) {
    cout << "usage: ./varySites <sequence> <deletion_loc> <size>" << endl;
    exit(0);
}

int seq_length;
int protein_length;
int protein_start;
int alt_protein_start;
int protein_end;
int alt_protein_end;
int alt_seq_length;
int deletion_loc;
int deletion_size;


int del_nts(char* seq_ptr, int deletion_loc, int deletion_size);

int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        usage();
    }

    //set parameters
    protein_length = 7;
    seq_length = strlen(argv[1]);
    deletion_loc = stoi(argv[2]);
    deletion_size = stoi(argv[3]);
    protein_start = deletion_loc + deletion_size;
    protein_end = seq_length - protein_length + 1;


    char seq[seq_length + 1];
    char* seq_ptr = seq;
    seq_ptr = argv[1];
    seq_ptr[seq_length] = 0;

    char alter[seq_length + 1];
    strcpy(alter, seq_ptr);
    char* alter_ptr;
    alter_ptr = alter;
    alter_ptr[seq_length] = 0;


     //////////////////////////////////////////////////////////////////
    
        
        strcpy(alter, seq_ptr);
        del_nts(alter_ptr, deletion_loc, deletion_size);

        alt_seq_length = strlen(alter_ptr);
        alt_protein_start = deletion_loc;
        alt_protein_end = alt_seq_length - protein_length + 1;

        cout << seq_ptr << endl << alter_ptr << endl;

        //initialize parameters for folding
        char* mfe_structure = (char*)vrna_alloc(sizeof(char) * (strlen(seq_ptr) + 1));
        char* prob_string = (char*)vrna_alloc(sizeof(char) * (strlen(seq_ptr) + 1));
        vrna_fold_compound_t* vc = vrna_fold_compound(seq_ptr, NULL, VRNA_OPTION_DEFAULT);

        char* mfe_structure_DEL = (char*)vrna_alloc(sizeof(char) * (strlen(alter_ptr) + 1));
        char* prob_string_DEL = (char*)vrna_alloc(sizeof(char) * (strlen(alter_ptr) + 1));
        vrna_fold_compound_t* vc_DEL = vrna_fold_compound(alter_ptr, NULL, VRNA_OPTION_DEFAULT);

        char* prob_string_protein = (char*)vrna_alloc(sizeof(char) * (strlen(seq_ptr) + 1));
        char* prob_string_protein_DEL = (char*)vrna_alloc(sizeof(char) * (strlen(alter_ptr) + 1));

        double mfe = (double)vrna_mfe(vc, mfe_structure);
        vrna_exp_params_rescale(vc, &mfe);
        double en = vrna_pf(vc, prob_string);

        double mfe_DEL = (double)vrna_mfe(vc_DEL, mfe_structure_DEL);
        vrna_exp_params_rescale(vc_DEL, &mfe_DEL);
        double en_DEL = vrna_pf(vc_DEL, prob_string_DEL);


        for (int x = protein_start; x <= protein_end; x += 1)
        {
            vrna_hc_init(vc);
            vrna_hc_init(vc_DEL);

            int diffB = alt_protein_start - protein_start;

            for (int i = x; i < (x + protein_length); i++)
            {
                vrna_hc_add_up(vc, i, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS);
                vrna_hc_add_up(vc_DEL, i + diffB, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS);
            }

            double en_protein = vrna_pf(vc, prob_string_protein);
            double en_protein_DEL = vrna_pf(vc_DEL, prob_string_protein_DEL);

            double ddG = en_DEL + en_protein - en - en_protein_DEL;

            cout << "  " << ddG << "  ";
        }

        cout << deletion_size << endl;

        free(prob_string);
        free(mfe_structure);
        vrna_fold_compound_free(vc);
        free(prob_string_DEL);
        free(mfe_structure_DEL);
        vrna_fold_compound_free(vc_DEL);
        free(prob_string_protein);
        free(prob_string_protein_DEL);
    

    return 0;
}


//function deletes nucleotides from the mid-section of the sequence. seq_ptr is a pointer to sequence, k is number of deletions 

int del_nts(char* seq_ptr, int deletion_loc, int deletion_size)
{
    for (int i = 0; i <= seq_length - (deletion_loc+deletion_size); i++)
    {
        seq_ptr[(protein_start - deletion_size + i) - 1] = seq_ptr[(protein_start + i) - 1];
    }

    for (int i = seq_length - deletion_size; i < seq_length; i++)
    {
        seq_ptr[i] = 0;
    }

    return 1;
}

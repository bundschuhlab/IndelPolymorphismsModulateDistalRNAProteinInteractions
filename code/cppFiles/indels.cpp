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

//compile with: g++ -I${HOME}/ViennaRNA/include -L${HOME}/ViennaRNA/lib -o DEL.x indels.cpp -lRNA

void usage(void) {
  cout << "usage: ./DEL.x <protein_size> <num_iterations> <random_sequence_seed>" << endl;
}

int rand_seq(int length, char* sequence);

int del_nts(char* seq_ptr, int k);

int seq_length; //length of the sequences to be folded
int alt_seq_length; // length of the sequence after a deletion
int const omit_region = 50;//region of the sequences where no binding happens. This is the region where deletions happen.
int nts_omitted[15]; // A list containing the types of deletions for each random sequence
int sub_seq; // length of portion of sequence where protein binding happens
int protein_startB; // first location protein will be bound on 2nd subsequence
int alt_protein_startB; // first location protein will be bound on 2nd subsequence of altered sequence
int protein_endB; //last location protein will be bound on 2nd subsequence
int alt_protein_endB; //last location protein will be bound on 2nd subsequence of altered sequence
int num_iterations; //number of sequences to be tested
int protein_length; //length of protein footprint
int nts_omissions; //number of bases omitted from mid-section
int random_sequence_seed;
bool header; // whether the file should contain a header or not


int main(int argc, char *argv[]) 
{
    if(argc != 6) 
    {
    cout << "usage: ./DEL.x <protein_size> <num_iterations> <seq_length> <random_sequence_seed> <header>" << endl << "seq_length has specifications and header is either 1 or 0: read README ";
    exit(0);
    }

  clock_t clocks; //record time to see how long it takes
  clocks = clock();
  

  //set parameters
  protein_length = atof(argv[1]);
  num_iterations = atof(argv[2]);
  seq_length = atof(argv[3]);
  random_sequence_seed = atof(argv[4]);
  header = atof(argv[5]);

  srand(random_sequence_seed); //set seed for random sequences

  for (int i = 1; i <= 5; i++)
      nts_omitted[i - 1] = i;

  for (int i = 1; i <= 9; i++)
      nts_omitted[i + 4] = 5*(i+1); //set the types of deletions that would occur

  sub_seq = seq_length - omit_region; //set length subsequence
  protein_startB = seq_length - (sub_seq/2) + 1; 
  protein_endB = seq_length - protein_length + 1;



  cout << " Protein_startB: " << protein_startB << ", Protein_endB: " << protein_endB << endl << endl;

  //create an output file
  ostringstream my_stream;
  my_stream << "/mnt/c/Users/ok/Desktop/currentStudies/indels/testData/DEL_" << protein_length << "_" << num_iterations << "_"<< omit_region <<"_"<<seq_length<<"_"<<random_sequence_seed<<".dat";
  ofstream DEL_out;
  DEL_out.open(my_stream.str().c_str());

  

 

  if(header)
  {
      //label columns of output file
      DEL_out << "#  ";

      for (int y = protein_startB; y <= protein_endB; y++) {
          DEL_out << "  " << y << "  ";
      }

      DEL_out << "  " << "dels" << "  " << "iteration" << "  ";

      DEL_out << endl;
  }

  //for each iteration generate a sequence and fold it with protein in each possible position
  for(int iteration = 1; iteration<=num_iterations; iteration++) {

    //generate random sequence
    
    char seq[seq_length+1]; //allocate space for original sequence. I put in 7 because when I allocate less space, i get an error. I am not sure why.
    char *seq_ptr; //define a pointer to the original sequence
    seq_ptr = seq; //set the pointer to point to the sequence
    rand_seq(seq_length, seq_ptr); //generate a the original sequence using the random sequence generator

    char alter[seq_length + 1]; //allocate space for the sequence with deleted nucleotides
    strcpy(alter, seq_ptr); // copy the original sequence to this memory location
    char* alter_ptr; // define a pointer to the altered sequence (initially the same as the original sequence)
    alter_ptr = alter; // set it to point to the altered sequence
    

 

    cout << "Sequence: " << seq_ptr << ", Sequence_length: " << strlen(seq_ptr) << ", Sequence_to_be_altered: " << alter_ptr << ", Length: " << strlen(alter_ptr) << endl << endl;

    for (int i=0; i <=13 ; i++) //create sequences with ommissions of different sizes
    {
        nts_omissions = nts_omitted[i]; //specify the size of the omission
        strcpy(alter, seq_ptr);  // copy the original sequence to the memory location of the altered sequence
        del_nts(alter_ptr, nts_omissions); //delete nucleotides from the mid-section of the altered sequence

        alt_seq_length = strlen(alter_ptr); //compute the new length of the altered sequence
        alt_protein_startB = alt_seq_length - (sub_seq / 2) + 1; //define the protein binding start position for the altered sequence
        alt_protein_endB = alt_seq_length - protein_length + 1; // define the protein binding end position for the altered sequence

        cout << "Sequence after " << nts_omissions << " omissions: " << alter_ptr; //printing out intermediary results
        
        cout << ", Length: " <<alt_seq_length << ", Protein_startB: "<< alt_protein_startB << ", Protein_endB: "<< alt_protein_endB << endl << endl; //printing out intermediary results


        //initialize parameters for folding
        char* mfe_structure = (char*)vrna_alloc(sizeof(char) * (strlen(seq_ptr)+1)); //allocate memory for mfe structure of original sequence
        char* prob_string = (char*)vrna_alloc(sizeof(char) * (strlen(seq_ptr)+1));  //allocate memory for probability string of original sequence
        vrna_fold_compound_t* vc = vrna_fold_compound(seq_ptr, NULL, VRNA_OPTION_DEFAULT); //struct that contains info on rna sequence

        char* mfe_structure_DEL = (char*)vrna_alloc(sizeof(char) * (strlen(alter_ptr)+1)); //allocate memory for mfe structure of altered sequence
        char* prob_string_DEL = (char*)vrna_alloc(sizeof(char) * (strlen(alter_ptr)+1)); //allocate memory for probability string of altered sequence
        vrna_fold_compound_t* vc_DEL = vrna_fold_compound(alter_ptr, NULL, VRNA_OPTION_DEFAULT); //struct that contains info on altered rna sequence

        char* prob_string_protein = (char*)vrna_alloc(sizeof(char) * (strlen(seq_ptr)+1)); //allocate mem for probability string after protein binding of original sequnce
        char* prob_string_protein_DEL = (char*)vrna_alloc(sizeof(char) * (strlen(alter_ptr)+1)); //allocate mem for probability string after protein binding of altered sequnce

        double mfe = (double)vrna_mfe(vc, mfe_structure); // returns minimum free energy of original sequence. mfe stucture contains mfe rna
        vrna_exp_params_rescale(vc, &mfe); // rescaling vc of original sequence using minimum free energy
        double en = vrna_pf(vc, prob_string); //returns ensemble free energy of original sequence. prob_string is modified to contain probabilities of base pairing

        double mfe_DEL = (double)vrna_mfe(vc_DEL, mfe_structure_DEL); // returns minimum free energy of altered sequence. mfe stucture contains mfe rna
        vrna_exp_params_rescale(vc_DEL, &mfe_DEL); // rescaling vc of altered sequence using minimum free energy
        double en_DEL = vrna_pf(vc_DEL, prob_string_DEL); //returns ensemble free energy of altered sequence. prob_string is modified to contain probabilities of base pairing

        cout <<  " MFE structure: " << mfe_structure << ", After deletion: " << mfe_structure_DEL; //prints out intermediate results
        cout << ", Minimum free energy: " << mfe << ", After deletion: " << mfe_DEL << ", Free energy: " << en << ", After deletion: " << en_DEL;
        cout << ", Probability string: " << prob_string << " , After deletion: " << prob_string_DEL << endl;

    

        for (int x = protein_startB; x <= protein_endB; x+=1)
        {
            vrna_hc_init(vc); // initialize propetrites of vc that deal with hard constraints for original sequence
            vrna_hc_init(vc_DEL); // initialize propetrites of vc that deal with hard constraints for the altered sequence

            int diffB = alt_protein_startB - protein_startB; //difference in the location of the first protein binding site between the altered sequence and the original sequence. To determine which nucleotides to constrain

            for (int i = x; i < (x + protein_length); i++) 
            {
                vrna_hc_add_up(vc, i, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS); //Introducing constraint no-base-pair at protein location in original sequence
                vrna_hc_add_up(vc_DEL, i + diffB, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS); //Introducing contraint no-base-pair at protein location in altered sequence

                cout << "Base " << i << " is " << i + diffB << " in altered sequence, "; // printing intermediate results
            }

            double en_protein = vrna_pf(vc, prob_string_protein); //returns ensemble free energy of original sequence after protein binding
            double en_protein_DEL = vrna_pf(vc_DEL, prob_string_protein_DEL); //returns ensemble free energy of altered sequence after protein binding


            // (en_protein-en)-(en_protein_DEL-en_DEL)
            double ddG = (en_DEL - en_protein_DEL) - (en - en_protein); //The change in this change (change in free energy due to no_paring constraint on protein binding site) due to an omission.
           

            cout << "Protein Bound at " << x << " and at " << x + diffB << " in altered sequence" << " -  Probability string: " << prob_string_protein << " , After deletion: " << prob_string_protein_DEL; //printing out intermediate results
            cout << " , Free Energy: " << en_protein << ", After deletion: " << en_protein_DEL << ", Change in free energy: " << ddG << endl;

            DEL_out << "  " << ddG << "  ";
        }

        cout << endl;

        DEL_out << "  " << nts_omissions << "  " << iteration <<"  "<< endl;

        //frees allocated memory
        free(prob_string);
        free(mfe_structure);
        vrna_fold_compound_free(vc);
        free(prob_string_DEL);
        free(mfe_structure_DEL);
        vrna_fold_compound_free(vc_DEL);
        free(prob_string_protein);
        free(prob_string_protein_DEL);
    }

    clocks = clock() - clocks;
    double diff = ((double)clocks) / CLOCKS_PER_SEC;
    cout << "indels_op5: " << diff << " seconds" << endl;

    
    }

  return 0;
}

//Functions

//function that generates a random sequence. Arguments are the length of the sequence and a pointer to the location where it should be stored.

int rand_seq(int length, char* sequence)  
{

  for(int s=0; s<length; s++) {

    double r = ((double)rand()/(RAND_MAX));

    if(r<.25)
      sequence[s] = 'A';
    if(r<.5 && r>=.25)
      sequence[s] = 'T';
    if(r<.75 && r>=.5)
      sequence[s] = 'C';
    if(r<1 && r>=.75)
      sequence[s] = 'G';
  }


  sequence[length] = 0;

  return 1;
}


//function deletes nucleotides from the mid-section of the sequence. seq_ptr is a pointer to sequence, k is number of deletions 

int del_nts(char* seq_ptr, int k) 
{
    for (int i = 0; i < (sub_seq / 2); i++)
    {
        seq_ptr[protein_startB + i - k - 1] = seq_ptr[protein_startB + i - 1];
    }

    for (int i = seq_length - k; i < seq_length; i++)
    {
        seq_ptr[i] = 0;
    }

    return 1;
}

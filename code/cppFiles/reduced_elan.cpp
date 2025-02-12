
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
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


#include <cmath>
#include <string>
#include <vector>
#include <time.h>
#include <algorithm>

extern "C"{
#include "stdlib.h"
#include "stdio.h"
#include <cstring>
#include "part_func.h"
#include "fold.h"
#include "fold_vars.h"
#include "read_epars.h"
}

using namespace std;

//compile with g++ -o findConstant -fopenmp -Wall reduced_elan.cpp -I../ViennaRNA-2.0.7_revise/H -L../ViennaRNA-2.0.7_revise/lib -lm -lgomp -lRNA

double kT;
double kT_RNAC; //kT in RNAcompete experiment

void k_mer_list(char* , char * , double, double); 
double RBPprofile(char *  , int , int, double);


void usage(void) {
  cout << "usage: ./RNAbind [-c <concentration>] <Protein file name> <sequence>\n";
}

int main(int argc, char * argv[])
{
  int FirstArg=1;
  float Concentration=-1.0;

  // interpret command line arguments: assign concentration to concentration variable and make the protein file the firstarg
  for(;FirstArg<argc && argv[FirstArg][0]=='-';FirstArg++) {
    switch(argv[FirstArg][1]) {
    case 'c':
      if (++FirstArg>=argc) {
        usage();
        return(20);
      }
      Concentration=atof(argv[FirstArg]);
      if (Concentration<0.0) {
        usage();
        return(21);
      }
      break;
    default:
      usage();
      return(10);
    }
  }

    if(argc != FirstArg+2)
    {
      usage();
      return 1;
    }

 
    temperature = 37;  //atof(argv[1]);//23.5;
    double temperature_RNAC = 37;//temperature in RNAcompete experiment 
    footprint = 7;
 
    kT = (temperature+273.15)*1.98717/1000.;
    kT_RNAC = (temperature_RNAC+273.15)*1.98717/1000.;
 
    pf_scale = 1.;
    do_backtrack = 0; //do not calculate base pair binding probability to save CPU time
    

    string sequence = argv[FirstArg+1]; //sequence is stored here
    string RBPname;    

    int seqlen = sequence.length();
    
    char seq[seqlen+1];
    strcpy(seq, sequence.c_str());

    //converting to RNA sequence
    for(int i = 0 ; i < seqlen ; i++)     
      {
	if(seq[i] == 'a')
	  seq[i] = 'A';
	if((seq[i] == 't') || (seq[i] == 'T') || (seq[i] == 'u'))
	  seq[i] = 'U';
	if(seq[i] == 'c')
	  seq[i] = 'C';
	if(seq[i] == 'g')
	  seq[i] = 'G';
      }
    
    //store scores for each 7mer in sequence in InvKd
    k_mer_list(argv[FirstArg], seq, 20., 18.2/1.2);//20nM, 18.2\mu g/1.2 \mu g
    
    RBPtype = 1;
     
    if (Concentration<0.0) {

     // calculate Kd and binding curve

     //=============== for Kd_{Vienna+P} ===============     
      
     //partition function when concentration is 0

      double C_a = 0;
      C_NC[1] = 0;//0nM
      double pf0 = pf_fold(seq, NULL);


      //partition function at concentration beyond Kd
      double C_b = 100;
      double pfc;
      double ratioc;
      do {
	        C_b*=10;
	        C_NC[1] = C_b;
	        pfc = pf_fold(seq, NULL);
	        ratioc = 1-exp(-(pf0-pfc)/kT);
        } while(ratioc<0.5);
      
      //partition function at middle     
      double C_c = (C_a+C_b)/2.;
      C_NC[1] = C_c;   
      pfc = pf_fold(seq, NULL);
      ratioc = 1-exp(-(pf0-pfc)/kT);
      
      //bisection solver
      while( abs(ratioc-0.5) > 0.0001 )
        {   
          if(ratioc > 0.5)
            {   
              C_b = C_c;
              
              C_c = (C_a + C_b)/2.;
              C_NC[1] = C_c;
              pfc = pf_fold(seq, NULL);
              ratioc = 1-exp(-(pf0-pfc)/kT); 
            } 
          else
            {
              C_a = C_c;
	    
              C_c = (C_a + C_b)/2.;
              C_NC[1] = C_c;
              pfc = pf_fold(seq, NULL);
              ratioc = 1-exp(-(pf0-pfc)/kT);
            }
        }
      
      cout << C_NC[1];
      
      // =============== Kd_{Vienna+P} end ===============  

    } else {
      
      // calculate protein binding profile

      C_NC[1]=Concentration;
      double pfbase = pf_fold(seq, NULL);

      for(int i = 1; i <= seqlen-footprint+1; i++) {
        InvKd[1][i] *= 1.01;
        cout << " " << (pfbase-pf_fold(seq,NULL))/kT/0.01;
        InvKd[1][i] /= 1.01;
      }

    } // end of: want binding curve or protein binding profile?
      
    return 0;
}


void k_mer_list(char* Motif_list, char * Realseq, double concentration, double alpha )
{ 
    ifstream Motiflist;
    //Motiflist.open(Motif_list.c_str());
    Motiflist.open(Motif_list);
    if(!Motiflist.is_open())
    {
        cout << "Fail to read Protein file: " << Motif_list << endl;
        return;
    }
    
    char deletehead[100];
    Motiflist.getline(deletehead,100);
    Motiflist.getline(deletehead,100);
           
    vector <string> motif;  //all 7-mer motifs
    vector <double> inverseKd; //corresponing scores
    
    string motiftemp;
    
    while(Motiflist >> motiftemp /*.good()*/)
    {
          
  	footprint = motiftemp.size();
     
   	//change all U to T for searching in motif list
    	for(int i = 0; i < footprint ; i++)
    	{
            if(motiftemp[i] == 'T')
                motiftemp[i] = 'U';
        }
        
        motif.push_back(motiftemp);
           
        double to_ratio;
        Motiflist >> to_ratio;  
        inverseKd.push_back(to_ratio); //store 1/K_{d,I}^{(0)}  

    }
    
    //int motif_list_length = motif.size();
    
    int seqlen = (int)strlen(Realseq);
    string * fraglist = new string[seqlen-footprint+2]; //store all k-mer constituents
    //double scorelist[seqlen-footprint+1]; //store corresponding zscores
    
    

    for(int i = 0; i < seqlen-footprint+1; i++)
    {
        string moving_read; //read each 7-mer constituents
        for(int j = i; j < i+footprint; j++)
        {
            moving_read.push_back(Realseq[j]);
        }
        fraglist[i+1] = moving_read;
        
            unsigned int search_motif = 0;
            while( search_motif < motif.size() &&
		   moving_read.compare(motif[search_motif]) != 0 )
            {
                search_motif += 1;
            }
	    if ( search_motif >= motif.size() ) {
	      cerr << "k-mer " << moving_read << " could not be found." << endl;
	    } else {
              InvKd[1][i+1] = inverseKd[search_motif];
	    }
            
    }


  delete[] fraglist;

  return ;
}



double RBPprofile(char * seq , int RBPsite, int Bsite, double F1)
{
    //double F1 = pf_fold(seq,NULL);
    double InvKd_ori = InvKd[RBPsite][Bsite];
    InvKd[RBPsite][Bsite] = 0;
    double F2 = pf_fold(seq,NULL);
    InvKd[RBPsite][Bsite] = InvKd_ori;
    double Profile = 1-exp(-(F2-F1)/kT);

    return Profile;
}



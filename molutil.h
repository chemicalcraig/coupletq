#ifndef MOLUTIL_H
#define MOLUTIL_H

/**************************************
 * This file contains routines
 * to parse NWChem output files
 * ***********************************/

#include <cstring>
#include <fstream>
#include <iostream>
#include "atom.h"
#include "molecule.h"
#include <iomanip>
using namespace std;

/* Search strings */
string geom_str=" Output coordinates";
string basislabel_str="    Basis function labels";
string overlap_str = " global array: AO ovl";
string nbasis_str="          This is a Direct SCF calculation.";
string tddft_str="  Convergence criterion met";
string tddft_state_stop_str="-----------------------------";
string tddft_task_str1="tddft";
string tddft_task_str2="TDDFT";
string tddft_task_str3="Tddft";
string mo_analysis_str="                       DFT Final Molecular Orbital Analysis";
string lindep_str = " !! The overlap matrix has";
string aobas_str = "          AO basis - number of functions:";
string com_str = " center of mass";
string mass_str = "      Atomic Mass";
/* Number of columns for printing of matrices */
int ncol = 6;

/**** Skip several lines in input file ****/
void getlines(ifstream &in, char *temp, int n, int length) {
  for (int i=0; i<n; i++) {
    in.getline(temp,length);
  }
}


/** Get TDDFT information from logfile **/
bool getTDDFT(string str, Molecule *mol) {
    ifstream infile;
    infile.precision(9);
    infile.open(str.c_str());

    char tempc[1000];
    int natoms=0;
    int nbasis=0;
    bool tddftstack = true;

    while(infile.getline(tempc, 1000)) {
      string temps(tempc);
     
    /** Get Stuff from TDDFT directive **/
      if ((temps.compare(0,5,tddft_task_str1,0,5) == 0 ||
        temps.compare(0,5,tddft_task_str2,0,5) == 0 ||
        temps.compare(0,5,tddft_task_str3,0,5) == 0) &&
        tddftstack) {
        
        while (temps.compare(0,3,"end",0,3) !=0 ) {
          //get nroots
          if (temps.compare(0,6,"nroots")==0 ||
            temps.compare(0,6,"NROOTS")==0 ||
            temps.compare(0,6,"Nroots")==0 ||
            temps.compare(0,6,"NRoots")==0) {
          
            temps = strtok(tempc," ");
            temps = strtok(NULL," ");
            mol->setnroots(atoi(temps.c_str()));
            tddftstack = false;
          } //end nroots

          //get spin state
          if ((temps.compare(0,9,"notriplet")==0) || 
              (temps.compare(0,9,"Notriplet")==0) ||
              (temps.compare(0,9,"NOTRIPLET")==0) ||
              (temps.compare(0,9,"NoTriplet")==0)) {
            mol->spinstate = 0; 
            cout<<"Spin state = "<<mol->spinstate<<endl;
          }

          if ((temps.compare(0,9,"nosinglet")==0) || 
              (temps.compare(0,9,"Nosinglet")==0) ||
              (temps.compare(0,9,"NOSINGLET")==0) ||
              (temps.compare(0,9,"triplet")==0)   ||
              (temps.compare(0,9,"NoSinglet")==0)) {
            mol->spinstate = 1;
            cout<<"Spin state = "<<mol->spinstate<<endl;
          }

          infile>>ws;
          infile.getline(tempc,1000);
          temps = tempc;
        }
      } //end tddft roots
     

      //Grab COM while we're here
      if (temps.compare(0,15,com_str,0,15) == 0) {
        getlines(infile,tempc,2,1000);
        temps = strtok(tempc,"=");
        temps = strtok(NULL," ");
        mol->com[0] = atof(temps.c_str());
        mol->icom[0] = mol->com[0];
        for (int i=0; i<3; i++) temps = strtok(NULL," ");
        mol->com[1] = atof(temps.c_str());
        mol->icom[1] = mol->com[1];
        for (int i=0; i<3; i++) temps = strtok(NULL," ");
        mol->com[2] = atof(temps.c_str());
        mol->icom[2] = mol->com[2];
      }

      //Get atomic masses while we're here
      if (temps.compare(0,17,mass_str,0,17) == 0) {
        
        infile>>ws;
        getlines(infile,tempc,3,1000);

        while(temps.compare(0,10,"Effective nuclear",0,10) != 0) {
          string temps2 = strtok(tempc," ");
          temps = strtok(NULL," ");
          mol->setAtomicMasses(temps2,temps);
          infile>>ws;
          infile.getline(tempc,1000);
          temps = tempc;
        }
      }
      
      //Get TDDFT information
      if (temps.compare(0,10,tddft_str,0,10) == 0) {
        mol->allocateMemTddft();
        
        //Get the ground state energy
        getlines(infile,tempc,2,1000);
        temps = strtok(tempc," ");
        for (int i=0; i<3; i++) temps = strtok(NULL," ");

        mol->groundenergy = atof(temps.c_str());
        getlines(infile,tempc,2,1000);

        for (int root=0; root<mol->nroots; root++) {

          infile.getline(tempc,1000);
          temps = strtok(tempc," ");
          
          for (int i=0; i<4; i++) temps = strtok(NULL," ");
          //get energy of root
          if (root == mol->target-1) {
            cout<<"Root "<<mol->target<<" has an energy "<<atof(temps.c_str())<<endl;
            //CTC change this for multiple excited states on each molecule
            //mol->excenergy[root+1] = atof(temps.c_str());
            mol->excenergy[1] = atof(temps.c_str());
          }

          getlines(infile,tempc,2,1000);
          //skip for triplet
        
          cout<<root<<" "<<mol->nroots<<endl;
        if (mol->spinstate == 0) {
          for (int k=0; k<3; k++) {
            temps = strtok(tempc," ");
            temps = strtok(NULL," ");
      
            for (int j=0; j<3; j++) {
              for (int i=0; i<2; i++) temps = strtok(NULL," ");
                mol->transmoment[k+j*3+root*9] = atof(temps.c_str());
                cout<<"getting transition moment "<<atof(temps.c_str())<<endl;
                /** set initial dipole moment **/
                //if (k==0) {
                //  mol->idip[j] = atof(temps.c_str());
                //  mol->dip[j] = atof(temps.c_str());
                //}
            }

            infile.getline(tempc,1000);
          }//end trans moment
          //get oscillator strength
          temps = strtok(tempc," ");
          for (int i=0; i<3; i++) temps=strtok(NULL," ");
          mol->oscstrength[root] = atof(temps.c_str());
          infile.getline(tempc,1000);
        } else {
          infile.getline(tempc,1000);
          infile.getline(tempc,1000);
        }

        infile.getline(tempc,1000);
        temps = tempc;
        
        //get CI coeffs
        //nrpa is for RPA
        int nrpa = 1;     
        double ycoeff;
        double sumci = 0.;
        double sumci2 = 0.;
        while (temps.compare(0,10,tddft_state_stop_str,0,10) != 0 &&
              temps.compare(0,10,"Target root",0,10) !=0 ) {
          temps = strtok(tempc," ");
          temps = strtok(NULL," ");
          int row = atoi(temps.c_str())-1;
          for (int i=0; i<4; i++) temps = strtok(NULL," .");
          int col = atoi(temps.c_str())-1;
          for (int i=0; i<2; i++) temps = strtok(NULL," ");
          
          nrpa++;
          infile>>ws;
          infile.getline(tempc,1000);
          temps=tempc;
        } //end CI coeffs
   
        } //end nroots
    }//end tddft
  //return true;
  }
}

/************************************
 * Get Number of Atoms
 * ********************************/
/* Get Natoms */
int getNatoms(string filename, int nmol, Molecule *mol) {
  ifstream infile;
  infile.open(filename.c_str());
  char tempc[1000];
  infile.getline(tempc,1000);
  string temps = strtok(tempc,": ");
  temps = strtok(NULL,": ");
  infile.close();

  return atoi(temps.c_str());
}

/*************************************
 * Retrieve Transition Charges 
 * Overloaded for projection w/ 
 * alpha/beta spins
 * **********************************/
void getCharges(string filename, Molecule *mol, int nstates, int a, int b) {
  ifstream infile;
  infile.open(filename.c_str());
  char tempc[1000];
  string temps;
  
  //skip natoms on first line
  infile.getline(tempc,1000);
  
  for (int j=0; j<mol->natoms; j++) {
    infile.getline(tempc,1000);
    temps = strtok(tempc," ");
    temps = strtok(NULL," ");
    mol->atoms[j].type = temps;

    /** Get atomic positions **/
    temps = strtok(NULL," ");
    mol->atoms[j].pos[0] = atof(temps.c_str());
    mol->atoms[j].ipos[0] = atof(temps.c_str());
    temps = strtok(NULL," ");
    mol->atoms[j].pos[1] = atof(temps.c_str());
    mol->atoms[j].ipos[1] = atof(temps.c_str());
    temps = strtok(NULL," ");
    mol->atoms[j].pos[2] = atof(temps.c_str());
    mol->atoms[j].ipos[2] = atof(temps.c_str());
    temps = strtok(NULL," ");
    
    /** Get the charges **/
    mol->atoms[j].charges[a+b*nstates] = atof(temps.c_str());
    mol->atoms[j].charges[b+a*nstates] = atof(temps.c_str());
  
    mol->atoms[j].spos[0] = sqrt(mol->atoms[j].pos[0]*mol->atoms[j].pos[0]
                                 + mol->atoms[j].pos[1]*mol->atoms[j].pos[1] 
                                 + mol->atoms[j].pos[2]*mol->atoms[j].pos[2]);
    mol->atoms[j].spos[1] = acos(mol->atoms[j].pos[2]/mol->atoms[j].spos[0]);
    mol->atoms[j].spos[2] = atan(mol->atoms[j].pos[1]/mol->atoms[j].pos[0]);
  
  }
}

void getCharges(string filename, Molecule *mol, int nstates, int a, int b, string spin) {
  ifstream infile;
  infile.open(filename.c_str());
  char tempc[1000];
  string temps;
  
  //skip natoms on first line
  infile.getline(tempc,1000);
  
  for (int j=0; j<mol->natoms; j++) {
    infile.getline(tempc,1000);
    temps = strtok(tempc," ");
    temps = strtok(NULL," ");
    mol->atoms[j].type = temps;

    /** Get atomic positions **/
    temps = strtok(NULL," ");
    mol->atoms[j].pos[0] = atof(temps.c_str());
    mol->atoms[j].ipos[0] = atof(temps.c_str());
    temps = strtok(NULL," ");
    mol->atoms[j].pos[1] = atof(temps.c_str());
    mol->atoms[j].ipos[1] = atof(temps.c_str());
    temps = strtok(NULL," ");
    mol->atoms[j].pos[2] = atof(temps.c_str());
    mol->atoms[j].ipos[2] = atof(temps.c_str());
    temps = strtok(NULL," ");
    
    /** Get the charges **/
    /** if a projection calculation, and alpha/beta spins are desired **/
    if (spin.compare(0,5,"alpha",0,5) == 0) {
      temps = strtok(NULL," ");
    } else if (spin.compare(0,4,"beta",0,4) == 0) {
      temps = strtok(NULL," ");
      temps = strtok(NULL," ");
    }
    mol->atoms[j].charges[a+b*nstates] = atof(temps.c_str());
    mol->atoms[j].charges[b+a*nstates] = atof(temps.c_str());
  
    mol->atoms[j].spos[0] = sqrt(mol->atoms[j].pos[0]*mol->atoms[j].pos[0]
                                 + mol->atoms[j].pos[1]*mol->atoms[j].pos[1] 
                                 + mol->atoms[j].pos[2]*mol->atoms[j].pos[2]);
    mol->atoms[j].spos[1] = acos(mol->atoms[j].pos[2]/mol->atoms[j].spos[0]);
    mol->atoms[j].spos[2] = atan(mol->atoms[j].pos[1]/mol->atoms[j].pos[0]);
  }
}


#endif

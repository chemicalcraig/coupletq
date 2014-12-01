#ifndef PARSE_H
#define PARSE_H

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


/*** Read in MO vectors from movec file ***/
/* str             - input, file name to read
 * mol            - input, Molecule object
 * mol->occupation - output, occupation numbers
 * mol->moeigenv   - output, MO eigenvalues
 * mol->mos        - output, MO in orthonormal basis
 */

bool readMOs(string str, Molecule *mol) {
  cout.precision(10);
  char tempc[1000];
  ifstream infile;
  infile.open(str.c_str());

  if (!infile.is_open()) {
    cout<<"Error opening MO file"<<endl;
    return false;
  } else {
    cout<<"Opening file: "<<str<<endl;
  }
  
  /* Skip right to the occupation numbers */
  getnlines(infile,tempc,14,500);
  mol->nocc = 0;
  mol->nuocc = 0;
  for (int i=0; i<mol->nmo; i++) {
    infile>>tempc;
    mol->occupation[i] = (int)atof(tempc);
    if (mol->occupation[i] == 0) {
      mol->nuocc++;
    } else {
      mol->nocc++;
    }
  }
cout<<"HOMO = "<<mol->nocc<<", # virt. orb. = "<<mol->nuocc<<endl;
  /* Read in eigenvalues */
  for (int i=0; i<mol->nmo; i++) {
    infile>>tempc;
    mol->moeigenv[i] = atof(tempc);
  }

  /* Read in MO vectors */
 /* for (int i=0; i<mol->nmo; i++)
    for (int j=0; j<mol->nbasis; j++) {
      infile>>tempc;
      mol->mos[i+j*mol->nmo] = atof(tempc);
    }
*/
  return true;
}

Molecule *parseComfile(ifstream &comfile) {
  
  //parse input file
  char tempc[1000];
  string temps;
  comfile.getline(tempc,1000);
  temps = strtok(tempc,":");
  temps = strtok(NULL,": ");
  
  //Get number of molecules
  int nmol = atoi(temps.c_str());
  Molecule *mol = new Molecule[nmol];
  for (int i=0; i<nmol; i++)
    mol[i].nmol = nmol;

  //Get order of interaction
  //limited to {1,2}
  //1 = FRET
  //2 = 2nd order perturbation theory
  comfile.getline(tempc,1000);
  temps = strtok(tempc,":");
  temps = strtok(NULL,": ");
  for (int i=0; i<nmol; i++)
    mol[i].interaction = atoi(temps.c_str());

  //Number of electronic states to consider
  //on each molecule
  for (int i=0; i<nmol; i++) {
    comfile.getline(tempc,1000);
    temps = strtok(tempc,":");
    temps = strtok(NULL,": ");
    mol[i].nstates = atoi(temps.c_str());
    cout<<"There are "<<mol[i].nstates<<" excited states for molecule "<<i+1<<endl;
    
    //get natoms and charge densities of each molecule, 
    //  each one requires a separate file generated
    //  using either denstq or aodens
    //this should change if inter-excited transitions are to 
    //be included
    for (int icharges=0; icharges<2*mol[i].nstates-1; icharges++) {
      comfile.getline(tempc,1000);
      temps = strtok(tempc,":");
      temps = strtok(NULL,": ");

      //Natoms
      if (icharges == 0) {
        //mol[i].natoms = getNatoms(temps,nmol,mol);
        mol[i].atoms = new Atom[mol[i].natoms];
        cout<<"Molecule "<<i+1<<" has "<<mol[i].natoms<<" atoms"<<endl;
        //Allocate atoms and all of their densities
        //We use lower triangular form for the couplings
        for (int j=0; j<mol[i].natoms; j++) {
          //this should change if inter-excited transitions are to 
          //be included
          mol[i].atoms[j].allocateCharges(mol[i].nstates*mol[i].nstates);
        }
      }
      
      //Get the transition charges.
      //The diagonal components (j==k) are the
      //state densities, while the off diagonal
      //components are the transition densities.
      //Since we are only considering transitions to/from
      //the ground state only the first column become 
      //populated with transition densities
      //We only need the lower triangle since the 

      //getCharges(temps,nmol,&mol[i],icharges,mol[i].nstates);
    }

    //symmetrize charges
    for (int iatom = 0; iatom<mol[i].natoms; iatom++) {
      for (int icharge = 0; icharge<mol[i].nstates; icharge++) {
        for (int icharge2 = icharge; icharge2<mol[i].nstates; icharge2++) {
          mol[i].atoms[iatom].charges[icharge2 + icharge*mol[i].nstates] = 
              mol[i].atoms[iatom].charges[icharge + icharge2*mol[i].nstates];
        }
      }
    }
  }
  
  //get energies and dipoles of each state of each molecule
  //energies and dipoles are gotten from a tddft calculation.
  //Energies are all relative to the ground state, and we 
  //take that to be zero for all molecules.
  
  for (int i=0; i<nmol; i++) {
    comfile.getline(tempc,1000);
    temps = strtok(tempc,":");
    temps = strtok(NULL,": ");
    
    //retrieves nroots and respective energies and
    //transition dipoles for all 
    //excited states in given log file

    //getTDDFT(temps,&mol[i]);

    cout<<"The transition dipole vector for S1 of molecule "<<i+1<<
    " is ("<<mol[i].transmoment[0+3*0]<<" x, "<<mol[i].transmoment[0+3*1]<<" y "<<
    mol[i].transmoment[0+3*2]<<" z) "<<endl;
    cout<<"The ground state energy of molecule "<<i+1<<" is "<<mol[i].groundenergy<<" a.u."<<endl;
  }

  //get output file name to which coupling will be written
  comfile.getline(tempc,1000);
  temps = strtok(tempc,":");
  temps = strtok(NULL,": ");
  mol[0].outputfilename = temps;
  cout<<"Couplings are written to : "<<mol[0].outputfilename<<endl;

  /** Set Molecular Mass **/
  for (int i=0; i<nmol; i++) {
    double mass = 0.;
    for (int j=0; j<mol[i].natoms; j++) {
      mass += mol[i].atoms[j].mass;
    }
    mol[i].setMass(mass);
  }
  return mol;
}

#endif // PARSE_H

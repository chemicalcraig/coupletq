/*************************************************************
 * Programs to calculate the three-body
 * coupling elements and energy transfer
 * rates in a perturbative framework using
 * transition charges.
 *
 * Input: transition charges of each molecule involved,
 *        calculated using either AO's or spatial densities
 * 
 * CTC 10-10-14
 * **********************************************************/

#include "main.h"

int main(int argc, char **argv) {
  //open input file
  comfile.open(argv[1]);

  //parse input file
  char tempc[1000];
  string temps;
  comfile.getline(tempc,1000);
  temps = strtok(tempc,":");
  temps = strtok(NULL,": ");
  
  //Get number of molecules
  nmol = atoi(temps.c_str());
  Molecule *mol = new Molecule[nmol];

  //Number of electronic states to consider
  //on each molecule
  for (int i=0; i<nmol; i++) {
    comfile.getline(tempc,1000);
    temps = strtok(tempc,":");
    temps = strtok(NULL,": ");
    mol[i].nstates = atoi(temps.c_str());
    cout<<"There are "<<mol[i].nstates<<" states for molecule "<<i+1<<endl;
    
    //get natoms and charge densities of each molecule, 
    //  each one requires a separate file generated
    //  using either denstq or aodens with the
    //  appropriate cube file as input
    //this should change if inter-excited transitions are to 
    //be included
    for (int icharges=0; icharges<2*mol[i].nstates-1; icharges++) {
      comfile.getline(tempc,1000);
      temps = strtok(tempc,":");
      temps = strtok(NULL,": ");

      //Natoms
      if (icharges == 0) {
        natoms = getNatoms(temps,nmol,mol);
        mol[i].natoms = natoms;
        mol[i].atoms = new Atom[natoms];
        cout<<"Molecule "<<i+1<<" has "<<natoms<<" atoms"<<endl;
      }

      
      //Allocate atoms and all of their densities
      //We use lower triangular form for the couplings
      for (int j=0; j<natoms; j++) {
        //this should change if inter-excited transitions are to 
        //be included
        mol[i].atoms[j].charges = new double[2*mol[i].nstates-1];
      }
    
      //Get the transition charges.
      //The diagonal components (j==k) are the
      //state densities, while the off diagonal
      //components are the transition densities.
      //Since we are only considering transitions to/from
      //the ground state only the first column become 
      //populated with transition densities
      //We only need the lower triangle since the 
      getCharges(temps,nmol,mol[i],icharges,mol[i].nstates);
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
    getTDDFT(temps,&mol[i]);

    cout<<"The transition dipole vector for S1 of molecule "<<i+1<<
    " is ("<<mol[i].transmoment[0+3*0]<<" x, "<<mol[i].transmoment[0+3*1]<<" y "<<
    mol[i].transmoment[0+3*2]<<" z) "<<endl;
    cout<<"The ground state energy of molecule "<<i+1<<" is "<<mol[i].groundenergy<<" a.u."<<endl;
  }
  return 0;
}

/********************************************************
 * *****************************************************
 * Functions in util.h
 * *****************************************************
 * ******************************************************/



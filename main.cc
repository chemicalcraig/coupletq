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
  
  //Number of molecules
  nmol = atoi(temps.c_str());
  Molecule *mol = new Molecule[nmol];

  //Number of electronic states to consider
  //on each molecule


  for (int i=0; i<nmol; i++) {
    comfile.getline(tempc,1000);
    temps = strtok(tempc,":");
    temps = strtok(NULL,": ");
    mol[i].nstates = atoi(temps.c_str());
    
    //get natoms and charge densities of each molecule, 
    //  each one requires a separate file generated
    //  using either denstq or aodens with the
    //  appropriate cube file as input
    int nupper = 0;
    for (int iupper = 0; iupper < mol[i].nstates; iupper++) nupper += iupper+1;
    
    for (int icharges=0; icharges<nupper; icharges++) {
      comfile.getline(tempc,1000);
      temps = strtok(tempc,":");
      temps = strtok(NULL,": ");

      //Natoms
      int natoms = getNatoms(temps,nmol,mol);
      mol[i].natoms = natoms;
      mol[i].atoms = new Atom[natoms];
    
      //Allocate atoms and all of their densities
      //We use upper triangular form for the couplings
      for (int j=0; j<natoms; j++) {
        mol[i].atoms[j].charges = new double[nstates*nstates];
      }
    
      //get first charges for mol i
      getCharges(temps,nmol,mol[i],0,0);

      //get the rest of the charges for molecule i
      //we skip (0,0) ie the ground state since we
      //just got that above.
      //The diagonal components (j==k) are the
      //state densities, while the off diagonal
      //components are the transition densities.
      //We only need the upper triangle since the 
      //density and transition charges are symmetric

      for (int j=0; j<mol[i].nstates; j++) {
        for (int k=j; k<mol[i].nstates; k++) {
          if (j==0 && k==0) continue;
            comfile.getline(tempc,1000);
            temps = strtok(tempc,":");
            temps = strtok(NULL,": ");
            getCharges(temps,nmol,mol[i],j,k);
        }
      }
    } 
  }
  
  //get energies of each state of each molecule
  //energies are gotten from a tddft calculation, which
  //are all relative to the ground state, and we take that to 
  //be zero for all molecules.
  
  for (int i=0; i<nmol; i++) {
    comfile.getline(tempc,1000);
    temps = strtok(tempc,":");
    temps = strtok(NULL,": ");
    
    //retrieves nroots and respective energies for all 
    //excited states in given log file
    getTDDFT(temps,mol[i]);
  }
  return 0;
}

/********************************************************
 * *****************************************************
 * Functions
 * *****************************************************
 * ******************************************************/

/* Get Natoms */
int getNatoms(string filename, int nmol, Molecule *mol) {
  ifstream infile;
  cout<<filename.c_str()<<endl;
  infile.open(filename.c_str());
  char tempc[1000];
  infile.geline(tempc,1000);
  string temps = strtok(tempc,": ");
  temps = strtok(NULL,": ");

  return atoi(temps.c_str());
}

/* Get the charges */
void getCharges(string filename, int nmol, Molecule mol, int ndens1, int ndens2) {
  ifstream infile;
  infile.open(filename.c_str());
  char tempc[1000];
  string temps;

  while(infile.good()) {
    infile.getline(tempc,1000);
    temps = strtok(tempc," ");
    for (int i=0; i<2; i++) temps = strtok(NULL," ");
    mol.atoms[j].x = atof(temps.c_str());
    temps = strtok(NULL," ");
    mol.atoms[j].y = atof(temps.c_str());
    temps = strtok(NULL," ");
    mol.atoms[j].z = atof(temps.c_str());
    temps = strtok(NULL," ");
    mol.atoms[j].charges[ndens1+nstates*ndens2] = atof(temps.c_str());
    cout<<temps<<endl;
  }
}

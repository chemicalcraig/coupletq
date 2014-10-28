#ifndef UTIL_H
#define UTIL_H

#define ang2au 1.889725989;

#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include "atom.h"
#include "molecule.h"
#include <iomanip>
using namespace std;

/********************************************************
 * Functions
 * *******************************************************/

/** Compute the coupling between two distributions of transition charges **/
/** dcharge and acharge specify which transition to use **/
double computeCoupling(Molecule mold, Molecule mola,int dcharge, int acharge) {
  double res = 0.;
  double interaction = 0.;
  for (int i=0; i<mold.natoms; i++) {
    for (int j=0; j<mola.natoms; j++) {
      double rda2 = 0.;
      for (int k=0; k<3; k++) {
        rda2 += (mold.atoms[i].pos[k] - mola.atoms[j].pos[k])
          *(mold.atoms[i].pos[k] - mola.atoms[j].pos[k]);
      }
      double rda = sqrt(rda2);
      interaction = mold.atoms[i].charges[dcharge] * mola.atoms[j].charges[acharge] / rda;
      res += interaction;
    }
  }
  return res;
}

/**** Skip several lines in input file ****/
void getnlines(ifstream &in, char *temp, int n, int length) {
  for (int i=0; i<n; i++) {
    in.getline(temp,length);
  }
}

/** Calculate molecular dipole moment
 * from transition charges **/
bool calcdip(Molecule *mol) {
  
  double dx = 0.;
  double dy = 0.;
  double dz = 0.;
  
  for (int i=0; i<mol->natoms; i++) {
    dx += mol->atoms[i].x * mol->atoms[i].tq;//indo;
    dy += mol->atoms[i].y * mol->atoms[i].tq;//indo;
    dz += mol->atoms[i].z * mol->atoms[i].tq;//indo;
  }

  dx *= ang2au;
  dy *= ang2au;
  dz *= ang2au;

  mol->dipx = dx;
  mol->dipy = dy;
  mol->dipz = dz;
  
  cout<<"Transition dipole moment: "<<dx<<" x, "<<dy<<" y, "<<dz<<" z"<<endl;
  cout<<"Magnitude = "<<sqrt(dx*dx+dy*dy+dz*dz)/0.39345<<endl;
  return true;
}

/*** Calculate the transition charges for a given root ***/
bool calctqindo(Molecule *mol, int root) {
 double sum = 0.; 
  for (int atom=0; atom<mol->natoms; atom++) {

    mol->atoms[atom].tqindo = 0.;

    for (int b=0; b<mol->nbasis; b++) { //basis functions
      if (mol->nbasisatom[b] != atom+1) continue;
      for (int i=0; i<mol->nocc; i++) { //MO's ci's
        for (int j=mol->nocc; j<mol->nmo; j++) { //MO ci's
              mol->atoms[atom].tqindo += mol->ci[root*mol->nmo*mol->nmo+i+j*mol->nmo]
                  * mol->mos[i+b*mol->nmo] * mol->mos[j+b*mol->nmo];
        }//end unoccupied MO
      }//end occupied MO
    } //end NAO on atom of interest
    cout<<"atom "<<atom<<endl;
    mol->atoms[atom].tqindo *= sqrt(2);
    cout<<mol->nbasisatom[atom]<<" "<<mol->atoms[atom].type<<" "<<mol->atoms[atom].tqindo<<" "<<endl;
    sum += mol->atoms[atom].tqindo;
  } //end atoms
  
  cout<<"Net charge: "<<sum<<endl;
  return true;
}



/*** Calculate the transition charges for a given root ***/
bool calctq(Molecule *mol, int root) {
 double sum = 0.; 
  for (int atom=0; atom<mol->natoms; atom++) {
      
    mol->atoms[atom].tq = 0.;

    for (int b=0; b<mol->nbasis; b++) {
      if (mol->nbasisatom[b] != atom+1) continue;
      for (int c=0; c<mol->nbasis; c++) {
        //if (mol->nbasisatom[c] == atom+1) continue;
        for (int i=0; i<mol->nocc; i++) {
          for (int j=mol->nocc; j<mol->nbasis; j++) {
              mol->atoms[atom].tq += mol->ci[i+j*mol->nmo+root*mol->nmo*mol->nmo]
                    * mol->mos[i+b*mol->nmo] * mol->mos[j+c*mol->nmo]
                    * mol->overlapm[c+b*mol->nbasis];
          }//end unoccupied MO
        }//end occupied MO
      }//end all other AOs
    
    } //end NAO on atom of interest
    mol->atoms[atom].tq *= sqrt(2);
    cout<<mol->atoms[atom].type<<" "<<mol->atoms[atom].tq<<" "<<endl;
    sum += mol->atoms[atom].tq;
  } //end atoms
  cout<<"Net charge: "<<sum<<endl;
  return true;
}

/** Print some stuff **/
void printStuff(Molecule *mol) {
  ofstream outfile;
  outfile.open("aodata.dat");
  
  /* print atomic numbers, elements, positions, and tq's */
  /* natoms first */
  for (int i=0; i<mol->natoms; i++) {
    outfile<<i+1<<" "<<mol->atoms[i].type.c_str()<<" "
      <<mol->atoms[i].x<<" "<<mol->atoms[i].y<<" "<<mol->atoms[i].z<<" "
      <<mol->atoms[i].tq<<endl;
  }
}

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

/* Get the charges */
void getCharges(string filename, int nmol, Molecule *mol, int icharge, int nstates) {
  ifstream infile;
  infile.open(filename.c_str());
  char tempc[1000];
  string temps;
  
  //skip natoms and energy on first two lines
  getnlines(infile,tempc,2,1000);
  
  for (int j=0; j<mol->natoms; j++) {
    infile.getline(tempc,1000);
    temps = strtok(tempc," ");
    for (int i=0; i<2; i++) temps = strtok(NULL," ");
    mol->atoms[j].x = atof(temps.c_str());
    mol->atoms[j].pos[0] = atof(temps.c_str());
    mol->atoms[j].ipos[0] = atof(temps.c_str());
    temps = strtok(NULL," ");
    mol->atoms[j].y = atof(temps.c_str());
    mol->atoms[j].pos[1] = atof(temps.c_str());
    mol->atoms[j].ipos[1] = atof(temps.c_str());
    temps = strtok(NULL," ");
    mol->atoms[j].z = atof(temps.c_str());
    mol->atoms[j].pos[2] = atof(temps.c_str());
    mol->atoms[j].ipos[2] = atof(temps.c_str());
    temps = strtok(NULL," ");
    //this should change if inter-excited transitions are to 
    //be included
    mol->atoms[j].charges[icharge] = atof(temps.c_str());
  }
}

int findCol(int x, const int rowrange, const int colrange) {
  for (int i=0; i<rowrange; i++) 
    for (int j=0; j<colrange; j++)
      if ((i+j*rowrange) == x)
        return j;
}

int findRow(int x, const int rowrange, const int colrange) {
  for (int i=0; i<rowrange; i++) 
    for (int j=0; j<colrange; j++)
      if ((i+j*rowrange) == x)
        return i;
}

#endif // UTIL_H

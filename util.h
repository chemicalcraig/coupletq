#ifndef UTIL_H
#define UTIL_H

#define ang2au 1.889725989;

#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include "atom.h"
#include "molecule.h"
#include "coulomb.h"
#include <iomanip>
#include "gsl/gsl_blas.h"
using namespace std;

/** Kronecker Delta **/
#define Kronecker(a,b) ((a == b ? 1 : 0))

/********************************************************
 * Functions
 * *******************************************************/
/** window function **/
double window(double center, double value, double delta, int which) {
  switch (which) {
    case 0: //square well
      double sval = value-center;
      sval *= sval;
      if (sval <= delta*delta)
        return 1;
      else
        return 0;
      break;
  }
}

/** Calculate Coulomb coupling between chromophores I and J
 * undergoing transisions between i and j, and k and l, respectively **/
double getCoulomb(Molecule *mol, int I, int J, int i, int j, int k, int l) {
  double res, sum, temp, pos[3];
  temp = 0.;
  for (int ii=0; ii<mol[I].natoms; ii++) {
    for (int jj=0; jj<mol[J].natoms; jj++) {
      double r12 = 0.;
      for (int kk=0; kk<3; kk++) {
        pos[kk] = mol[I].atoms[ii].pos[kk] - mol[J].atoms[jj].pos[kk];
      }
      r12 = cblas_ddot(3,pos,1,pos,1);
      r12 = sqrt(r12);

      temp += mol[I].atoms[ii].charges[i+j*mol[I].nstates] 
            * mol[J].atoms[jj].charges[k+l*mol[J].nstates]/r12;
     // if (I==0 && J==1 && i==0 && j==1 && k==1 && l==0){
     //   cout<<"charges "<<temp<<" temp -< "<<mol[I].atoms[ii].charges[i+j*mol[I].nstates]<<" "
     //             <<mol[J].atoms[jj].charges[k+l*mol[J].nstates]<<endl;
     // }

    }
  }
  return temp;
}

void createCoulomb3(Molecule *mol, Coulomb coul) {
  for (int i=0; i<mol[0].nstates; i++) {
    for (int j=0; j<mol[1].nstates; j++) {
      for (int k=0; k<mol[2].nstates; k++) {
        for (int l=0; l<mol[0].nstates; l++) {
          for (int m=0; m<mol[1].nstates; m++) {
            for (int n=0; n<mol[2].nstates; n++) {
              int index = i + j*2 + k*4 + l*8 + m*16 + n * 32;
              int index3 = l+m*2+n*4;
              int index2 = i+j*2+k*4;
              coul.int3[index2+index3*8] = getCoulomb(mol,0,1,i,l,j,m) * Kronecker(k,n)
                                          + getCoulomb(mol,0,2,i,l,k,n) * Kronecker(j,m)
                                          + getCoulomb(mol,1,2,j,m,k,n) * Kronecker(i,l);
   //             if (i==1 && j==0 && l==0 && m == 1 && k==n)
   //               cout<<"temp get "<<getCoulomb(mol,0,1,i,l,j,m)<<endl;

 cout<<"making coulomb "<<index2<<" "<<index3<<" "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<n<<" "
        <<coul.int3[index2+index3*8]<<" "<<getCoulomb(mol,0,1,i,l,j,m)*Kronecker(k,n)<<" "
        <<getCoulomb(mol,0,2,i,l,k,n) * Kronecker(j,m)<<" " 
        <<getCoulomb(mol,1,2,j,m,k,n) * Kronecker(i,l)<<endl;
                //cout<<index2<<" "<<index3<<" "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<n<<" "<<index<<" index "<<endl;
            }
          }
        }
      }
    }
  }
}
/** Wrap DGEMM **/
void multmm(double *one,double *two,double *three,int m)
{
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
				m,m,m,1.0,one,
	            m,two,m,0,three,m);
};
 
/** Create Coulomb Matrix for three chormophores**/
void createCoulomb3(Molecule *mol, double *mat) {
  double sum;

  for (int i=0; i<mol[0].nstates; i++) {
    for (int j=0; j<mol[1].nstates; j++) {
      for (int k=0; k<mol[2].nstates; k++) {
        for (int l=0; l<mol[0].nstates; l++) {
          for (int m=0; m<mol[1].nstates; m++) {
            for (int n=0; n<mol[2].nstates; n++) {
              int index = i + j*2 + k*4 + l*8 + m*16 + n * 32;
              int index3 = l+m*2+n*4;
              int index2 = i+j*2+k*4;
              mat[index2+index3*8] = getCoulomb(mol,0,1,i,l,j,m) * Kronecker(k,n)
                + getCoulomb(mol,0,2,i,l,k,n) * Kronecker(j,m)
                + getCoulomb(mol,1,2,j,m,k,n) * Kronecker(i,l);
            }
          }
        }
      }
    }
  }
}


/** Arrange molecules to desired configuration **/
void arrangeMol(Molecule *mol) {
  for (int i=0; i<mol[0].nmol; i++) {
    mol[i].arrangeMol();
  }
}

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
      //cout<<i<<" "<<j<<" "<<rda<<" r12 "<<endl;
      //cout<<"charges "<<mold.atoms[i].charges[dcharge]<<" "<<mola.atoms[j].charges[acharge]<<endl;
      interaction = mold.atoms[i].charges[dcharge] * mola.atoms[j].charges[acharge] / rda;
      res += interaction;
    }
  }
  return res;
}

/** Orientation factor between two molecules **/
double orient(Molecule mold, Molecule mola,double *r12) {
  
  /** get unit vector in direction of each transition dipole vector **/
  double dipunit1[3],dipunit2[3],temp[3];
  
  for (int i=0; i<3; i++) {
    temp[i] = mold.idip[i] + mold.com[i];
  }

  for (int i=0; i<3; i++) {
    dipunit1[i] = mold.dip[i] / sqrt(mold.dipmag);
    dipunit2[i] = mola.dip[i] / sqrt(mola.dipmag);
  }

  double res = cblas_ddot(3,dipunit1,1,dipunit2,1) 
              -3* cblas_ddot(3,dipunit1,1,r12,1)*cblas_ddot(3,dipunit2,1,r12,1);
  return res;
}


/** PDA coupling **/
double pdaCalc(Molecule *mol, double &res) {
  
  //Get unit vector connecting the com's
  double diff[3];
  double sum=0.;
  
  //construct vector connecting COM's
  for (int i=0; i<3; i++) {
    diff[i] = mol[0].com[i] - mol[1].com[i];
    sum += diff[i]*diff[i];
  }

  //normalize it
  for (int i=0; i<3; i++) {
    diff[i] /= sqrt(sum);
  }

  /** Calculate orientation factor **/
  double kappa = orient(mol[0],mol[1],diff);

  res = kappa * mol[1].dipmag*mol[1].dipmag / (sqrt(sum)*sum);
  //cout<<kappa<<" "<<mol[1].dipmag<<" "<<sqrt(sum)<<endl;
}

/**** Skip several lines in input file ****/
void getnlines(ifstream &in, char *temp, int n, int length) {
  for (int i=0; i<n; i++) {
    in.getline(temp,length);
  }
}

/** Calculate molecular dipole moment
 * from transition charges **/
bool calcdip(Molecule &mol) {
  
  double dx = 0.;
  double dy = 0.;
  double dz = 0.;
  
  for (int i=0; i<mol.natoms; i++) {
    dx += mol.atoms[i].pos[0] * mol.atoms[i].charges[1];
    dy += mol.atoms[i].pos[1] * mol.atoms[i].charges[1];
    dz += mol.atoms[i].pos[2] * mol.atoms[i].charges[1];
  }

  dx *= ang2au;
  dy *= ang2au;
  dz *= ang2au;

  mol.dip[0] = dx;
  mol.dip[1] = dy;
  mol.dip[2] = dz;
  mol.dipmag = sqrt(dx*dx+dy*dy+dz*dz)/0.39345;
  
  cout<<"Transition dipole moment: "<<dx<<" x, "<<dy<<" y, "<<dz<<" z"<<endl;
  cout<<"Magnitude = "<<sqrt(dx*dx+dy*dy+dz*dz)<<endl;
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
    temps = strtok(NULL," ");
    mol->atoms[j].type = temps;
    temps = strtok(NULL," ");
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
    cout<<"Atom pos "<<j<<" from "<<filename<<" "<<mol->atoms[j].pos[0]<<" "<<mol->atoms[j].pos[1]<<" "<<mol->atoms[j].pos[2]<<endl;
    //this should change if inter-excited transitions are to 
    //be included
    if (icharge < mol->nstates) { //diagonal terms (state densities)
      mol->atoms[j].charges[icharge+icharge*mol->nstates] = atof(temps.c_str());

    } else {
      mol->atoms[j].charges[icharge] = atof(temps.c_str());
    }
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

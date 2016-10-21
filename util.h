#ifndef UTIL_H
#define UTIL_H

#define ang2au 1.889725989;

#include <cstring>
#include <cmath>
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

/** Print level **/
bool DEBUG = false;

/********************************************************
 * Functions
 * *******************************************************/

/** Calculates the distance-independent Coulomb interaction
 * for projection calculations
 */
double getCoulombNoDist(Molecule *mol, int I, int J, int i, int j, int k, int l) {
  double res, sum, temp, pos[3];
  temp = 0.;
  for (int ii=0; ii<mol[I].natoms; ii++) {
    for (int jj=0; jj<mol[J].natoms; jj++) {
      temp += mol[I].atoms[ii].charges[i+j*mol[I].nstates] 
            * mol[J].atoms[jj].charges[k+l*mol[J].nstates];
    }
  }

  return temp;
}

/**** Scale transition charges for one electron transfer ****/
void scaletq(Molecule *mol) {
  cout<<"Scaling transition charges"<<endl;
  double sump1 = 0.;
  double summ1 = 0.;
  double sump2 = 0.;
  double summ2 = 0.;

  for (int i=0; i<mol[0].natoms; i++) {
    if (mol[1].atoms[i].charges[0+mol[1].fstate*mol[1].nstates] < 0) {
      summ1 += mol[1].atoms[i].charges[0+mol[1].fstate*mol[1].nstates];
    }
    else if (mol[1].atoms[i].charges[0+mol[1].fstate*mol[1].nstates] > 0)
      sump1 += mol[1].atoms[i].charges[0+mol[1].fstate*mol[1].nstates];
    if (mol[0].atoms[i].charges[0+mol[0].fstate*mol[0].nstates] < 0)
      summ2 += mol[0].atoms[i].charges[0+mol[0].fstate*mol[0].nstates];
    else if (mol[0].atoms[i].charges[0+mol[0].fstate*mol[0].nstates] > 0)
      sump2 += mol[0].atoms[i].charges[0+mol[0].fstate*mol[0].nstates];

  }

  cout<<"charges "<<summ1<<" "<<sump1<<" "<<summ2<<" "<<sump2<<endl;
  double delta = 1.0e-5;
  
  if (((summ1+sump1) >= delta) || ((summ2+sump2) >= delta)) {
    cout<<"positive and negative charges unequal, exiting"<<endl;
    exit(0);
  }

/** Check sign for scaling to make sure electron/holes match for diff mol **/
  for (int i=0; i<mol[0].natoms; i++) {
    mol[0].atoms[i].charges[0+mol[0].fstate*mol[0].nstates] /= sqrt(2)/sump2;
    mol[1].atoms[i].charges[0+mol[1].fstate*mol[1].nstates] /= -sqrt(2)/sump1;
  }

/*
 * Checking sum
 *
  sump1=0.;
  summ1=0.;
  sump2=0.;
  summ2=0.;
  for (int i=0; i<mol[0].natoms; i++) {
    if (mol[1].atoms[i].charges[0+mol[1].fstate*mol[1].nstates] < 0) {
      summ1 += mol[1].atoms[i].charges[0+mol[1].fstate*mol[1].nstates];
    }
    else if (mol[1].atoms[i].charges[0+mol[1].fstate*mol[1].nstates] > 0)
      sump1 += mol[1].atoms[i].charges[0+mol[1].fstate*mol[1].nstates];
    if (mol[0].atoms[i].charges[0+mol[0].fstate*mol[0].nstates] < 0)
      summ2 += mol[0].atoms[i].charges[0+mol[0].fstate*mol[0].nstates];
    else if (mol[0].atoms[i].charges[0+mol[0].fstate*mol[0].nstates] > 0)
      sump2 += mol[0].atoms[i].charges[0+mol[0].fstate*mol[0].nstates];

  }

  cout<<"charges "<<summ1<<" "<<sump1<<" "<<summ2<<" "<<sump2<<endl;
*/
}



/** Calculates the absolute transition charge difference between
 * different excited states */

double *projecttq(Molecule *mol) {
  /** atom by atom **/
  double sum = 0;
  double *temp = new double[mol[0].natoms];
  for (int i=0; i<mol[0].natoms; i++) {
    double proj;

      proj = abs((mol[1].atoms[i].charges[0+mol[1].fstate*mol[1].nstates]
                      - mol[0].atoms[i].charges[0+mol[0].fstate*mol[0].nstates])
                   )  ;//  / mol[0].atoms[i].charges[0+mol[0].fstate*mol[0].nstates]);
    temp[i] = proj;
    sum += proj;
  }
  cout<<"projection = "<<1.-sum<<endl;
  return temp;
}

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

      //CTC scaling
      //r12 /= 3.375;
      
      //convert from angstroms to a0
      r12 /= 0.529177249;
      
      temp += mol[I].atoms[ii].charges[i+j*mol[I].nstates] 
            * mol[J].atoms[jj].charges[k+l*mol[J].nstates]
            /r12;
    }
  }

  /** print stuff **/
  if (DEBUG) {
   cout<<"coulomb for <"<<i<<k<<"|V_"<<I<<J<<"|"<<j<<l<<"> = "<<temp<<endl;
  }
  
  return temp;
}

/** Create coulomb matrix, overloaded for windowing **/
void createCoulomb3(Molecule *mol, Coulomb coul, double *en,Reader r) {
  for (int i=0; i<mol[0].nindices; i++) {
    for (int j=0; j<mol[0].nindices; j++) {
      coul.int3[i+j*mol[0].nindices] = 0.;
      for (int m1=0; m1<mol[0].nmol; m1++) {
        for (int m2=0; m2<mol[0].nmol; m2++) {
          if (m1 == m2) continue;
          int kron=1;
          for (int k=0; k<mol[0].nmol; k++) {
            if ((k==m1) || (k==m2)) continue;
            kron *= Kronecker(mol[0].indices[k+i*mol[0].nmol],mol[0].indices[k+j*mol[0].nmol]);
          }
          coul.int3[i+j*mol[0].nindices] += getCoulomb(mol,m1,m2,
                               mol[0].indices[m1+i*mol[0].nmol],
                               mol[0].indices[m1+j*mol[0].nmol],
                               mol[0].indices[m2+i*mol[0].nmol],
                               mol[0].indices[m2+j*mol[0].nmol])*kron
                                *window(en[i],en[j],r.calc.ewindow,0);
        }//end molecule 2
      }//end molecule 1
      //fix over counting
      coul.int3[i+j*mol[0].nindices] /= 2.;
      if (DEBUG) {
        cout<<"making coulomb with window "<<i<<" "<<j<<" "<<coul.int3[i+j*mol[0].nindices]<<endl;
      }
    }//end index column
  }//end index row
}

void createCoulomb3(Molecule *mol, Coulomb coul) {
  cout<<92<<" "<<mol[0].nindices<<endl;
  for (int i=0; i<mol[0].nindices; i++) {
    for (int j=0; j<mol[0].nindices; j++) {
      coul.int3[i+j*mol[0].nindices] = 0.;
      for (int m1=0; m1<mol[0].nmol; m1++) {
        for (int m2=m1; m2<mol[0].nmol; m2++) {
          if (m1 == m2) continue;
          int kron=1;
          for (int k=0; k<mol[0].nmol; k++) {
            if ((k==m1) || (k==m2)) continue;
            kron *= Kronecker(mol[0].indices[k+i*mol[0].nmol],mol[0].indices[k+j*mol[0].nmol]);
          }
          double dum = getCoulomb(mol,m1,m2,
                               mol[0].indices[m1+i*mol[0].nmol],
                               mol[0].indices[m1+j*mol[0].nmol],
                               mol[0].indices[m2+i*mol[0].nmol],
                               mol[0].indices[m2+j*mol[0].nmol]);

          coul.int3[i+j*mol[0].nindices] += dum*kron;
        }//end molecule 2
      }//end molecule 1
      //adjust for over counting
      //coul.int3[i+j*mol[0].nindices] /= 2.;
      if (DEBUG) {
        cout<<"making coulomb "<<i<<" "<<j<<" "<<coul.int3[i+j*mol[0].nindices]<<endl;
      }
    }//end index column
  }//end index row  
}


/** Wrap DGEMM **/
void multmm(double *one,double *two,double *three,int m)
{
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
				m,m,m,1.0,one,
	            m,two,m,0,three,m);
};
 

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
      
      //CTC scaling
      //rda /= 3.375;
      
      //convert from angstroms to a0
      rda *= 1.889725988579;
      
      interaction = //85 * 
                    mold.atoms[i].charges[0+mold.fstate*mold.nstates]
                      * mola.atoms[j].charges[0+mola.fstate*mola.nstates] 
                      / (rda );
      //interaction = mold.atoms[i].charges[dcharge] * mola.atoms[j].charges[acharge] / rda;
      res += interaction;
    }
  }
  return res;
}

/** Orientation factor between two molecules **/
double orient(Molecule mold, Molecule mola,double *r12,double dist) {
  
  /** convert distance to au **/
  dist *= ang2au;

  /** get unit vector in direction of each transition dipole vector **/
  double dipunit1[4],dipunit2[4],temp[3],temp2[3];
  double sum=0;
  double sum2=0;

  /** temp/temp2 are the transition dipoles **/
  for (int i=0; i<3; i++) {
    temp[i] = mold.transmoment[i*3+(mold.target-1)*9];
    temp2[i] = mola.transmoment[i*3+(mola.target-1)*9] - mola.com[i];
    
    /** convert r12 to au **/
    r12[i] *= ang2au;
  }

  /** sum & sum2 are the transition dipole magnitudes **/
  for (int i=0; i<3; i++) {
    sum += temp[i]*temp[i];
    sum2 += temp2[i]*temp2[i];
  }
  //sum = sqrt(sum);
  //sum2 = sqrt(sum2);

  for (int i=0; i<3; i++) {
    dipunit1[i+1] = temp[i];
    dipunit2[i+1] = temp2[i];
  }

  double res = (cblas_ddot(3,dipunit1,1,dipunit2,1)/(dist*dist*dist))
              -3* cblas_ddot(3,dipunit1,1,r12,1)*cblas_ddot(3,dipunit2,1,r12,1)/(dist*dist*dist*dist*dist);
  
  return res;
}


/** PDA coupling **/
double pdaCalc(Molecule *mol, double &res) {
  
  //Get unit vector connecting the com's
  double diff[3];
  double sum=0.;
  
  //construct vector connecting COM's
  for (int i=0; i<3; i++) {
    diff[i] = mol[1].com[i] - mol[0].com[i];
    sum += diff[i]*diff[i];
  }
  
  //normalize it
  //for (int i=0; i<3; i++) {
  //  diff[i] = sqrt(diff[i]*diff[i]);
  //}

  /** Calculate orientation factor **/
  double kappa = orient(mol[1],mol[0],diff,sqrt(sum));
  
  double temp[3],temp2[3];
  double sum3,sum4;
  sum3=0.;
  sum4=0.;
  /** temp/temp2 are the transition dipoles translated to the COM **/
  for (int i=0; i<3; i++) {
    temp[i] = mol[0].transmoment[i*3+(mol[0].target-1)*9];
    temp2[i] = mol[1].transmoment[i*3+(mol[1].target-1)*9];
  }

  /** sum & sum2 are the transition dipole magnitudes **/
  for (int i=0; i<3; i++) {
    sum3 += temp[i]*temp[i];
    sum4 += temp2[i]*temp2[i];
  }
  sum3 = sqrt(sum3);
  sum4 = sqrt(sum4);

  res = kappa ;// / (sqrt(sum)*sum);
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
    dx += mol.atoms[i].pos[0] * mol.atoms[i].charges[0];
    dy += mol.atoms[i].pos[1] * mol.atoms[i].charges[0];
    dz += mol.atoms[i].pos[2] * mol.atoms[i].charges[0];
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

/*
int *indicesFromIndex(const int ind, Molecule *mol) {
  int prod = 0;
  int *a;
  for (int i=0; i<mol[0].nmol; i++) {
    prod *= mol[i].nstates;
  }
  a=new int[prod];
  for (int i=0; i<prod; i++) {
       
  }
}
*/

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

/** indices from index, hard-coded for up to 7 molecules for now **/
int *indicesFromIndex(const int ind, Molecule *mol) {
  int *a;
  a=new int[mol[0].nmol];
  
  switch(mol[0].nmol) {
    case 2:
      for (int i=0; i<mol[0].nstates; i++)
        for (int j=0; j<mol[1].nstates; j++) {
          int x = i+j*mol[0].nstates;
          if (x == ind) {
            a[0] = i;
            a[1] = j;
            return a;
          }
        }
      break;
    case 3:
      for (int i=0; i<mol[0].nstates; i++)
        for (int j=0; j<mol[1].nstates; j++) 
          for (int k=0; k<mol[2].nstates; k++) {
            int x = i+j*mol[0].nstates+k*mol[0].nstates*mol[1].nstates;
            if (x == ind) {
              a[0] = i;
              a[1] = j;
              a[2] = k;
              return a;
            }
          }
      break;
    case 4:
      for (int i=0; i<mol[0].nstates; i++)
        for (int j=0; j<mol[1].nstates; j++) 
          for (int k=0; k<mol[2].nstates; k++) 
            for (int l=0; l<mol[3].nstates; l++) {
              int x = i + j*mol[0].nstates + k*mol[0].nstates*mol[1].nstates + l*mol[0].nstates*mol[1].nstates*mol[2].nstates;
              if (x == ind) {
                a[0] = i;
                a[1] = j;
                a[2] = k;
                a[3] = l;
                return a;
              }
          }
      break;
    case 5:
      for (int i=0; i<mol[0].nstates; i++)
        for (int j=0; j<mol[1].nstates; j++) 
          for (int k=0; k<mol[2].nstates; k++) 
            for (int l=0; l<mol[3].nstates; l++) 
              for (int m=0; m<mol[4].nstates; m++) {
                int x = i + j*mol[0].nstates + k*mol[0].nstates*mol[1].nstates 
                    + l*mol[0].nstates*mol[1].nstates*mol[2].nstates
                    + m*mol[0].nstates*mol[1].nstates*mol[2].nstates*mol[3].nstates;
                if (x == ind) {
                  a[0] = i;
                  a[1] = j;
                  a[2] = k;
                  a[3] = l;
                  a[4] = m;
                  return a;
                }
              }
      break;
    case 6:
      for (int i=0; i<mol[0].nstates; i++)
        for (int j=0; j<mol[1].nstates; j++) 
          for (int k=0; k<mol[2].nstates; k++) 
            for (int l=0; l<mol[3].nstates; l++) 
              for (int m=0; m<mol[4].nstates; m++) 
                for (int n=0; n<mol[5].nstates; n++) {
                  int x = i + j*mol[0].nstates + k*mol[0].nstates*mol[1].nstates 
                      + l*mol[0].nstates*mol[1].nstates*mol[2].nstates
                      + m*mol[0].nstates*mol[1].nstates*mol[2].nstates*mol[3].nstates
                      + n*mol[0].nstates*mol[1].nstates*mol[2].nstates*mol[3].nstates*mol[4].nstates;
                  if (x == ind) {
                    a[0] = i;
                    a[1] = j;
                    a[2] = k;
                    a[3] = l;
                    a[4] = m;
                    a[5] = n;
                    return a;
                }
              }
      break;
    case 7:
      for (int i=0; i<mol[0].nstates; i++)
        for (int j=0; j<mol[1].nstates; j++) 
          for (int k=0; k<mol[2].nstates; k++) 
            for (int l=0; l<mol[3].nstates; l++) 
              for (int m=0; m<mol[4].nstates; m++) 
                for (int n=0; n<mol[5].nstates; n++) 
                  for (int o=0; o<mol[6].nstates; o++) {
                  int x = i + j*mol[0].nstates + k*mol[0].nstates*mol[1].nstates 
                      + l*mol[0].nstates*mol[1].nstates*mol[2].nstates
                      + m*mol[0].nstates*mol[1].nstates*mol[2].nstates*mol[3].nstates
                      + n*mol[0].nstates*mol[1].nstates*mol[2].nstates*mol[3].nstates*mol[4].nstates
                      
                      + o*mol[0].nstates*mol[1].nstates*mol[2].nstates*mol[3].nstates
                            *mol[4].nstates*mol[5].nstates;
                  if (x == ind) {
                    a[0] = i;
                    a[1] = j;
                    a[2] = k;
                    a[3] = l;
                    a[4] = m;
                    a[5] = n;
                    a[6] = o;
                    return a;
                }
              }
      break;
  }
}
/****************************************
 * Initialize the index matrix
 * **************************************/
void setIndices(Molecule *mol, const int x, const int y) {
  mol[0].indices = new int[x*y];
  int *a;
  for (int i=0; i<y; i++) {
    a = indicesFromIndex(i,mol);
    for (int j=0; j<x; j++) {
      mol[0].indices[j+i*x] = a[j];
    }
  }
}


#endif // UTIL_H

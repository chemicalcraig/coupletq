#ifndef PERT_H
#define PERT_H

/**************************************
 * This file contains routines
 * to perform a perturbative 3-body calculation.
 * ***********************************/

#include <cstring>
#include <fstream>
#include <iostream>
#include "atom.h"
#include "molecule.h"
#include "coulomb.h"
#include <iomanip>
using namespace std;

/** Perturbation Calculation in eigenbasis of the 3-bit 
 * Coulomb operator **/
void pertCalc(Molecule *mol, Coulomb coul, double *energies,double *int3,double dum) {
  cout<<" *** in pertCalc *** "<<endl;
  for (int i=0; i<8; i++) coul.int3[i+i*8] = 0.;
  for (int i=0; i<8; i++) { //initial state
    
    for (int j=0; j<8; j++) { //final state
cout<<i<<" "<<j<<" "<<coul.int3[i+j*8]<<" "<<int3[i+j*8]<<" "<<coul.evecs3[i+j*8]<<endl;
      coul.int3[i+j*8] *= 10.;

    }
  }
  double res[64];
  //Make diagonal hamiltonian with exciton energies on the diagonal
  double ham[64],tildeint[64];
  for (int i=0; i<8; i++) {
    ham[i+i*8] = energies[i];
  }
  //Make C.ham.C^T
  double ham2[64],tempm[64],tempm2[64];
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,8,8,8,1.,
              coul.evecs3,8,ham,8,0,tempm,8);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,8,8,8,1.,
              tempm,8,coul.evecs3,8,0,ham2,8);
  //Make C.V.C^T
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,8,8,8,1.,
              coul.int3,8,coul.int3,8,0,tempm,8);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,8,8,8,1.,
              tempm,8,coul.evecs3,8,0,tildeint,8);
  
  double numerator, denominator;
  for (int i=0; i<8; i++) { //initial state
    for (int j=0; j<8; j++) { //final state
      numerator = 0.;
      denominator = 0.;
    /*  for (int k=0; k<8; k++) { //intermediate state
        double sum1 = 0.;
        double sum2 = 0.;
        double sum3 = 0.;
        double sum4 = 0.;
        if (i==k) continue;
        for (int a=0; a<8; a++) {
          for (int b=0; b<8; b++) {
            if (a!=b) continue;
            sum1 += coul.evecs3[i+a*8]*coul.evecs3[b+k*8]*coul.int3[a+b*8];//[b+a*8];
          }
        }
        for (int g=0; g<8; g++) {
          for (int d=0; d<8; d++) {
            if (g!=d) continue;
            sum2 += coul.evecs3[k+g*8]*coul.evecs3[d+j*8]*coul.int3[g+d*8];//tildeint[g+d*8];     
          }
        }
        numerator += sum1*sum2;
        for (int a=0; a<8; a++) {
          for (int b=0; b<8; b++) {
            sum3 += coul.evecs3[i+a*8]*coul.evecs3[b+i*8]*ham2[a+b*8];
          }
        }
        for (int a=0; a<8; a++) {
          for (int b=0; b<8; b++) {
            sum4 += coul.evecs3[k+a*8]*coul.evecs3[b+k*8]*ham2[a+b*8];
          }
        }
        denominator += sum3-sum4;//coul.evals3[i] - coul.evals3[k];//sum3-sum4;
      }
      */
      res[i+j*8] = tempm[i+j*8];//numerator/(denominator );
    }
  }
  cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,8,8,8,1.,
              coul.evecs3,8,res,8,0,tempm,8);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,8,8,8,1.,
              tempm,8,coul.evecs3,8,0,tildeint,8);
 for (int i=0; i<8; i++)
   for (int j=0; j<8; j++) cout<<i<<" "<<j<<" "<<res[i+j*8]<<" "<<tildeint[i+j*8]<<" "<<coul.int3[i+j*8]<<endl;
}

void pertCalc(Molecule *mol, Coulomb coul,double *intham,double *energies) {

  double inv, temp,temp2,temp3,r12,r13,r23,sum,sum2,pos[3];
  temp2 = 0.;
  temp = 0.;
  sum = 0.;
  sum2 = 0.;
  temp3 = 0.;

  for (int i=0; i<8; i++) { //initial
    for (int j=0; j<8; j++) {//final
      temp2 = 0.;
      temp3 = 0.;
      for (int k=0; k<8; k++) {
        temp2 += energies[i]*coul.evecs3[i+k*8]*coul.evecs3[i+k*8];
        temp3 += energies[j]*coul.evecs3[j+k*8]*coul.evecs3[j+k*8];
      }
      sum = 0.;
      sum2 = 0.;
      double energy;
      //if (i != j) {
        for (int k=0; k<8; k++) {//intermediate zero order basis
          if (i==k ) continue;
            for (int l=0; l<8; l++) {
              //i term
              energy = (coul.evals3[i]-coul.evals3[k]);//*(coul.evals3[j]-coul.evals3[k]);
              temp = coul.evecs3[i+l*8]*coul.evecs3[k+l*8];
              sum += temp*coul.evals3[l];//energy;
              
              //j term
                     //if ((i == 7 && j==7)) cout<<i<<" "<<j<<" "<<k<<" "<<l<<" sum adding "<<sum<<" "<<sum2<<" "<<temp<<endl;
            }
            for (int l=0; l<8; l++) {
              temp = coul.evecs3[j+l*8]*coul.evecs3[k+l*8];
              sum2 += temp*coul.evals3[l];
 
            }
          }
     // } else {
      /*  for (int k=0; k<8; k++) {
          if (i==k) continue;
          energy = (coul.evals3[k]-coul.evals3[i])*(coul.evals3[k]-coul.evals3[i]);
          sum += coul.evecs3[i+k*8]*coul.evecs3[i+k*8]
            *coul.evecs3[i+k*8]*coul.evecs3[i+k*8]
            *coul.evals3[k]*coul.evals3[k]/energy;
        }
        sum *= -0.5;*/
      //}
      
      intham[i+j*8] = 2*sum*sum2*window(energies[i],energies[j],1,0);
          }
  }
  double mat[64],mat2[64];
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,8,8,8,1.,
              coul.evecs3,8,intham,8,0,mat,8);
      cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,8,8,8,1.,
              mat,8,coul.evecs3,8,0,mat2,8);
      for (int i=0; i<8; i++) 
        for (int j=0; j<8; j++)
        cout<<i<<" "<<j<<" "<<intham[i+j*8]<<" "<<mat2[i+j*8]<<endl;

      double dumd[64];
      for (int i=0; i<8; i++) {
        for (int j=0; j<8; j++) {
          double dum = 0.;
          for (int k=0; k<8; k++) {
            //if (k==i || k==j) continue;
            dum += coul.int3[i+k*8]*coul.int3[k+j*8];
          }
          dumd[i+j*8] = dum;
          cout<<i<<" "<<j<<" coulomb "<<dum<<" "<<coul.int3[i+j*8]<<endl;
        }
      }
 double term1 = 0.;
 temp = 0.;
  for (int i=0; i<8; i++) {
    //if (i==1) continue;
    term1 += coul.evecs3[1+i*8]*coul.evecs3[2+i*8]*coul.evals3[i];
    temp += energies[7]*coul.evecs3[7+i*8]*coul.evecs3[7+i*8];
    }
  cout<<"term 1 "<<term1<<" temp "<<temp<<endl;
}


void pertCalc(Molecule *mol) {
  //define molecular positions
  //calculate couplings
  double e100 = mol[0].excenergy[1];// - mol[0].groundenergy;
  double energies[2][2][2];

  for (int i=0; i<2; i++) 
    for (int j=0; j<2; j++) 
      for (int k=0; k<2; k++) {
        energies[i][j][k] = (mol[0].excenergy[i] )//- mol[0].groundenergy)
                            + (mol[1].excenergy[j]*.99)// - mol[1].groundenergy)
                            + (mol[2].excenergy[k]*1.1);// - mol[2].groundenergy);
      }

  double inv, temp,temp2,r12,r13,r23,sum,pos[3];
  temp2 = 0.;
  temp = 0.;
  sum = 0.;
  
  //I term
  for (int I=0; I<2; I++) {
    inv = 1./(energies[1][0][0] - energies[I][1][0]);

    for (int i=0; i<mol[0].natoms; i++) {
      for (int j=0; j<mol[1].natoms; j++) {
        //get r12
        for (int p=0; p<3; p++) {
          pos[p] = mol[0].atoms[i].pos[p] - mol[1].atoms[j].pos[p];
        }
        r12 = cblas_ddot(3,pos,1,pos,1);
        r12 = sqrt(r12);

//CTC test start
//        temp += 10*mol[0].atoms[i].charges[2] * mol[1].atoms[j].charges[2] / r12;

        for (int k=0; k<mol[0].natoms; k++) {
          for (int l=0; l<mol[2].natoms; l++) {
            //get r13
            for (int p=0; p<3; p++) {
              pos[p] = mol[0].atoms[k].pos[p] - mol[2].atoms[l].pos[p];
            }
            r13 = cblas_ddot(3,pos,1,pos,1);
            r13 = sqrt(r13);

            temp += mol[0].atoms[i].charges[2-I] * mol[1].atoms[j].charges[2]
                    * mol[0].atoms[k].charges[I*2] * mol[2].atoms[l].charges[2]
                    / (r12*r13);
          } //end A2 atoms
        } //end D' atoms
      } //end A1 atoms
    } //end D atoms
//    cout<<temp<<endl;
//    exit(0);
//CTC test end
    temp *= inv;



    inv = 1./(energies[1][0][0] - energies[I][0][1]);
    //cout<<I<<" "<<temp<<" "<<r12*r13<<" "<<inv<<endl;
    for (int i=0; i<mol[0].natoms; i++) {
      for (int j=0; j<mol[2].natoms; j++) {
        //get r13
        for (int p=0; p<3; p++) {
          pos[p] = mol[0].atoms[i].pos[p] - mol[2].atoms[j].pos[p];
        }

        r13 = cblas_ddot(3,pos,1,pos,1);
        r13 = sqrt(r13);
        for (int k=0; k<mol[0].natoms; k++) {
          for (int l=0; l<mol[1].natoms; l++) {
            //get r12
          for (int p=0; p<3; p++) {
            pos[p] = mol[0].atoms[k].pos[p] - mol[1].atoms[l].pos[p];
          }

            r12 = cblas_ddot(3,pos,1,pos,1);
            r12 = sqrt(r12);
 
            temp2 += mol[0].atoms[i].charges[2-I] * mol[2].atoms[j].charges[2]
                    * mol[0].atoms[k].charges[I*2] * mol[1].atoms[l].charges[2]
                    /(r12*r13);

          } //end A2 atoms
        } //end D' atoms

      } //end A1 atoms
    } //end D atoms
    temp2 *= inv;
    sum += temp + temp2;
    cout<<"Sum = "<<sum<<" "<<temp<<" "<<temp2<<endl;
  } //end I

  cout<<"Sum after I = "<<sum<<endl;

  //J term
  temp = 0.; temp2 = 0.;
  for (int J=0; J<2; J++) {
    inv = 1./(energies[1][0][0] - energies[0][J][0]);

    for (int i=0; i<mol[0].natoms; i++) {
      for (int j=0; j<mol[1].natoms; j++) {
        //get r12
        for (int p=0; p<3; p++) {
          pos[p] = mol[0].atoms[i].pos[p] - mol[1].atoms[j].pos[p];
        }

        r12 = cblas_ddot(3,pos,1,pos,1);
        r12 = sqrt(r12);
        for (int k=0; k<mol[1].natoms; k++) {
          for (int l=0; l<mol[2].natoms; l++) {
            //get r23
            for (int p=0; p<3; p++) {
              pos[p] = mol[1].atoms[k].pos[p] - mol[2].atoms[l].pos[p];
            }
            r23 = cblas_ddot(3,pos,1,pos,1);
            r23 = sqrt(r23);
 
            temp += mol[0].atoms[i].charges[2] * mol[1].atoms[j].charges[J*2]
                    * mol[1].atoms[k].charges[2-J] * mol[2].atoms[l].charges[2]
                    / (r12*r23);
          } //end A2 atoms
        } //end A1' atoms
      } //end A1 atoms
    } //end D atoms
    temp *= inv;
  
    inv = 1./(energies[1][0][0] - energies[1][J][1]);
    for (int i=0; i<mol[1].natoms; i++) {
      for (int j=0; j<mol[2].natoms; j++) {
        //get r23
        for (int p=0; p<3; p++) {
          pos[p] = mol[1].atoms[i].pos[p] - mol[2].atoms[j].pos[p];
        }

        r23 = cblas_ddot(3,pos,1,pos,1);
        r23 = sqrt(r23);
        for (int k=0; k<mol[0].natoms; k++) {
          for (int l=0; l<mol[1].natoms; l++) {
            //get r12
            for (int p=0; p<3; p++) {
              pos[p] = mol[0].atoms[k].pos[p] - mol[1].atoms[l].pos[p];
            }

            r12 = cblas_ddot(3,pos,1,pos,1);
            r12 = sqrt(r12);
 
            temp2 += mol[1].atoms[i].charges[2*J] * mol[2].atoms[j].charges[2]
                    * mol[0].atoms[k].charges[2] * mol[1].atoms[l].charges[2-J]
                    / (r12*r23);
          } //end A2 atoms
        } //end D' atoms
      } //end A1 atoms
    } //end D atoms
    temp2 *= inv;
    sum += temp + temp2;
  } //end J

  cout<<"Sum after J = "<<sum<<endl;
  
  //K term
  temp = 0.; temp2 = 0.;
  for (int K=0; K<2; K++) {
    inv = 1./(energies[1][0][0] - energies[0][0][K]);

    for (int i=0; i<mol[0].natoms; i++) {
      for (int j=0; j<mol[2].natoms; j++) {
        //get r13
        for (int p=0; p<3; p++) {
          pos[p] = mol[0].atoms[i].pos[p] - mol[2].atoms[j].pos[p];
        }

        r13 = cblas_ddot(3,pos,1,pos,1);
        r13 = sqrt(r13);
        for (int k=0; k<mol[1].natoms; k++) {
          for (int l=0; l<mol[2].natoms; l++) {
            //get r23
            for (int p=0; p<3; p++) {
              pos[p] = mol[1].atoms[k].pos[p] - mol[2].atoms[l].pos[p];
            }

            r23 = cblas_ddot(3,pos,1,pos,1);
            r23 = sqrt(r23);
 
            temp += mol[0].atoms[i].charges[2] * mol[2].atoms[j].charges[2*K]
                    * mol[1].atoms[k].charges[2] * mol[2].atoms[l].charges[2-K]
                    / (r13*r23);

          } //end A2 atoms
        } //end A1' atoms
      } //end A2 atoms
    } //end D atoms
    temp *= inv;

    inv = 1./(energies[1][0][0] - energies[1][1][K]);
    for (int i=0; i<mol[1].natoms; i++) {
      for (int j=0; j<mol[2].natoms; j++) {
        //get r23
        for (int p=0; p<3; p++) {
          pos[p] = mol[1].atoms[i].pos[p] - mol[2].atoms[j].pos[p];
        }

        r23 = cblas_ddot(3,pos,1,pos,1);
        r23 = sqrt(r23);
        for (int k=0; k<mol[0].natoms; k++) {
          for (int l=0; l<mol[2].natoms; l++) {
            //get r13
            for (int p=0; p<3; p++) {
              pos[p] = mol[0].atoms[k].pos[p] - mol[2].atoms[l].pos[p];
            }

            r13 = cblas_ddot(3,pos,1,pos,1);
            r13 = sqrt(r13);
 
            temp2 += mol[1].atoms[i].charges[2] * mol[2].atoms[j].charges[2*K]
                    * mol[0].atoms[k].charges[2] * mol[2].atoms[l].charges[2-K]
                    / (r13*r23);
          } //end A2' atoms
        } //end D' atoms
      } //end A2 atoms
    } //end A1 atoms
    temp2 *= inv;
    sum += temp + temp2;
  } //end K

  cout<<"Sum after K = "<<sum<<endl;

  //calculate rate constants

  //CTC Temporary stuff
  //<a|b> + <b|a> = - <a|a> - <b|b> 
  
/*  double sum = 0.;
  for (int i=0; i<mol[0].natoms; i++) {
    sum = -1*mol[0].atoms[i].charges[0] 
      +1*mol[0].atoms[i].charges[1];
    cout<<i<<" "<<mol[0].atoms[i].type<<" "<<sum<<endl;
  }
*/


}
#endif //PERT_H

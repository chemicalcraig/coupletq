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
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_blas.h"
using namespace std;

double project(int mol, int state, const double *re, const double *im, const double *evecs) {
  double res = 0;
  double sum = 0.;
  double vecre[8];
  double vecim[8];

  //convert psi to index-basis
  cblas_dgemv(CblasColMajor,CblasNoTrans,8,8,1.,evecs,8,re,1,0.,vecre,1);
  cblas_dgemv(CblasColMajor,CblasNoTrans,8,8,1.,evecs,8,im,1,0.,vecim,1);

  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      int index;
      if (mol == 0)
        index = state+i*2+j*4;
      else if (mol == 1)
        index = i+state*2+j*4;
      else if (mol == 2)
        index  = i+j*2+state*4;
      sum += vecre[index]*vecre[index] + vecim[index]*vecim[index];
    }
  }

  return sum;
}

void propagateTime(Molecule *mol, Coulomb coul, double *energies, double tstart, double tfinish,
                  double dtt, double *intham) {
  double *psire = new double[8];
  double *psiim = new double[8];
  double *dpsire = new double[8];
  double *dpsiim = new double[8];
  double ttime = tstart;

  //construct full hamiltonian
  double *ham, *diagham;
  ham = new double[64];
  diagham = new double[64];
  for (int i=0; i<8; i++) {
    ham[i+i*8] = energies[i]; //coul.int3 already has energies in it
    for (int j=0; j<8; j++) {
      ham[i+j*8] += intham[i+j*8] + coul.int3[i+j*8];
      cout<<i<<" "<<j<<" hams "<<intham[i+j*8]<<" "<<coul.int3[i+j*8]<<" "<<ham[i+j*8]<<endl;
    }
  }

  //diagonalize full hamiltonian
  gsl_matrix_view m = gsl_matrix_view_array(ham,8,8);
  gsl_vector *eval = gsl_vector_alloc(8);
  gsl_matrix *evec = gsl_matrix_alloc(8,8);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(8);
  gsl_eigen_symmv(&m.matrix,eval,evec,w);
  gsl_eigen_symmv_free(w);

  //coul.diagonalize(8,coul.evecs3,coul.evals3,coul.int3);
  double vec[8],vec2[8],evecs[64];
  for (int i=0; i<8; i++) {
    cout<<"Full evals "<<gsl_vector_get(eval,i)<<endl;
  }
  for (int i=0; i<8; i++) {
    for (int j=0; j<8; j++) {
      evecs[i+j*8] = gsl_matrix_get(evec,i,j);
      cout<<"Full evecs "<<i<<" "<<j<<" "<<evecs[i+j*8]<<endl;
    }
  }
  
  //convert initial state to eigenbasis
  for (int i=0; i<8; i++) {
    psire[i] = 0.;
    psiim[i] = 0.;
  }
  psire[1] = 1/sqrt(2);
  psire[0] = 1/sqrt(2);
  //psire[3] = 1;

  double ham2[64],tempm[64],tempm2[64];
  cblas_dgemv(CblasColMajor,CblasTrans,8,8,1.,evecs,8,psire,1,0.,dpsire,1);
  cblas_dgemv(CblasColMajor,CblasTrans,8,8,1.,evecs,8,psiim,1,0.,dpsiim,1);
  
  //project onto selected state
  int state = 0;
  int site = 0;

  double population = project(site,state,dpsire,dpsiim,evecs);
  
  ofstream outfile;
  outfile.open("populations.dat");
  outfile.precision(16);
  
  //propagate wavefunction in time
  double initenergy = 0.;
  double norm = 0.;
  for (int i=0; i<8; i++) {
    initenergy += energies[i]*(psire[i]*psire[i]+psiim[i]*psiim[i]);
 }
  
  int counter = 1;
  //Planck's constant in ev fs
  double planck = 4.135667516;
  
  while (ttime < tfinish) {

    for (int i=0; i<8; i++) {
      double ev = gsl_vector_get(eval,i);
      dpsire[i] = cos(ev*dtt/planck)*dpsire[i] + sin(ev*dtt/planck)*dpsiim[i];
      dpsiim[i] = cos(ev*dtt/planck)*dpsiim[i] - sin(ev*dtt/planck)*dpsire[i];
    }

    if (counter%10000 == 0) {
    population = project(0,1,dpsire,dpsiim,evecs);
    cout<<"population of donor excited state = "<<population<<endl;
    double energy = 0.;
    norm = 0.;
    for (int i=0; i<8; i++) {
      energy += gsl_vector_get(eval,i)*(dpsire[i]*dpsire[i]+dpsiim[i]*dpsiim[i]);
      norm += (dpsire[i]*dpsire[i]+dpsiim[i]*dpsiim[i]);
    }

    outfile<<ttime<<" "<<project(0,0,dpsire,dpsiim,evecs)
      <<" "<<project(1,0,dpsire,dpsiim,evecs)
      <<" "<<project(2,0,dpsire,dpsiim,evecs)
      <<" "<<project(0,1,dpsire,dpsiim,evecs)
      <<" "<<project(1,1,dpsire,dpsiim,evecs)
      <<" "<<project(2,1,dpsire,dpsiim,evecs)
      <<" "<<norm<<endl;
    }
    ttime += dtt;
    counter ++;
  }
}

/** Perturbation Calculation in eigenbasis of the 3-bit 
 * Coulomb operator **/
void pertCalcNonDegen(Molecule *mol, Coulomb coul, double *energies,double *int3,double &dum) {
  cout<<" *** in pertCalc, using constant energy differences *** "<<endl;
  
  /** Scale Coulomb interaction NB: Fix this! **/
  for (int i=0; i<8; i++) { //initial state
    coul.int3[i+i*8] -= energies[i];
    for (int j=0; j<8; j++) { //final state
      if (i!=j)
      coul.int3[i+j*8] /= 1.;
    }
  }
  double res[64];
  
  //Make diagonal hamiltonian with exciton energies on the diagonal
  double sum = 0.;
  for (int i=0; i<8; i++) { //initial state
    for (int j=0; j<8; j++) { //final state
      sum = 0.;
      for (int k=0; k<8; k++) { //intermediate state
        if (i==k) continue;
        if (energies[i] != energies[k])
          sum += coul.int3[i+k*8]*coul.int3[k+j*8]/(energies[i]-energies[k]);
        else //CTC This is an appoximation!
          sum += coul.int3[i+k*8]*coul.int3[k+j*8]/5;//(energies[i]-energies[k]);
          
          //Pring coupling elements
          if (i==3 && j==5) {
            cout<<"<"<<i<<"|V|"<<k<<"><"<<k<<"|V|"<<j<<"> = "
            <<sum<<" "<<coul.int3[i+k*8]<<" "<<coul.int3[k+j*8]<<endl;
          }
        } //end intermediate state
      res[i+j*8] = sum*window(energies[i],energies[j],1,0);
    } //end final state
  } //end initial state
 for (int i=0; i<8; i++)
   for (int j=0; j<8; j++) {
      cout<<i<<" "<<j<<" "<<res[i+j*8]<<"    <"<<i<<"|V|"<<j<<"> = "<<coul.int3[i+j*8]<<endl;
  }
}

void pertCalcDegen(Molecule *mol, Coulomb coul, double *energies,double *int3,double &dum,double *intham) {
  cout<<" *** in pertCalc*** "<<endl;
  
  /** Scale Coulomb interaction NB: Fix this! **/
  for (int i=0; i<8; i++) { //initial state
    //coul.int3[i+i*8] += energies[i];
    for (int j=0; j<8; j++) { //final state
    //if (i!=j)
      //coul.int3[i+j*8] *= 10.;
    }
  }
  double res[64];
  
  //Make diagonal hamiltonian with exciton energies on the diagonal
  double ham[64],tildeint[64];
  for (int i=0; i<64; i++) ham[i] = 0.;
  for (int i=0; i<8; i++) {
    ham[i+i*8] = energies[i];
  }

  //Make C.ham.C^T
  double ham2[64],tempm[64],tempm2[64];
  cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,8,8,8,1.,
              coul.evecs3,8,ham,8,0,tempm,8);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,8,8,8,1.,
              tempm,8,coul.evecs3,8,0,ham2,8);
  //Make C.V.C^T
  cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,8,8,8,1.,
              coul.evecs3,8,coul.int3,8,0,tempm,8);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,8,8,8,1.,
              tempm,8,coul.evecs3,8,0,tildeint,8);

  double sum = 0.;
  double numerator, denominator;
  for (int i=0; i<8; i++) { //initial state
    for (int j=0; j<8; j++) { //final state
      numerator = 0.;
      denominator = 0.;
      for (int k=0; k<8; k++) { //intermediate state
        double sum1 = 0.;
        double sum2 = 0.;
        double sum3 = 0.;
        double sum4 = 0.;
        if (i==k) continue;
        for (int a=0; a<8; a++) {
          //for (int b=0; b<8; b++) {
            //if (a!=b) continue;
            //sum1 += coul.evecs3[i+a*8]*coul.evecs3[b+k*8]*tildeint[a+b*8];
            sum1 += coul.evecs3[i+a*8]*coul.evecs3[k+a*8]*coul.evals3[a];
          //}
        }
        for (int g=0; g<8; g++) {
          //for (int d=0; d<8; d++) {
            //if (g!=d) continue;
            //sum2 += coul.evecs3[k+g*8]*coul.evecs3[d+j*8]*tildeint[g+d*8];     
            sum2 += coul.evecs3[k+g*8]*coul.evecs3[j+g*8]*coul.evals3[g];
          //}
        }

        for (int a=0; a<8; a++) {
          for (int b=0; b<8; b++) {
            sum3 += coul.evecs3[i+a*8]*coul.evecs3[i+b*8]*ham2[a+b*8];
          }
        }
        for (int a=0; a<8; a++) {
          for (int b=0; b<8; b++) {
            sum4 += coul.evecs3[k+a*8]*coul.evecs3[k+b*8]*ham2[a+b*8];
          }
        }
        
        //numerator += coul.int3[i+k*8]*coul.int3[k+j*8]/(coul.evals3[i]-coul.evals3[k]); 
        
        numerator += sum1*sum2;//(coul.evals3[i]-coul.evals3[k]); //(sum3-sum4);
        denominator +=sum3-sum4;// coul.evals3[i] - coul.evals3[k];//sum3-sum4;
      //Pring coupling elements
          if (i==1 && j==6) {
            cout<<"<"<<i<<"|V|"<<k<<"><"<<k<<"|V|"<<j<<"> = "
            <<numerator<<" "<<" "<<denominator<<" "<<tildeint[i+k*8]<<" "<<tildeint[k+j*8]<<endl;
          }
        }
      

      res[i+j*8] = numerator/denominator;
    }
  }
  cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,8,8,8,1.,
              coul.evecs3,8,res,8,0,tempm,8);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,8,8,8,1.,
              tempm,8,coul.evecs3,8,0,tildeint,8);
 for (int i=0; i<8; i++)
   for (int j=0; j<8; j++) {
   cout<<i<<"   "<<j<<"   "<<res[i+j*8]<<"   "
      <<tildeint[i+j*8]*window(energies[i],energies[j],1,0)<<"       <"<<i<<"|V|"<<j<<"> = "<<coul.int3[i+j*8]<<endl;
    
    intham[i+j*8] = res[i+j*8]*window(energies[i],energies[j],0.4,0);
  }
  dum = res[1+6*8];
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

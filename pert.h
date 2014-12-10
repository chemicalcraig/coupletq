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
#include "gsl/gsl_complex.h"
using namespace std;

double project(Molecule *mol,int nmol, int state, const double *re, const double *im, const double *evecs) {
  double res = 0;
  double sum = 0.;
  double vecre[mol[0].nindices];
  double vecim[mol[0].nindices];

  //convert psi to index-basis
  cblas_dgemv(CblasColMajor,CblasNoTrans,mol[0].nindices,
      mol[0].nindices,1.,evecs,mol[0].nindices,re,1,0.,vecre,1);
  cblas_dgemv(CblasColMajor,CblasNoTrans,mol[0].nindices,
      mol[0].nindices,1.,evecs,mol[0].nindices,im,1,0.,vecim,1);

  for (int i=0; i<mol[0].nindices; i++) {
    if (mol[0].indices[nmol+i*mol[0].nmol] == state) {
      sum += vecre[i]*vecre[i] + vecim[i]*vecim[i];
    }
  }


  return sum;
}

void propagateTime(Molecule *mol, Coulomb coul, double *energies, double tstart, double tfinish,
                  double dtt, double *intham, Reader r) {
  double *psire = new double[mol[0].nindices];
  double *psiim = new double[mol[0].nindices];
  double *dpsire = new double[mol[0].nindices];
  double *dpsiim = new double[mol[0].nindices];
  double ttime = tstart;

  //construct zeroth order hamiltonian in site basis
  double *ham;
  ham = new double[mol[0].nindices*mol[0].nindices];
  for (int i=0; i<mol[0].nindices*mol[0].nindices; i++)
    ham[i] = 0.;
  
  for (int i=0; i<mol[0].nindices; i++) {
    ham[i+i*mol[0].nindices] = energies[i];
  }

  //Make C.ham.C^T=ham2
  double ham2[mol[0].nindices*mol[0].nindices];
  double tempm[mol[0].nindices*mol[0].nindices];
  double tildeint[mol[0].nindices*mol[0].nindices];
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,mol[0].nindices,
          mol[0].nindices,mol[0].nindices,1.,coul.evecs3,mol[0].nindices,
          ham,mol[0].nindices,0,tempm,mol[0].nindices);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,mol[0].nindices,
          mol[0].nindices,mol[0].nindices,1.,tempm,mol[0].nindices,
          coul.evecs3,mol[0].nindices,0,ham2,mol[0].nindices);
  //Make C.V.C^T
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,mol[0].nindices,
          mol[0].nindices,mol[0].nindices,1.,coul.evecs3,
          mol[0].nindices,coul.int3,mol[0].nindices,0,tempm,mol[0].nindices);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,mol[0].nindices,
          mol[0].nindices,mol[0].nindices,1.,tempm,mol[0].nindices,
          coul.evecs3,mol[0].nindices,0,tildeint,mol[0].nindices);
  
   for (int i=0; i<mol[0].nindices; i++) {
    for (int j=0; j<mol[0].nindices; j++) {
     // cout<<i<<" "<<j<<" "<<ham2[i+j*mol[0].nindices]<<" "<<tildeint[i+j*mol[0].nindices]<<endl;
    }
  }
  //Zero out ham
  for (int i=0; i<mol[0].nindices*mol[0].nindices; i++)
    ham[i] = 0.;

  //Construct ham in |I> basis
  for (int i=0; i<mol[0].nindices; i++) {
    for (int j=0; j<mol[0].nindices; j++) {
      ham[i+j*mol[0].nindices] += ham2[i+j*mol[0].nindices] + intham[i+j*mol[0].nindices] 
            + tildeint[i+j*mol[0].nindices];
    }
  }
  
  /** Convert ham back into energy basis **/
  cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,mol[0].nindices,
          mol[0].nindices,mol[0].nindices,1.,coul.evecs3,
          mol[0].nindices,ham,mol[0].nindices,0,tempm,mol[0].nindices);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,mol[0].nindices,
          mol[0].nindices,mol[0].nindices,1.,tempm,mol[0].nindices,
          coul.evecs3,mol[0].nindices,0,ham2,mol[0].nindices);
  
  //Filter energies
  for (int i=0; i<mol[0].nindices; i++) {
    for (int j=0; j<mol[0].nindices; j++) {
      ham2[i+j*mol[0].nindices] *= window(energies[i],energies[j],r.calc.ewindow,0);
      cout<<i<<" "<<j<<" hams "<<intham[i+j*mol[0].nindices]<<" "<<tildeint[i+j*mol[0].nindices]<<" "<<ham2[i+j*mol[0].nindices]<<endl;//" "<<ham[i+j*mol[0].nmol]<<endl;
    }
  }
 
  //diagonalize full hamiltonian
  gsl_matrix_view m = gsl_matrix_view_array(ham2,mol[0].nindices,mol[0].nindices); 
  gsl_vector *eval = gsl_vector_alloc(mol[0].nindices);
  gsl_matrix *evec = gsl_matrix_alloc(mol[0].nindices,mol[0].nindices);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(mol[0].nindices);
  gsl_eigen_symmv(&m.matrix,eval,evec,w);
  gsl_eigen_symmv_free(w);

  
  //CTC start complex
/*  gsl_matrix_complex *m1 = gsl_matrix_complex_alloc(mol[0].nindices,mol[0].nindices);
  for (int i=0; i<mol[0].nindices; i++) {
      for (int j=0; j<mol[0].nindices; j++) {
        //m1->data[i].dat[0] = m.matrix.data[i];
        gsl_complex comp;
        comp.dat[0] = gsl_matrix_get(&m.matrix,i,j);
        gsl_matrix_complex_set(m1,i,j,comp);
      }
  }
  gsl_vector_complex *evalc = gsl_vector_complex_alloc(mol[0].nindices);
  gsl_matrix_complex *evec1 = gsl_matrix_complex_alloc(mol[0].nindices,mol[0].nindices);
  gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(mol[0].nindices);
  gsl_eigen_nonsymmv_params(0,w);
  gsl_eigen_nonsymmv(&m.matrix,evalc,evec1,w);
  gsl_eigen_nonsymmv_free(w);
  
  for (int i=0; i<mol[0].nindices; i++) {
    gsl_complex comp;
    for (int j=0; j<mol[0].nindices; j++) {
        //m1->data[i].dat[0] = m.matrix.data[i];
        comp = gsl_matrix_complex_get(evec1,i,j);
        cout<<"re,im part of ev_{"<<i<<","<<j<<"} = ("<<comp.dat[0]<<", "<<comp.dat[1]<<")"<<endl;
        gsl_matrix_complex_set(m1,i,j,comp);
        gsl_matrix_set(evec,i,j,comp.dat[0]);
      }
      comp = gsl_vector_complex_get(evalc,i);
      gsl_vector_set(eval,i,comp.dat[0]);
  }
*/  //CTC end complex
  
  double vec[mol[0].nindices],vec2[mol[0].nindices],evecs[mol[0].nindices*mol[0].nindices];
  for (int i=0; i<mol[0].nindices; i++) {
    cout<<"Full evals "<<gsl_vector_get(eval,i)<<endl;
  }
  for (int i=0; i<mol[0].nindices; i++) {
    for (int j=0; j<mol[0].nindices; j++) {
      evecs[i+j*mol[0].nindices] = gsl_matrix_get(evec,i,j);
      cout<<"Full evecs "<<i<<" "<<j<<" "<<evecs[i+j*mol[0].nindices]<<endl;
    }
  }
  
  //convert initial state to eigenbasis
  for (int i=0; i<mol[0].nindices; i++) {
    psire[i] = 0.;
    psiim[i] = 0.;
  }
  psire[1] = 1.;
  //psire[0] = 1/sqrt(2);
  //psire[3] = 1;

   cblas_dgemv(CblasColMajor,CblasTrans,mol[0].nindices,
       mol[0].nindices,1.,evecs,mol[0].nindices,psire,1,0.,dpsire,1);
  cblas_dgemv(CblasColMajor,CblasTrans,mol[0].nindices,
      mol[0].nindices,1.,evecs,mol[0].nindices,psiim,1,0.,dpsiim,1);

  for (int i=0; i<mol[0].nindices; i++) {
    cout<<"initial vector "<<i<<" "<<dpsire[i]<<" "<<dpsiim[i]<<endl;
  }
  
  //project onto selected state
  int state = 0;
  int site = 0;
  double population;
  
  ofstream outfile;
  outfile.open("populations.dat");
  outfile.precision(16);
  
  //propagate wavefunction in time
  double initenergy = 0.;
  double norm = 0.;
  for (int i=0; i<mol[0].nindices; i++) {
    initenergy += energies[i]*(psire[i]*psire[i]+psiim[i]*psiim[i]);
 }
  
  int counter = 1;
  //Planck's constant in ev fs
  double planck = 1;//4.135667516;
  
  //convert dt from au to fs
  //dtt *= .02418884326505;

  while (ttime < tfinish) {

    for (int i=0; i<mol[0].nindices; i++) {
      double ev = gsl_vector_get(eval,i);
      dpsire[i] = cos(ev*dtt/planck)*dpsire[i] + sin(ev*dtt/planck)*dpsiim[i];
      dpsiim[i] = cos(ev*dtt/planck)*dpsiim[i] - sin(ev*dtt/planck)*dpsire[i];
    }

    if (counter%10000 == 0) {
      population = project(mol,0,1,dpsire,dpsiim,evecs);
      cout<<"population of donor excited state = "<<population<<endl;
      double energy = 0.;
      norm = 0.;
      for (int i=0; i<mol[0].nindices; i++) {
        energy += gsl_vector_get(eval,i)*(dpsire[i]*dpsire[i]+dpsiim[i]*dpsiim[i]);
        norm += (dpsire[i]*dpsire[i]+dpsiim[i]*dpsiim[i]);
      }

      outfile<<ttime*.02418884326505<<" ";
      for (int m=0; m<mol[0].nmol; m++) {
        for (int st=0; st<mol[m].nstates; st++) {
          outfile<<project(mol,m,st,dpsire,dpsiim,evecs)<<" ";
        }
      }
      outfile<<" "<<norm<<" "<<energy<<endl;
    }
    ttime += dtt;
    counter ++;
  }
}

/** Perturbation Calculation in eigenbasis of the 3-bit 
 * Coulomb operator **/
void pertCalcEigen(Molecule *mol, Coulomb coul, double *energies,double *int3,double *intham) {
  
  cout<<" *** in pertCalc, using eigenbasis of V *** "<<endl;
  double sum = 0.;
  double sum1=0.;
  
  cout<<"/** Second order corrections in eigenbasis of V **/"<<endl;
  cout<<244<<endl;
  cout<<mol[0].nindices<<endl;
  for (int i=0; i<mol[0].nindices; i++) { //initial state
    for (int j=0; j<mol[0].nindices; j++) { //final state
      sum = 0.;
      for (int k=0; k<mol[0].nindices; k++) { //intermediate state
        if (i==k)  continue;
        cout<<i<<" "<<j<<" "<<mol[0].nindices<<endl;
        sum1=0.;
        for (int a=0; a<mol[0].nindices; a++) {
          for (int b=0; b<mol[0].nindices; b++) {
            for (int c=0; c<mol[0].nindices; c++) {
              for (int d=0; d<mol[0].nindices; d++) {
                double dum = coul.evecs3[a+i*mol[0].nindices]
                    *coul.int3[a+b*mol[0].nindices]
                    *coul.evecs3[b+k*mol[0].nindices]
                    *coul.evecs3[c+k*mol[0].nindices]
                    *coul.int3[c+d*mol[0].nindices]
                    *coul.evecs3[d+j*mol[0].nindices];

                sum1 += dum;
              }
            }
          }
        }
        double en = coul.evals3[i] - coul.evals3[k];
        sum += sum1/en;
      }
      intham[i+j*mol[0].nindices] = sum;
    }//end final state
  }//end initial state
}

void pertCalcDegen(Molecule *mol, Coulomb coul, double *energies,double *int3,double &dum,double *intham, Reader r) {
  cout<<" *** in pertCalc*** "<<endl;
  
  /** Scale Coulomb interaction NB: Fix this! **/
  double res[mol[0].nindices*mol[0].nindices];
  
  //Make diagonal hamiltonian with exciton energies on the diagonal
  double ham[mol[0].nindices*mol[0].nindices],tildeint[mol[0].nindices*mol[0].nindices];
  
  for (int i=0; i<mol[0].nindices*mol[0].nindices; i++) ham[i] = 0.;
  
  for (int i=0; i<mol[0].nindices; i++) {
    ham[i+i*mol[0].nindices] = energies[i];
  }

  //Make C.ham.C^T
  double ham2[mol[0].nindices*mol[0].nindices],tempm[mol[0].nindices*mol[0].nindices],tempm2[mol[0].nindices*mol[0].nindices];
  cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,mol[0].nindices,mol[0].nindices,mol[0].nindices,1.,
              coul.evecs3,mol[0].nindices,ham,mol[0].nindices,0,tempm,mol[0].nindices);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,mol[0].nindices,mol[0].nindices,mol[0].nindices,1.,
              tempm,mol[0].nindices,coul.evecs3,mol[0].nindices,0,ham2,mol[0].nindices);
  //Make C.V.C^T
  cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,mol[0].nindices,mol[0].nindices,mol[0].nindices,1.,
              coul.evecs3,mol[0].nindices,coul.int3,mol[0].nindices,0,tempm,mol[0].nindices);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,mol[0].nindices,mol[0].nindices,mol[0].nindices,1.,
              tempm,mol[0].nindices,coul.evecs3,mol[0].nindices,0,tildeint,mol[0].nindices);

  double sum = 0.;
  double numerator, denominator;
  for (int i=0; i<mol[0].nindices; i++) { //initial state
    for (int j=0; j<mol[0].nindices; j++) { //final state
      numerator = 0.;
      denominator = 0.;
      for (int k=0; k<mol[0].nindices; k++) { //intermediate state
        double sum1 = 0.;
        double sum2 = 0.;
        double sum3 = 0.;
        double sum4 = 0.;
        if (i==k) continue;
        for (int a=0; a<mol[0].nindices; a++) {
          //for (int b=0; b<mol[0].nindices; b++) {
            //if (a!=b) continue;
            //sum1 += coul.evecs3[i+a*mol[0].nindices]*coul.evecs3[b+k*mol[0].nindices]*tildeint[a+b*mol[0].nindices];
            sum1 += coul.evecs3[i+a*mol[0].nindices]*coul.evecs3[k+a*mol[0].nindices]*coul.evals3[a];
          //}
        }
        for (int g=0; g<mol[0].nindices; g++) {
          //for (int d=0; d<mol[0].nindices; d++) {
            //if (g!=d) continue;
            //sum2 += coul.evecs3[k+g*mol[0].nindices]*coul.evecs3[d+j*mol[0].nindices]*tildeint[g+d*mol[0].nindices];     
            sum2 += coul.evecs3[k+g*mol[0].nindices]*coul.evecs3[j+g*mol[0].nindices]*coul.evals3[g];
          //}
        }

        for (int a=0; a<mol[0].nindices; a++) {
          for (int b=0; b<mol[0].nindices; b++) {
            sum3 += coul.evecs3[i+a*mol[0].nindices]*coul.evecs3[i+b*mol[0].nindices]*ham2[a+b*mol[0].nindices];
          }
        }
        for (int a=0; a<mol[0].nindices; a++) {
          for (int b=0; b<mol[0].nindices; b++) {
            sum4 += coul.evecs3[k+a*mol[0].nindices]*coul.evecs3[k+b*mol[0].nindices]*ham2[a+b*mol[0].nindices];
          }
        }
        
        //numerator += coul.int3[i+k*mol[0].nindices]*coul.int3[k+j*mol[0].nindices]/(coul.evals3[i]-coul.evals3[k]); 
        
        numerator += (sum1*sum2)/(sum3-sum4);//(coul.evals3[i]-coul.evals3[k]); //(sum3-sum4);
        denominator +=sum3-sum4;// coul.evals3[i] - coul.evals3[k];//sum3-sum4;
      //Pring coupling elements
          if (i==2 && j==0) {
            cout<<"<"<<i<<"|V|"<<k<<"><"<<k<<"|V|"<<j<<"> = "
            <<numerator<<" "<<" "<<denominator<<" "<<numerator/denominator<<" "<<sum1<<" "<<sum2<<" "<<sum3<<" "<<sum4<<endl;
          }
        }
      

      res[i+j*mol[0].nindices] = numerator;///denominator;
    }
  }
 
  for (int i=0; i<mol[0].nindices; i++)
   for (int j=0; j<mol[0].nindices; j++) {
   cout<<i<<"   "<<j<<"   "<<res[i+j*mol[0].nindices]<<"   <"
      <<i<<"|V|"<<j<<"> = "<<coul.int3[i+j*mol[0].nindices]<<endl;
    
    intham[i+j*mol[0].nindices] = res[i+j*mol[0].nindices];
                     // *window(energies[i],energies[j],r.calc.ewindow,0);
  }
  dum = res[1+6*mol[0].nindices];
}

void pertCalc(Molecule *mol, Coulomb coul,double *intham,double *energies) {

  double inv, temp,temp2,temp3,r12,r13,r23,sum,sum2,pos[3];
  temp2 = 0.;
  temp = 0.;
  sum = 0.;
  sum2 = 0.;
  temp3 = 0.;

  for (int i=0; i<mol[0].nindices; i++) { //initial
    for (int j=0; j<mol[0].nindices; j++) {//final
      temp2 = 0.;
      temp3 = 0.;
      for (int k=0; k<mol[0].nindices; k++) {
        temp2 += energies[i]*coul.evecs3[i+k*mol[0].nindices]*coul.evecs3[i+k*mol[0].nindices];
        temp3 += energies[j]*coul.evecs3[j+k*mol[0].nindices]*coul.evecs3[j+k*mol[0].nindices];
      }
      sum = 0.;
      sum2 = 0.;
      double energy;
      //if (i != j) {
        for (int k=0; k<mol[0].nindices; k++) {//intermediate zero order basis
          if (i==k ) continue;
            for (int l=0; l<mol[0].nindices; l++) {
              //i term
              energy = (coul.evals3[i]-coul.evals3[k]);//*(coul.evals3[j]-coul.evals3[k]);
              temp = coul.evecs3[i+l*mol[0].nindices]*coul.evecs3[k+l*mol[0].nindices];
              sum += temp*coul.evals3[l];//energy;
              
              //j term
                     //if ((i == 7 && j==7)) cout<<i<<" "<<j<<" "<<k<<" "<<l<<" sum adding "<<sum<<" "<<sum2<<" "<<temp<<endl;
            }
            for (int l=0; l<mol[0].nindices; l++) {
              temp = coul.evecs3[j+l*mol[0].nindices]*coul.evecs3[k+l*mol[0].nindices];
              sum2 += temp*coul.evals3[l];
 
            }
          }
     // } else {
      /*  for (int k=0; k<mol[0].nindices; k++) {
          if (i==k) continue;
          energy = (coul.evals3[k]-coul.evals3[i])*(coul.evals3[k]-coul.evals3[i]);
          sum += coul.evecs3[i+k*mol[0].nindices]*coul.evecs3[i+k*mol[0].nindices]
            *coul.evecs3[i+k*mol[0].nindices]*coul.evecs3[i+k*mol[0].nindices]
            *coul.evals3[k]*coul.evals3[k]/energy;
        }
        sum *= -0.5;*/
      //}
      
      intham[i+j*mol[0].nindices] = 2*sum*sum2*window(energies[i],energies[j],1,0);
          }
  }
  double mat[mol[0].nindices*mol[0].nindices],mat2[mol[0].nindices*mol[0].nindices];
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,mol[0].nindices,mol[0].nindices,mol[0].nindices,1.,
              coul.evecs3,mol[0].nindices,intham,mol[0].nindices,0,mat,mol[0].nindices);
      cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,mol[0].nindices,mol[0].nindices,mol[0].nindices,1.,
              mat,mol[0].nindices,coul.evecs3,mol[0].nindices,0,mat2,mol[0].nindices);
      for (int i=0; i<mol[0].nindices; i++) 
        for (int j=0; j<mol[0].nindices; j++)
        cout<<i<<" "<<j<<" "<<intham[i+j*mol[0].nindices]<<" "<<mat2[i+j*mol[0].nindices]<<endl;

      double dumd[mol[0].nindices*mol[0].nindices];
      for (int i=0; i<mol[0].nindices; i++) {
        for (int j=0; j<mol[0].nindices; j++) {
          double dum = 0.;
          for (int k=0; k<mol[0].nindices; k++) {
            //if (k==i || k==j) continue;
            dum += coul.int3[i+k*mol[0].nindices]*coul.int3[k+j*mol[0].nindices];
          }
          dumd[i+j*mol[0].nindices] = dum;
          cout<<i<<" "<<j<<" coulomb "<<dum<<" "<<coul.int3[i+j*mol[0].nindices]<<endl;
        }
      }
 double term1 = 0.;
 temp = 0.;
  for (int i=0; i<mol[0].nindices; i++) {
    //if (i==1) continue;
    term1 += coul.evecs3[1+i*mol[0].nindices]*coul.evecs3[2+i*mol[0].nindices]*coul.evals3[i];
    temp += energies[7]*coul.evecs3[7+i*mol[0].nindices]*coul.evecs3[7+i*mol[0].nindices];
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

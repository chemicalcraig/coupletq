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
  /** Convert time parameters from fs to au **/
  tstart /= .02418884326505;
  tfinish /= .02418884326505;
  dtt /= .02418884326505;
  
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
    for (int j=0; j<mol[0].nindices; j++) {
      coul.int3[i+j*mol[0].nindices] *= window(energies[i],energies[j],r.calc.ewindow,0);
    }
  }

  double ham2[mol[0].nindices*mol[0].nindices];
  double tempm[mol[0].nindices*mol[0].nindices];
  double tildeint[mol[0].nindices*mol[0].nindices];
  
  /** convert W into site basis **/
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,mol[0].nindices,
          mol[0].nindices,mol[0].nindices,1.,coul.evecs3,
          mol[0].nindices,intham,mol[0].nindices,0,tempm,mol[0].nindices);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,mol[0].nindices,
          mol[0].nindices,mol[0].nindices,1.,tempm,mol[0].nindices,
          coul.evecs3,mol[0].nindices,0,tildeint,mol[0].nindices);
   
   for (int i=0; i<mol[0].nindices; i++) {
    for (int j=0; j<mol[0].nindices; j++) {
      tildeint[i+j*mol[0].nindices] *= window(energies[i],energies[j],r.calc.ewindow,0);
    }
  }
  
  //Zero out ham
  for (int i=0; i<mol[0].nindices*mol[0].nindices; i++)
    ham[i] = 0.;
    
  for (int i=0; i<mol[0].nindices; i++) 
      ham[i+i*mol[0].nindices] = energies[i];

  //Construct ham in |I> basis
  for (int i=0; i<mol[0].nindices; i++) {
    for (int j=0; j<mol[0].nindices; j++) {
      ham[i+j*mol[0].nindices] += coul.int3[i+j*mol[0].nindices] + tildeint[i+j*mol[0].nindices];
    }
  }
  
  //Filter energies
  for (int i=0; i<mol[0].nindices; i++) {
    for (int j=0; j<mol[0].nindices; j++) {
      cout<<i<<" "<<j<<" hams "<<tildeint[i+j*mol[0].nindices]
            <<" "<<ham[i+j*mol[0].nindices]<<endl;
    }
  }
 
  //diagonalize full hamiltonian
  gsl_matrix_view m = gsl_matrix_view_array(ham,mol[0].nindices,mol[0].nindices); 
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
  
  ofstream outfile,outfile2;
  outfile.open(r.dyn.out[0].file.c_str());
  outfile2.open("evpops.dat");
  outfile.precision(16);
  outfile2.precision(16);
  
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

    /** Write stuff **/
    if (counter%r.dyn.wstep == 0) {
      population = project(mol,0,1,dpsire,dpsiim,evecs);
      double energy = 0.;
      norm = 0.;
      for (int i=0; i<mol[0].nindices; i++) {
        energy += gsl_vector_get(eval,i)*(dpsire[i]*dpsire[i]+dpsiim[i]*dpsiim[i]);
        norm += (dpsire[i]*dpsire[i]+dpsiim[i]*dpsiim[i]);
      }

      outfile<<ttime*.02418884326505<<" ";
      outfile2<<ttime*.02418884326505<<" ";

      if (r.dyn.printAll) {
        for (int m=0; m<mol[0].nmol; m++) {
          for (int st=0; st<mol[m].nstates; st++) {
            outfile<<project(mol,m,st,dpsire,dpsiim,evecs)<<" ";
          }
        }
        outfile<<" "<<norm<<" "<<energy<<endl;
      } else {
        for (int m=0; m<r.dyn.noutput; m++) {
          outfile<<project(mol,r.dyn.out[m].mol,r.dyn.out[m].state,
                          dpsire,dpsiim,evecs)<<" ";
        }
        outfile<<" "<<norm<<" "<<energy<<endl;
      }
      
      /** Write eigen state populations **/
    /*  for (int m=0; m<mol[0].nindices; m++) {
        double dum = (dpsire[m]*dpsire[m]+dpsiim[m]*dpsiim[m]);
        outfile2<<dum<<" ";
      }
      outfile2<<endl;
    */
    }
    ttime += dtt;
    counter ++;
  }
  outfile.close();
  outfile2.close();
}

/** Get individual elements of the perturbation operator **/
double pertCalcElements(Molecule *mol, Coulomb coul, double *int3, double *energies) {
  cout<<" *** Elements of W2_SSSF *** "<<endl;
  double sum=0;
  double dum;
  /** Sum over intermediate states **/
  for (int i=0; i<mol[0].nindices; i++) { //qb basis
    if (i==1) continue;
    if (mol[0].spinAllowed) {
      if (i==2 || i==3 || i==4 || i==5 )
        continue;
    }
      dum = coul.int3[i+1*mol[0].nindices]
          *coul.int3[6+i*mol[0].nindices];
      sum += dum/(energies[1]-energies[i]);
      cout<<"<1|V|"<<i<<"><"<<i<<"|V|6>/dE = "
        <<27211*dum/(energies[1]-energies[i])
        <<", <1|V"<<i<<"><"<<i<<"|V|6> = "<<27211*dum<<" "
        <<coul.int3[i+1*mol[0].nindices]*27211<<" "
        <<coul.int3[6+i*mol[0].nindices]*27211<<endl;
    }
  cout<<sum*27211<<endl;
  return sum;
}

/** Perturbation Calculation in eigenbasis of the 3-bit 
 * Coulomb operator **/
void pertCalcEigen(Molecule *mol, Coulomb coul, double *energies,double *int3,double *intham) {
  
  cout<<" *** in pertCalcEigen, using eigenbasis of V *** "<<endl;
  double sum = 0.;
  double sum1=0.;
  double temp[mol[0].nindices*mol[0].nindices];

  cout<<"/** Second order corrections in eigenbasis of V **/"<<endl;
  for (int i=0; i<mol[0].nindices; i++) { //initial state
    for (int j=0; j<mol[0].nindices; j++) { //final state
      sum = 0.;
      for (int k=0; k<mol[0].nindices; k++) { //intermediate state
        if (i==k)  continue;
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
        if (i==1 && j==6)
        cout<<"<1|V|"<<k<<"> = "<<sum<<" , <"<<k<<"|V|6> = "<<27211*sum1/en<<endl;
        sum += sum1/en;
     }

      temp[i+j*mol[0].nindices] = sum;
    }//end final state
  }//end initial state

  /** Build W2 in qubit basis **/
  double dumb;
  for (int i=0; i<mol[0].nindices; i++) { //initial state - qb
    for (int j=0; j<mol[0].nindices; j++) { //final state - qb
      sum = 0.;
      for (int k=0; k<mol[0].nindices; k++) { //initial state - ev
        for (int l=0; l<mol[0].nindices; l++) { //final state - ev
          dumb = temp[k+l*mol[0].nindices]
                  *coul.evecs3[i+k*mol[0].nindices]
                  *coul.evecs3[j+l*mol[0].nindices];
          if (i==1 && j==6)
            cout<<"making W2-SSSF "<<k<<" "<<l<<" "<<dumb<<endl;
          sum += dumb;
        }
      }
      intham[i+j*mol[0].nindices] = sum;
      if (i==1 && j==6) cout<<"W2-SSSF "<<sum*27211<<" meV"<<endl;
    }
  }  
  //exit(0);
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


#endif //PERT_H

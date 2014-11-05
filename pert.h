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
#include <iomanip>
using namespace std;

void pertCalc(Molecule *mol) {
  //define molecular positions
  //calculate couplings
  double e100 = mol[0].excenergy[1];// - mol[0].groundenergy;
  double energies[2][2][2];

  for (int i=0; i<2; i++) 
    for (int j=0; j<2; j++) 
      for (int k=0; k<2; k++) {
        energies[i][j][k] = (mol[0].excenergy[i] )//- mol[0].groundenergy)
                            + (mol[1].excenergy[j]*1.1)// - mol[1].groundenergy)
                            + (mol[2].excenergy[k]*.9);// - mol[2].groundenergy);
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
    temp *= inv;

//CTC test end

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

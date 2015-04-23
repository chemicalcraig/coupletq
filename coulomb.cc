#include "coulomb.h"

Coulomb::Coulomb(const int n) {
  this->n2d = 4;
  n3d = n;

  //2-bit stuff
  this->int2 = new double[this->n2d*this->n2d];
  this->evals2 = new double[this->n2d*this->n2d];
  this->evecs2 = new double[this->n2d*this->n2d];
  this->dint2 = new double[this->n2d*this->n2d];
  
  //3-bit stuff
  int3 = new double[n*n];
  this->evals3 = new double[n*n];
  evecs3 = new double[n*n];
  dint3 = new double[n*n];
}

void Coulomb::reinitialize(const int n) {
  this->n2d = 4;
  n3d = n;

  //2-bit stuff
  this->int2 = new double[this->n2d*this->n2d];
  this->evals2 = new double[this->n2d*this->n2d];
  this->evecs2 = new double[this->n2d*this->n2d];
  this->dint2 = new double[this->n2d*this->n2d];
  
  //3-bit stuff
  int3 = new double[n*n];
  this->evals3 = new double[n*n];
  evecs3 = new double[n*n];
  dint3 = new double[n*n];

}

/** Diagonalize the coulomb matrix **/
void Coulomb::diagonalizenonsymm(int nd, double *evc, double *evl, double *mat) {
cout.precision(10);
  /** Copy 3-bit interaction into evecs **/
  for (int i=0; i<nd; i++) {
    for (int j=0; j<nd; j++) {
      //if (i==j)
        evc[i+j*nd] = mat[i+j*nd];
      //else
      //  this->evecs3[i+j*this->n3d] = 0;
    }
  }
  
  //Stuff for lapack diagonalization
  char balanc = 'B';
  char jobvl = 'N';
  char jobvr = 'V';
  char sense = 'V';
  int lwork = 100000000;
  double *evalr = new double [nd];
  double *evali = new double [nd];
  double *work = new double[lwork];
  int liwork = 10000000;
  int *iwork = new int [liwork];
  int info; 
  double leftevec[nd*nd];
  double rightevec[nd*nd];
  int ilo,ihi;
  double scale[nd];
  double abnrm;
  double rcond[nd];
  double rcondv[nd];


  /** do the diagonalizing, with evals and evecs **/
//  dgeevx_(&balanc,&jobvl,&jobvr,&sense,&nd,evc,&nd,evalr,evali,leftevec,&nd,rightevec,
//          &nd,&ilo,&ihi,scale,&abnrm,rcond,rcondv,work,&lwork,iwork,&info);
  //dsyevd_(&jobz, &uplo, &nd, evc, &nd,evl, work, &lwork,
	//          iwork, &liwork, &info);
    for (int i=0; i<nd; i++) {
    for (int j=0; j<nd; j++) {
        evc[i+j*nd] = rightevec[i+j*nd];
    }
  }
  

  for (int i=0; i<nd; i++) {
  cout<<"eval "<<i<<" = "<<evalr[i]<<" "<<evali[i]<<endl;
//    for (int j=0; j<nd; j++) {
//        cout<<i<<" "<<j<<" "<<evc[i+j*nd]<<endl;
//    }
  }

  delete[] work;
  delete[] iwork;
  

}


/** Diagonalize the coulomb matrix **/
void Coulomb::diagonalize(int nd, double *evc, double *evl, double *mat) {
cout.precision(10);
  /** Copy 3-bit interaction into evecs **/
  for (int i=0; i<nd; i++) {
    for (int j=0; j<nd; j++) {
      //if (i==j)
        evc[i+j*nd] = mat[i+j*nd];
      //else
      //  this->evecs3[i+j*this->n3d] = 0;
    }
  }
  
  //Stuff for lapack diagonalization
  char jobz = 'V';
  char uplo = 'U';
  int lwork = 100000000;
  double *work = new double [lwork];
  int liwork = 10000000;
  int *iwork = new int [liwork];
  int info; 

  /** do the diagonalizing, with evals and evecs **/
//  dsyevd_(&jobz, &uplo, &nd, evc, &nd,evl, work, &lwork,
//	          iwork, &liwork, &info);
  delete[] work;
  delete[] iwork;
  for (int i=0; i<nd; i++) {
    double sum = 0.;
    for (int j=0; j<8; j++) {
      sum += mat[i+j*8]*evc[0+j*8];
    }
    cout<<"eval "<<i<<" = "<<evl[i]<<", info = "<<info<<" "<<sum<<endl;
//    for (int j=0; j<nd; j++) {
//        cout<<i<<" "<<j<<" "<<evc[i+j*nd]<<endl;
//    }
  }
}

/** Calculate Coulomb coupling between chromophores I and J
 * undergoing transisions between i and j, and k and l, respectively **/
double Coulomb::getCoulomb(Molecule *mol, int I, int J, int i, int j, int k, int l) {
  double res, sum, temp, pos[3];
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
    }
  }

  return temp;
}

/** Create Coulomb Matrix for two chormophores**/
void Coulomb::createCoulomb2(Molecule *mol) {
  double sum;

  for (int i=0; i<mol[0].nstates; i++) {
    for (int j=0; j<mol[1].nstates; j++) {
      for (int k=0; k<mol[0].nstates; k++) {
        for (int l=0; l<mol[1].nstates; l++) {
          int index = i + j*2 + k*4 + l*8;
          if ((i+j) != 1 || (k+l) != 1)
            this->int2[index] = 0;
            else
            this->int2[index] = getCoulomb(mol,0,1,i,k,j,l);  
          this->int2[15] = getCoulomb(mol,0,1,1,1,1,1);
          this->int2[0] = getCoulomb(mol,0,1,0,0,0,0);
          cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<index<<" "<<int2[index]<<endl;
        }
      }
    }
  }
}

/** Create Coulomb Matrix for three chormophores**/
void Coulomb::createCoulomb3(Molecule *mol) {
  double sum;

  for (int i=0; i<mol[0].nstates; i++) {
    for (int j=0; j<mol[1].nstates; j++) {
      for (int k=0; k<mol[2].nstates; k++) {
        for (int l=0; l<mol[0].nstates; l++) {
          for (int m=0; m<mol[1].nstates; m++) {
            for (int n=0; n<mol[2].nstates; n++) {
              int index = i + j*2 + k*4 + l*8 + m*16 + n * 32;
                int3[index] = getCoulomb(mol,0,1,i,l,j,m) * Kronecker(k,n);
                int3[index] += getCoulomb(mol,0,2,i,l,k,n) * Kronecker(j,m);
                int3[index] += getCoulomb(mol,1,2,j,m,k,n) * Kronecker(i,l);
              cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<n<<" "<<index<<" index"<<endl;
            }
          }
        }
      }
    }
  }
}


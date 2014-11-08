#ifndef COULOMB_H
#define COULOMB_H

#include "molecule.h"

/** Kronecker Delta **/
#define Kronecker(a,b) ((a == b ? 1 : 0))
extern "C"  int dsyevd_(char *jobz, char *uplo, int *n, double *
    a, int *lda, double *w, double *work, int *lwork,
     int *iwork, int *liwork, int *info);

using namespace std;

class Coulomb {
  public:
    int n2d,n3d;
    double *int2, *int3, *dint3, *dint2; //Coulomb coupling matrices, and diagonal
    double *evals3, *evecs3; //evecs and vals for 3-bit system
    double *evals2, *evecs2; //evecs and vals for 2-bit system

    /** Functions **/
    double getCoulomb(Molecule *mol, int I, int J, int i, int j, int k, int l);
    void createCoulomb3(Molecule *mol);
    void createCoulomb2(Molecule *mol);
    void diagonalize(int nd, double *evc, double *evl, double *mat);
    
    Coulomb();
    ~Coulomb() {};

};

#endif

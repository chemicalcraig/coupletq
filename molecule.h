#ifndef MOLECULE_H
#define MOLECULE_H

#include "atom.h"
using namespace std;

class Molecule {
  public:
  int nstates;
  int nao,nbasis;
  int nmo;
  int nroots;
  int nocc,nuocc;
  double dipx,dipy,dipz;
  int nlindep;
  string excMethod;
  double *activeCharges;

  double *ci; //ci coefficients
  int* nbasisatom; //number of basis functions on an atom
  string *nbasisatomorbitals;
  string *nbasisatomelements; //which element basis function belongs to
  double *overlapm;
  double *excenergy;
  double *transmoment;
  double *oscstrength;
  int * occupation;
  double *moeigenv;
  double *mos;

  double *posx,*posy,*posz;
  Atom *atoms;
  int natoms;
  int nx,ny,nz;
  double dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3;
  double ox,oy,oz;
  double ***dens, ***densrad;

  void setnatoms(int n) {natoms = n;};
  void setxyz(int x, int y, int z) {nx=x;ny=y;nz=z;};

  Molecule();
  ~Molecule() {};

  /* Memory allocation routines */
  void allocateMem(const int nb); 
  void allocateMemLindep(const int nb, const int nlindep);
  void allocateMemAtoms(const int na);
  void allocateMemTddft();
  void setnroots(const int n) {this->nroots=n;};
};

#endif // MOLECULE_H

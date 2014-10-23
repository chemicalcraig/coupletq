#ifndef ATOM_H
#define ATOM_H
#include <cstring>
#include <string>
using namespace std;

class Atom {
  public:
  int num,nao,nocc,nuocc,nbasis;
  int atomicnum;
  int *basisfuncs;
  string type;
  double x,y,z;
  double r;
  double tq,tqindo;
  double charge;
  double *charges;
  double mass;
  double pos[3];
  double ipos[3];

  void allocateCharges(const int n) {this->charges = new double[n];};
  
  Atom();
  Atom(int nbf);

};


#endif // ATOM_H

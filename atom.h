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

  Atom();
  Atom(int nbf);
};


#endif // ATOM_H

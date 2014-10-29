#include "atom.h"

Atom::Atom()
{
}

Atom::Atom(int nbf) {
  nbasis = nbf;
  basisfuncs = new int[nbf];
}

/** Set atomic mass **/
void Atom::setMass() {
  for (int i=0; i<119; i++) {
    if (this->type == str_Atom[i]) {
      this->mass = atomicMasses[i];
      continue;
    }
  }
}

void Atom::setMass(string m) {
  this->mass = atof(m.c_str());

}

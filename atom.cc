#include "atom.h"

Atom::Atom()
{
}

Atom::Atom(int nbf) {
  nbasis = nbf;
  basisfuncs = new int[nbf];
}

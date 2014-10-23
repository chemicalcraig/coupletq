#ifndef FRET_H
#define FRET_H

/**************************************
 * This file contains routines
 * to perform a FRET calculation.
 * ***********************************/

#include <cstring>
#include <fstream>
#include <iostream>
#include "atom.h"
#include "molecule.h"
#include <iomanip>
#include "util.h"
using namespace std;

void fretCalc(Molecule *mol, double &couple) {
  //define molecular positions either according to 
  //preferred dipole orientation or cartesian molecular axes

  //calculate couplings
  couple = computeCoupling(mol[0], mol[1], 1,1);

  //calculate rate constants
}
#endif //FRET_H

#ifndef PRINT_H
#define PRINT_H

/**************************************
 * This file contains routines
 * to print data f
 * ***********************************/

#include <cstring>
#include <fstream>
#include <iostream>
#include "atom.h"
#include "molecule.h"
#include <iomanip>
using namespace std;

class Print {
  public:
    Molecule *mol; 


  Print(Molecule *mol);
  ~Print() {}
};

#endif

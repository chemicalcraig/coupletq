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
  void positions(string filename,const int which); //print current positions
  void appendData3d(ofstream &outfile, double x, double y, double z); //append 3d data to file
  void appendData2d(ofstream &outfile, double x, double y); //append 2d data to file
  Print(Molecule *mol);
  ~Print() {}
};

#endif

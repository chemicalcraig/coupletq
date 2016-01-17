#ifndef MAIN_H
#define MAIN_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <list>
#include "util.h"
#include "parse.h"
#include "fret.h"
#include "pert.h"
#include "grid.h"
#include "print.h"
#include "coulomb.h"
#include "reader.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_eigen.h"
using namespace std;

/*****************************
 * Some variables
 * ***************************/

//Input file
ifstream comfile;

//number of molecules, electronic states to consider, atoms, and interaction order
int nmol, nstates, natoms, interactionOrder;

//get number of atoms in a molecule
int getNatoms(string filename, int nmol, Molecule *mol);

//This is the main molecule object
Molecule *mol;

//exciton coupling constant
double coupling;

//Second order interaction matrix in 3-bit eigenbasis
double *intham;
double *int3;

#endif // MAIN_H

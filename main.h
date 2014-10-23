#ifndef MAIN_H
#define MAIN_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "util.h"
#include "parse.h"
#include "fret.h"
#include "pert.h"
#include "grid.h"
using namespace std;

/*****************************
 * Some variables
 * ***************************/

const string str_Atom[] = { " X",
                            " H","He","Li","Be"," B"," C"," N"," O"," F","Ne",
                            "Na","Mg","Al","Si"," P"," S","Cl","Ar"," K","Ca",
                            "Sc","Ti"," V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
                            "Ga","Ge","As","Se","Br","Kr","Rb","Sr"," Y","Zr",
                            "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
                            "Sb","Te"," I","Xe","Cs","Ba","La","Ce","Pr","Nd",
                            "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
                            "Lu","Hf","Ta"," W","Re","Os","Ir","Pt","Au","Hg",
                            "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
                            "Pa"," U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
                            "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds",
                            "Rg","Cn","Uut","Fl","Uup","Lv","Uus","Uuo",
                          };


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
#endif // MAIN_H

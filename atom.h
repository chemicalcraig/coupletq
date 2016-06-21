#ifndef ATOM_H
#define ATOM_H
#include <cstring>
#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

const string str_Atom[] = { " X"
                            "H","He","Li","Be","B","C","N","O","F","Ne",
                            "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
                            "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
                            "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
                            "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
                            "Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
                            "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
                            "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
                            "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
                            "Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
                            "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds",
                            "Rg","Cn","Uut","Fl","Uup","Lv","Uus","Uuo",
                          };
//CTC maybe interface with NWChem?
const double atomicMasses[] = { 0.0,
                            1.01,4,6.94,9.01,10.81,12.01,14.01,16,19,20.18,22.99,24.31,
                            26.98,28.09,30.97,32.07,35.45,39.95,39.1,40.08,44.96,
                            47.88,50.94,52,54.94,55.85,58.93,58.69,63.46,65.39,69.72,
                            72.61,74.92,78.96,79.9,83.8,85.47,87.62,88.91,91.22,92.91,
                            95.94,98,101.07,102.91,106.42,107.87,112.41,114.82,118.71,
                            121.76,127.6,126.9,131.29,132.91,137.33,138.91,140.12,140.91,
                            144.24,145,150.36,151.97,157.25,158.93,162.5,164.93,167.26,
                            168.93,173.04,174.97,178.49,180.95,183.85,186.21,190.2,192.22,
                            195.08,196.97,200.59,204.38,207.2,208.98,209,210,222,223,
                            226.03,227,232.04,231.04,238.03,237.05,244,243,247,247,251,
                            252,257,258,259,260,261,262,263,262,265,266};




class Atom {
  public:
  int num,nao,nocc,nuocc,nbasis;
  int atomicnum;
  int *basisfuncs;
  string type;
  double x,y,z;
  double spos[3]; //spherical coords: r, theta, phi
  double tq,tqindo;
  double charge;
  double *charges;
  double mass;
  double pos[3];
  double ipos[3];
  double tqoverlap;

  void allocateCharges(const int n) {this->charges = new double[n];};
  
  /** Set Mass
   * Overloaded, first one relies on arrays in this file
   * the one that takes the string is from the NWChem output
   */
  void setMass();
  void setMass(string m);

  Atom();
  Atom(int nbf);

};


#endif // ATOM_H

#ifndef READER_H
#define READER_H
/****************************************************
 * The reader class parses an coupletq com file
 * and sets up the appropriate calculation
 *
 * See README for appropriate format and 
 * options
 *
 * **************************************************/

#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <list>
using namespace std;

/************************
 * Data structures to
 * store input settings
 * *********************/
struct ChargeFile {
  int i,f;
  string file;
  string spin;
};

/** Scale a molecule **/
struct Scale {
  double origin[3];
  double scale;
  bool scale_bool;
};

/** Move a molecule **/
struct Move {
  string axis;
  int iaxis;
  double min,max;
  int steps;
};

/** Rotate a molecule **/
struct Rot {
  string axis;
  double theta;
};

/** Initially populated state, for dynamics **/
struct InitPop {
  int mol, state;
  double population;
};

struct Output {
  int mol,state;
  string file;
};

/** Calculation details **/
struct Calc{
/** Calculation Stack **/
  string type,configuration,outfile;
  int itype;
  int istate,fstate;
  double ewindow;
  int molecules;
  bool spin;

  //Configurational stuff
  bool C1_;//=false;
  bool C2_;//=false;
  bool C3_;//=false;
};

struct Mol{
  /** Molecule Stack **/
  int number;
  int nstates;
  int ncharges;
  int nmov;
  int nrot;
  int nscale;
  int target;
  ChargeFile *cf;
  Move *mv;
  Rot *rot;
  Scale *sc;
  string tddftfile;
  bool scaleMol;
  double sf;
};

struct Dyn{
  /** Dynamics Stack **/
  double tstart,tfinish,increment,winc;
  int tsteps,wstep;
  int noutput;
  InitPop *pop;
  Output *out;
  bool popSet,printAll;
};

struct Fret{
};

/*****************************
 * Reader class definition
 * **************************/
class Reader {

  public:
  Calc calc;
  Mol *mol;
  Dyn dyn;
  Fret fret;

  char* filename;

  /** Functions **/
  void readBlock(string s, ifstream &in, int molcount);
  void readSubBlock(string w,string s, ifstream &in,int molcount);
  
  /** Constructor and Destructor **/
  Reader(string f);
  ~Reader() {};
};

#endif

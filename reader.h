#ifndef READER_H
#define READER_H
/****************************************************
 * The reader class parses an exciton com file
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
};

struct Move {
  string axis;
  int iaxis;
  double min,max;
  int steps;
};

struct Rot {
  string axis;
  double theta;
};

struct InitPop {
  int mol, state;
  double population;
};

struct Output {
  int mol,state;
  string file;
};

struct Calc{
/** Calculation Stack **/
  string type;
  int itype;
  int istate,fstate;
  double ewindow;
  int molecules;
};

struct Mol{
  /** Molecule Stack **/
  int number;
  int nstates;
  int ncharges;
  int nmov;
  int nrot;
  ChargeFile *cf;
  Move *mv;
  Rot *rot;
  string tddftfile;
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

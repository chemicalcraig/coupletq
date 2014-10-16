#include "molecule.h"

Molecule::Molecule()
{
}

void Molecule::allocateMem(const int nb) {
  this->nbasis = nb;
  this->nbasisatom = new int[nb];
  this->nbasisatomelements = new string[nb];
  this->nbasisatomorbitals = new string[nb];
  this->overlapm = new double[nb*nb];
  this->occupation = new int[nb];
  this->moeigenv = new double[nb];
  this->mos = new double[nb*nb];
  this->nlindep = 0;
  this->nmo = nb;//= mol->nbasis;
}

void Molecule::allocateMemLindep(const int nb, const int nlindep) {
  //delete[] this->mos;
  this->nlindep = nlindep;
  this->nmo = nb-nlindep;
  //this->mos = new double[this->nmo*nb];
}

void Molecule::allocateMemAtoms(const int na) {
  this->natoms = na;
  this->atoms = new Atom[na];
  this->activeCharges = new double[na];
}

void Molecule::allocateMemTddft() {
  int ntrans = 9*this->nroots;

  this->excenergy = new double[this->nroots];
  this->transmoment = new double[ntrans];
  this->oscstrength = new double[this->nroots];
  //skipping CI vectors for now since we'll be
  //focusing on spatial densities
  //this->ci = new double[this->nroots*this->nmo*this->nmo];
  //for (int i=0; i<this->nroots*this->nmo*this->nmo; i++)
  //  this->ci[i] = 0.;
}

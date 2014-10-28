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

//subroutines to rotate molecule, theta in radians
void Molecule::rotateTheta(double theta, int axis) {
  double sum = 0.;
  double pos[3];
  switch(axis) {
    //rotate about x-axis
    case 0:
      this->rot[0] = 1.;
      this->rot[1+3*1] = cos(theta);
      this->rot[1+3*2] = -1.*sin(theta);
      this->rot[2+3*1] = sin(theta);
      this->rot[2+3*2] = cos(theta);
      break;
    
    //y-axis
    case 1:
      this->rot[1+1*3] = 1.;
      this->rot[0+3*0] = cos(theta);
      this->rot[0+3*2] = sin(theta);
      this->rot[2+3*0] = -1.*sin(theta);
      this->rot[2+3*2] = cos(theta);
      break;

    //z-axis
    case 2:
      this->rot[2+2*3] = 1.;
      this->rot[0+3*0] = cos(theta);
      this->rot[1+3*0] = sin(theta);
      this->rot[0+3*1] = -1.*sin(theta);
      this->rot[1+3*1] = cos(theta);
      break;
  }
  for (int i=0; i<this->natoms; i++) {
    for (int j=0; j<3; j++) {
      sum = 0.;
      for (int k=0; k<3; k++) {
        sum += rot[j+k*3] * this->atoms[i].pos[k];
      }
      pos[j] = sum;
    }
    for (int j=0; j<3; j++)
      this->atoms[i].pos[j] = pos[j];
  }
}

//rotate about vector connecting centers of mass
void Molecule::rotateCom(double theta, double *cm) {

  //Get unit vector connecting the com's
  double diff[3];
  double sum=0.;
  
  //construct vector connecting COM's
  for (int i=0; i<3; i++) {
    diff[i] = this->com[i] - cm[i];
    sum += diff[i]*diff[i];
  }
  
  //normalize it
  for (int i=0; i<3; i++) {
    diff[i] /= sqrt(sum);
  }

  //construct rotation matrix
  this->rotmatcom[0 + 3*0 ] = cos(theta) + diff[0]*diff[0]*(1-cos(theta));
  this->rotmatcom[0 + 3*1 ] = diff[0]*diff[1]*(1-cos(theta)) - diff[2]*sin(theta);
  this->rotmatcom[0 + 3*2 ] = diff[0]*diff[2]*(1-cos(theta)) + diff[1]*sin(theta);
  this->rotmatcom[1 + 3*0 ] = diff[0]*diff[1]*(1-cos(theta)) + diff[2]*sin(theta);
  this->rotmatcom[1 + 3*1 ] = cos(theta) + diff[1]*diff[1]*(1-cos(theta));
  this->rotmatcom[1 + 3*2 ] = diff[1]*diff[2]*(1-cos(theta)) - diff[0]*sin(theta);
  this->rotmatcom[2 + 3*0 ] = diff[2]*diff[0]*(1-cos(theta)) - diff[1]*sin(theta);
  this->rotmatcom[2 + 3*1 ] = diff[1]*diff[2]*(1-cos(theta)) + diff[0]*sin(theta);
  this->rotmatcom[2 + 3*2 ] = cos(theta) + diff[2]*diff[2]*(1-cos(theta));

  //apply matrix to atomic positions
  for (int atms=0; atms<this->natoms; atms++) {
    double pos[3],pos2[3];
    pos[0] = this->atoms[atms].x;
    pos[1] = this->atoms[atms].y;
    pos[2] = this->atoms[atms].z;
    for (int i=0; i<3; i++) {
      double sum = 0.;
      for (int j=0; j<3; j++) {
        sum += this->rotmatcom[i + 3*j] * pos[j];
      }
      pos2[i] = sum;
    }
    this->atoms[atms].pos[0] = pos2[0];
    this->atoms[atms].pos[1] = pos2[1];
    this->atoms[atms].pos[2] = pos2[2];
  }
}

void Molecule::setcom() {
  for (int i=0; i<this->natoms; i++) {
  }
}

/** Translate molecule along cartesian axis **/
void Molecule::translate(const int which, double howmuch) {
  this->com[which] += howmuch;
  switch (which) {
    case 0: //x-axis
      for (int i=0; i<this->natoms; i++) {
          this->atoms[i].pos[0] += howmuch;
        }
      break;
    case 1: //y-axis
      for (int i=0; i<this->natoms; i++) {
          this->atoms[i].pos[1] += howmuch;
        }
      break;
    case 2: //z-axis
      for (int i=0; i<this->natoms; i++) {
        this->atoms[i].pos[2] += howmuch;
      }
      break;
  }
}

void Molecule::resetall() {
  for (int i=0; i<this->natoms; i++) {
    for (int j=0; j<3; j++) {
      this->atoms[i].pos[j] = this->atoms[i].ipos[j];
    }
    this->atoms[i].x = this->atoms[i].pos[0];
    this->atoms[i].y = this->atoms[i].pos[1];
    this->atoms[i].z = this->atoms[i].pos[2];
  }
}

void Molecule::resetExcept(int keep) {
  for (int i=0; i<this->natoms; i++) {
    for (int j=0; j<3; j++) {
      if (j != keep)
        this->atoms[i].pos[j] = this->atoms[i].ipos[j];
    }
  }
}

void Molecule::setPostoInit() {
  for (int i=0; i<this->natoms; i++) {
    for (int j=0; j<3; j++)
      this->atoms[i].ipos[j] = this->atoms[i].pos[j];
  }
}
/*******************************************
 * Operators
 * *****************************************/
/********************************************
 * Copy operator
 */
Molecule Molecule::operator=(const Molecule& other) {
  this->natoms = other.natoms;
  this->  interaction = other.interaction;
  this->  groundenergy = other.groundenergy;
  this->nstates = other.nstates;
  this->nao = other.nao;
  this->nroots = other.nroots;
  this->nocc = other.nocc;
  this->nuocc = other.nuocc;
  this->dipx = other.dipx;
  this->dipy = other.dipy;
  this->dipz = other.dipz;
  this-> nlindep = other.nlindep;
  this->excMethod = other.excMethod;
  this->activeCharges = other.activeCharges;

  //Copy Atoms
  this->atoms = other.atoms;
  //Copy grid
  this->grid = other.grid;
}


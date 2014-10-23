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

//subroutines to rotate molecule
void Molecule::rotateTheta(double theta, int axis) {
  switch(axis) {
    //rotate about x-axis
    case 1:
      this->rotx[0] = 1.;
      this->rotx[1+3*1] = cos(theta);
      this->rotx[1+3*2] = -1.*sin(theta);
      this->rotx[2+3*1] = sin(theta);
      this->rotx[2+3*2] = cos(theta);
    
    //y-axis
    case 2:
      this->roty[1+1*3] = 1.;
      this->roty[0+3*0] = cos(theta);
      this->roty[0+3*2] = sin(theta);
      this->roty[2+3*0] = -1.*sin(theta);
      this->roty[2+3*2] = cos(theta);

    //z-axis
    case 3:
      this->rotz[2+2*3] = 1.;
      this->rotz[0+3*0] = cos(theta);
      this->rotz[1+3*0] = sin(theta);
      this->rotz[0+3*1] = -1.*sin(theta);
      this->rotz[1+3*1] = cos(theta);
  }
}

//rotate about vector connecting centers of mass
void Molecule::rotateCom(double theta, double *cm) {
  ofstream rotout;
  rotout.open("rotcoords");

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
        //cout<<i<<" "<<j<<" "<<this->rotmatcom[i+j*3]<<endl;
        sum += this->rotmatcom[i + 3*j] * pos[j];
      }
      pos2[i] = sum;
    }
    this->atoms[atms].x = pos2[0];
    this->atoms[atms].y = pos2[1];
    this->atoms[atms].z = pos2[2];
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
          this->atoms[i].x += howmuch;
          this->atoms[i].pos[0] += howmuch;
        }
    case 1: //y-axis
      for (int i=0; i<this->natoms; i++) {
          this->atoms[i].y += howmuch;
          this->atoms[i].pos[1] += howmuch;
        }
    case 2: //z-axis
      for (int i=0; i<this->natoms; i++) {
        this->atoms[i].z += howmuch;
        this->atoms[i].pos[2] += howmuch;
      }
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

void Molecule::reset(int keep) {
  for (int i=0; i<this->natoms; i++) {
    for (int j=0; j<3; j++) {
      if (j != keep)
        this->atoms[i].pos[j] = this->atoms[i].ipos[j];
    }
    //this->atoms[i].x = this->atoms[i].pos[0];
    //this->atoms[i].y = this->atoms[i].pos[1];
    //this->atoms[i].z = this->atoms[i].pos[2];
  }
}

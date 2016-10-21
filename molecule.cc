#include "molecule.h"
#include "molutil.h"

/****************************************
 * Set initial geometric configuration
 * **************************************/
void Molecule::setInit(Reader r, int i) {
  int ta,ra;

    /** Translations, sets up grid **/
    for (int j=0; j<r.mol[i].nmov; j++) {
      if (r.mol[i].mv[j].axis.compare(0,1,"x",0,1)==0) {
        ta = 0;
      } else if (r.mol[i].mv[j].axis.compare(0,1,"y",0,1)==0) {
        ta = 1;
      } else if (r.mol[i].mv[j].axis.compare(0,1,"z",0,1)==0) {
        ta = 2;
      }
      
      this->grid[ta].setParams(r.mol[i].mv[j].min,r.mol[i].mv[j].max,r.mol[i].mv[j].steps);
    } //end translations

     /** Rotations, rotates molecule before translations **/
    for (int j=0; j<r.mol[i].nrot; j++) {
      cout<<"rotating molecule "<<i<<endl;
      if (r.mol[i].rot[j].axis.compare(0,1,"x",0,1)==0) {
        ra = 0;
      } else if (r.mol[i].rot[j].axis.compare(0,1,"y",0,1)==0) {
        ra = 1;
      } else if (r.mol[i].rot[j].axis.compare(0,1,"z",0,1)==0) {
        ra = 2;
      }
      this->rotateTheta(r.mol[i].rot[j].theta,ra);
    } //end rotations

    /** Scaling **/
    if (r.mol[i].scaleMol) {
      cout<<"scaling molecule "<<i<<" by "<<r.mol[i].sf<<" X"<<endl;
      this->scaleR(r.mol[i].sf);
    }
}

/****************************************
 * Initialize the Molecule Object
 * *************************************/
Molecule *initialize(Reader r) {

  /** Declare Molecule **/

  Molecule * mol = new Molecule[r.calc.molecules];
  mol[0].outputfilename = r.calc.outfile;

  /** Get number of atoms and densities **/
  for (int i=0; i<r.calc.molecules; i++) {
    mol[i].rot = new double[9];
    mol[i].rotmatcom = new double[9];
    mol[i].nstates = r.mol[i].nstates;
    mol[i].nmol = r.calc.molecules;
    mol[i].natoms = getNatoms(r.mol[i].cf[0].file,r.calc.molecules,mol);
    mol[i].atoms = new Atom[mol[i].natoms];
        cout<<"Molecule "<<i+1<<" has "<<mol[i].natoms<<" atoms"<<endl;
    if (r.mol[i].target != 0) {
      mol[i].target = r.mol[i].target;
      cout<<"looking for tddft target state "<<mol[i].target<<" on mol "<<i<<endl;
    } else {
      mol[i].target = 1;
    }
    if (r.calc.itype==1 || r.calc.itype==5) {
      mol[i].istate = r.mol[i].cf[0].i;
      mol[i].fstate = r.mol[i].cf[0].f;
    }
    //Get spin allowed/forbidden condition
    mol[i].spinAllowed = r.calc.spin;
    
    //Allocate atoms and all of their densities
    //We use lower triangular form for the couplings
    for (int j=0; j<mol[i].natoms; j++) {
      mol[i].atoms[j].allocateCharges(r.mol[i].nstates*r.mol[i].nstates);
    }
    cout<<"Getting transition charges for molecule "<<i+1<<endl;
    /** Get the transition charges **/
    /** This also sets initial positions **/
    for (int j=0; j<r.mol[i].ncharges; j++) {
      if (r.calc.itype != 5) {
        getCharges(r.mol[i].cf[j].file,&mol[i],
                r.mol[i].nstates,r.mol[i].cf[j].i,
                r.mol[i].cf[j].f);
      } else {
        cout<<"Spin-specific charges, "<<r.mol[i].cf[0].spin<<endl;
        getCharges(r.mol[i].cf[j].file,&mol[i],
                r.mol[i].nstates,r.mol[i].cf[j].i,
                r.mol[i].cf[j].f,r.mol[i].cf[0].spin);
       }
    } //end tq retrieval
  
    /** Get Excitation energies from TDDFT calc **/
    getTDDFT(r.mol[i].tddftfile,&mol[i]);
    cout<<"Done retrieving tddft information for molecule "<<i<<endl;
  }

  /** Set Molecular Mass **/
  for (int i=0; i<r.calc.molecules; i++) {
    double mass = 0.;
    for (int j=0; j<mol[i].natoms; j++) {
      mass += mol[i].atoms[j].mass;
    }
    mol[i].setMass(mass);
  }
    
  return mol;
}

/*************************************
  *** Molecule class initializer ***
  ************************************/
Molecule::Molecule()
{
  /** Set some defaults **/
  this->spinAllowed=true;
  this->spinstate = 0;
  this->dip = new double[3];
  this->idip = new double[3];
  this->com = new double[3];
  this->icom = new double[3];
 
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

  this->excenergy = new double[this->nroots+1];
  this->transmoment = new double[ntrans];
  this->oscstrength = new double[this->nroots];
  //skipping CI vectors for now since we'll be
  //focusing on spatial densities
  //this->ci = new double[this->nroots*this->nmo*this->nmo];
  //for (int i=0; i<this->nroots*this->nmo*this->nmo; i++)
  //  this->ci[i] = 0.;
}

/** Set atomic properties **/

void Molecule::setAtomicMasses(string t, string m) {
  for (int i=0; i<this->natoms; i++) {
    if (this->atoms[i].type == t) {
      this->atoms[i].setMass(m);
      continue;
    }
  }
}

/*
 * Subroutines to move and rotate the molecule 
 */

void Molecule::arrangeMol() {
  for (int i=0; i<this->griddim; i++) {
    this->translate(i,this->grid[i].min);
  }
}


//subroutines to rotate molecule, theta in radians
void Molecule::rotateTheta(double theta, int axis) {
  
  double sum = 0.;
  double pos[3],pos2[3],pos3[3];

  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      this->rot[i+j*3] = 0.;

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

  cblas_dgemv(CblasColMajor,CblasNoTrans,
              3,3,1,this->rot,3,this->com,1,0,pos2,1);
  cblas_dgemv(CblasColMajor,CblasNoTrans,
              3,3,1,this->rot,3,this->dip,1,0,pos3,1);
  for (int i=0; i<this->natoms; i++) {
    cblas_dgemv(CblasColMajor,CblasNoTrans,
              3,3,1,this->rot,3,this->atoms[i].pos,1,0,pos,1);
    for (int j=0; j<3; j++) {
      this->atoms[i].pos[j] = pos[j];
    }
  }

  for (int j=0; j<3; j++) {
    this->com[j] = pos2[j];
    this->dip[j] = pos3[j];
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
    double pos[3],pos2[3],pos3[3];
    for (int i=0; i<3; i++) {
      sum = 0.;
      double sum2 = 0.;
      double sum3 = 0.;
      for (int j=0; j<3; j++) {
        sum += this->rotmatcom[i + 3*j] * this->atoms[atms].pos[j];
        sum2 += this->rotmatcom[i+3*j] * this->com[j];
        sum3 += this->rotmatcom[i+3*j] * this->dip[j];
      }
      pos[i] = sum;
      pos2[i] = sum2;
      pos3[i] = sum3;
    }
    for (int i=0; i<3; i++) {
      this->atoms[atms].pos[i] = pos[i];
      this->com[i] = pos2[i];
      this->dip[i] = pos3[i];
    }
  }
}

/** Set molecule's mass **/
void Molecule::setMass(double m) {
  this->mass = m;
}

/** Set center of mass **/
void Molecule::setCom() {
  double sum[3];
  for (int i=0; i<3; i++) {
    sum[i] = 0.;
    for (int j=0; j<natoms; j++) {
      sum[i] += atoms[j].pos[i] * atoms[j].mass;
    }
    sum[i] /= this->mass;
    this->com[i] = sum[i];

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
        //this->com[0] += howmuch;
        //this->dip[0] += howmuch;
      break;
    case 1: //y-axis
      for (int i=0; i<this->natoms; i++) {
          this->atoms[i].pos[1] += howmuch;
        }
        //this->com[1] += howmuch;
        //this->dip[1] += howmuch;
      break;
    case 2: //z-axis
      for (int i=0; i<this->natoms; i++) {
        this->atoms[i].pos[2] += howmuch;
      }
      //this->com[2] += howmuch;
      //this->dip[2] += howmuch;
      break;
  }
  this->setCom();
}

/** Scale molecule by a factor 's' **/
void Molecule::scaleR(double s) {
  for (int i=0; i<this->natoms; i++) {
    this->atoms[i].spos[0] *= s;
    //set cartesian from new r
    //adj adjusts for negative x in transformation
    int adj = 1;
    if (atoms[i].pos[0] < 0)
      adj = -1;

    this->atoms[i].pos[0] = adj*this->atoms[i].spos[0]
                            * sin(this->atoms[i].spos[1])
                            * cos(this->atoms[i].spos[2]);
    this->atoms[i].pos[1] = adj*this->atoms[i].spos[0]
                            * sin(this->atoms[i].spos[1])
                            * sin(this->atoms[i].spos[2]);
    this->atoms[i].pos[2] = this->atoms[i].spos[0]
                            * cos(this->atoms[i].spos[1]);
    this->atoms[i].ipos[0] = this->atoms[i].pos[0]; 
    this->atoms[i].ipos[1] = this->atoms[i].pos[1]; 
    this->atoms[i].ipos[2] = this->atoms[i].pos[2];
  
  }
  this->setCom();
}

void Molecule::resetall() {

  for (int j=0; j<3; j++) {
    for (int i=0; i<this->natoms; i++) {
      this->atoms[i].pos[j] = this->atoms[i].ipos[j];
    }
    //cout<<"position "<<this->icom[j]<<endl;
    this->com[j] = this->icom[j];
    this->dip[j] = this->idip[j];
  }
}

void Molecule::moveTo(const int which, const double where) {
  double howmuch = where - this->com[which];
  for (int i=0; i<this->natoms; i++) {
    this->atoms[i].pos[which] += howmuch;
  }
  this->com[which] += howmuch;
}

void Molecule::resetExcept(int keep) {
  for (int i=0; i<this->natoms; i++) {
    for (int j=0; j<3; j++) {
      if (j != keep) {
        this->atoms[i].pos[j] = this->atoms[i].ipos[j];
        //this->com[j] = this->icom[j];
        this->dip[j] = this->idip[j];
      }
    }
  }
  this->setCom();
}

void Molecule::setPostoInit() {
  for (int i=0; i<this->natoms; i++) {
    for (int j=0; j<3; j++) {
      this->atoms[i].ipos[j] = this->atoms[i].pos[j];
      this->icom[j] = this->com[j];
      this->idip[j] = this->dip[j];
    }
  }
}

/*******************************************
 * Operators
 * *****************************************/
/********************************************
 * Copy operator
 */
Molecule Molecule::operator=(const Molecule& other) {
  this->spinAllowed = other.spinAllowed;
  this->natoms = other.natoms;
  this->interaction = other.interaction;
  this->groundenergy = other.groundenergy;
  this->nstates = other.nstates;
  this->nao = other.nao;
  this->nroots = other.nroots;
  this->nocc = other.nocc;
  this->nuocc = other.nuocc;
  this->dipx = other.dipx;
  this->dipy = other.dipy;
  this->dipz = other.dipz;
  this->activeCharges = other.activeCharges;
  this->dip = other.dip;
  this->idip = other.idip;
  this->dipmag = other.dipmag;
  this->fstate = other.fstate;
  this->istate = other.istate;
  this->nindices = other.nindices;
  this->spinstate = other.spinstate;
  this->interaction = other.interaction;
  this->groundenergy = other.groundenergy;
  this->mass = other.mass;
  this->nstates = other.nstates;
  this->target = other.target;
  this->nlindep = other.nlindep;
  this->excMethod = other.excMethod;
  this->rot = other.rot;
  this->rotmatcom = other.rotmatcom;
  this->griddim = other.griddim;
  this->ci = other.ci;
  this->nbasisatom = other.nbasisatom;
  this->nbasisatomorbitals = other.nbasisatomorbitals;
  this->nbasisatomelements = other.nbasisatomelements;
  this->overlapm = other.overlapm;
  this->excenergy = other.excenergy;
  this->transmoment = other.transmoment;
  this->oscstrength = other.oscstrength;
  this->occupation = other.occupation;
  this->moeigenv = other.moeigenv;
  this->mos = other.mos;
  this->posx = other.posx;
  this->posy = other.posy;
  this->posz = other.posz;
  this->nmol = other.nmol;
  this->nx = other.nx;
  this->ny = other.ny;
  this->nz = other.nz;
  this->com = other.com;
  this->icom = other.icom;
  this->indices = other.indices;
  //Copy Atoms
  this->atoms = other.atoms;
  //Copy grid
  this->grid = other.grid;
  return *this;
}


#ifndef PARSE_H
#define PARSE_H

/**************************************
 * This file contains routines
 * to parse NWChem output files
 * ***********************************/

#include <cstring>
#include <fstream>
#include <iostream>
#include "atom.h"
#include "molecule.h"
#include <iomanip>
using namespace std;

/* Search strings */
string geom_str=" Output coordinates";
string basislabel_str="    Basis function labels";
string overlap_str = " global array: AO ovl";
string nbasis_str="          This is a Direct SCF calculation.";
string tddft_str="  Convergence criterion met";
string tddft_state_stop_str="-----------------------------";
string tddft_task_str1="tddft";
string tddft_task_str2="TDDFT";
string tddft_task_str3="Tddft";
string mo_analysis_str="                       DFT Final Molecular Orbital Analysis";
string lindep_str = " !! The overlap matrix has";
string aobas_str = "          AO basis - number of functions:";
string com_str = " center of mass";
string mass_str = "      Atomic Mass";
/* Number of columns for printing of matrices */
int ncol = 6;

/*** Function to parse an NWChem output file ***/
/* Returns true if successful, false otherwise */
/* 
 * str        - input, log file name
 * mol       - input, Molecule object
 *
 * mol->nroots    -output, number of tddft roots
 * mol->natoms    -output, number of atoms in molecule
 * mol->atom     -output, individual atom objects
 * mol->atom.type  - output, element
 * mol->atom.x,y,z - output, geometry
 * mol->atom.charge  output, charge on atom
 */

bool parseLog(string str, Molecule *mol) {
    
    ifstream infile;
    infile.precision(9);
    infile.open(str.c_str());
    cout<<"Opening file: "<<str.c_str()<<endl;

    char tempc[1000];
    int natoms=0;
    int nbasis=0;
    bool tddftstack = true;
    infile>>ws;
    while(infile.getline(tempc, 1000)) {
      string temps(tempc);
      int current = infile.tellg();
      
/** Get Number of Roots in TDDFT **/
      if ((temps.compare(0,5,tddft_task_str1,0,5) == 0 ||
        temps.compare(0,5,tddft_task_str2,0,5) == 0 ||
        temps.compare(0,5,tddft_task_str3,0,5) == 0) &&
        tddftstack) {
        cout<<"getting nroots "<<temps<<endl;
        while (temps.compare(0,3,"end",0,3) !=0 ) {
          if (temps.compare(0,6,"nroots")==0 ||
            temps.compare(0,6,"NROOTS")==0 ||
            temps.compare(0,6,"Nroots")==0 ||
            temps.compare(0,6,"NRoots")==0) {

            temps = strtok(tempc," ");
            temps = strtok(NULL," ");
            mol->setnroots(atoi(temps.c_str()));
            cout<<mol->nroots<<endl;
            tddftstack = false;
          } //end nroots
          infile>>ws;
          infile.getline(tempc,1000);
          temps = tempc;
        }
      } //end tddft


/** Geometry **/     
      if (temps.compare(0,10,geom_str,0,10) == 0) {
      cout<<"getting geometry"<<endl;
        /* Skip a few lines to get to the geometry */
        getnlines(infile,tempc,3,1000);

        /* Read geometry */
        //get position in file
/*        int pos = infile.tellg();
        while (temps != "Atomic") {
          infile>>ws;
          infile.getline(tempc,1000);
          temps = strtok(tempc," ");
          //get natoms
          if (temps != "Atomic") {
            natoms++;
          }
        }

        /* Make an molecule entry and declare natoms for the molecule */
//        mol->allocateMemAtoms(natoms);

        //return to top of geometry
//        infile.seekg(pos);
        infile.getline(tempc,1000);
        temps = strtok(tempc," ");
        int i=0;
        while(temps!= "Atomic") {
          //Assign atom numbers
          mol->atoms[i].num = atoi(temps.c_str());
          temps = strtok(NULL," ");

          //Get atom type
          mol->atoms[i].type = temps;
          
          //Get atom charge
          temps = strtok(NULL," ");
          mol->atoms[i].charge = atof(temps.c_str());

          //X-coordinate
          temps = strtok(NULL," ");
          mol->atoms[i].x = atof(temps.c_str());
          
          //Y-coordinate
          temps = strtok(NULL," ");
          mol->atoms[i].y = atof(temps.c_str());

          //Z-coordinate
          temps = strtok(NULL," ");
          mol->atoms[i].z = atof(temps.c_str());

          infile>>ws;
          infile.getline(tempc,1000);
          temps = strtok(tempc," ");
          i++;
        }
        cout<<"Temps "<<temps<<endl;
      } //end geometry

/** Number of AO basis functions **/
      
      if (temps.compare(0,40,nbasis_str,0,40)==0) {
        cout<<"getting nAO"<<endl;
        infile.getline(tempc,1000);
        temps = strtok(tempc,":");
        temps = strtok(NULL,":");

        //allocate memory for molecule's properties
        mol->allocateMem(atoi(temps.c_str()));
      } //end ao basis

/** Number linearly dependent vectors **/
    if (temps.compare(0,25,lindep_str,0,25)==0) {
      cout<<"getting nlindep vectors"<<endl;
      temps = strtok(tempc," ");
      for (int i=0; i<5; i++) temps = strtok(NULL," ");
      mol->allocateMemLindep(mol->nbasis,atoi(temps.c_str()));
    }


/** Basis Function Labels **/
      
      if (temps.compare(0,20,basislabel_str,0,20) == 0) {
      cout<<"getting basis function labels"<<endl;
        getnlines(infile,tempc,4,1000);
        temps = tempc;
        
        for (int i=0; i<mol->natoms; i++) 
          mol->atoms[i].nao = 0;

        while (temps.compare(0,10,overlap_str,1,10) != 0) {
          temps = tempc;
          temps = strtok(tempc," ");

          //get function number
          int num = atoi(temps.c_str());

          //get atom number
          temps = strtok(NULL," ");
          mol->nbasisatom[num-1] = atoi(temps.c_str());

          //increment number of AOs on atom
          mol->atoms[atoi(temps.c_str())-1].nao ++;
          
          //get atom element
          temps = strtok(NULL," ");
          mol->nbasisatomelements[num-1] = temps;

          //get type of AO
          temps = strtok(NULL," ");
          mol->nbasisatomorbitals[num-1] = temps;
          
          infile>>ws;
          infile.getline(tempc,1000);
          temps = tempc;
        } //end import basis functions

      } //end basis label

/** Basis function overlap matrix **/
      if (temps.compare(0,20,overlap_str,1,20) == 0) {
        cout<<"getting overlap matrix"<<endl;
        int nblocks = mol->nbasis/ncol;

        for (int block=0; block<nblocks; block++) {
          if (block == 0)
            getnlines(infile,tempc,4,1000);
          else
            getnlines(infile,tempc,3,1000);
          for (int i=0; i<mol->nbasis; i++) {
            temps = strtok(tempc," ");
            for (int j=0; j<ncol; j++) {
              temps = strtok(NULL," ");
              mol->overlapm[i+(j+ncol*block)*mol->nbasis]
                = atof(temps.c_str());
                //if (i==189 && j==5) cout<<mol->overlapm[i+(j+ncol*block)*mol->nbasis]<<endl;
            } // end col
              infile.getline(tempc,1000);
          }  // end rows
        } //end full blocks
        int remain = mol->nbasis%ncol;
        cout<<"remain = "<<remain<<endl;
        if (remain > 0)
          getnlines(infile,tempc,3,1000);
        for (int i=0; i<mol->nbasis; i++) {
          temps = strtok(tempc," ");
          for (int j=0; j<remain; j++) {
            temps = strtok(NULL," ");
            mol->overlapm[i+(j+nblocks*ncol)*mol->nbasis] 
              = atof(temps.c_str());
            if (i==189 && j==remain-1) cout<<mol->overlapm[i+(j+nblocks*ncol)*mol->nbasis]<<endl;
          } //end col
          infile.getline(tempc,1000);
        } //end row
      } //end overlap matrix

/** Final MO analysis **/
    if (temps.compare(0,40,mo_analysis_str,0,40) == 0) {
      cout<<"getting MO's, nmo= "<<mol->nmo<<endl;
      int current;
      infile.getline(tempc,1000);
      infile>>ws;
      infile.getline(tempc,1000);
      temps = tempc;

      for (int i=0; i<mol->nmo; i++) {

        temps = strtok(tempc," ");
        temps = strtok(NULL," ");

        //which MO
        int nmo = atoi(temps.c_str())-1;

        temps = strtok(NULL,"=");
        temps = strtok(NULL,"=");
        
        //occupation number
        mol->occupation[nmo] = (int)atof(temps.c_str());

        getnlines(infile,tempc,3,1000);
        
        bool first=true;
        //vectors
        while (temps != "Vector" && temps != "center"
            && temps != "---------------------------------------------") {

          if (first) {
            infile>>temps;
            first = false;
          }
          int nbas = atoi(temps.c_str())-1;
          infile>>temps;
          mol->mos[nmo+nbas*mol->nmo] = atof(temps.c_str());
          current = infile.tellg();
          for (int j=0; j<3; j++) infile>>temps;
          if (temps == "d") {
            infile>>temps;
          }
          infile>>temps;
          if (temps != "Vector" && temps != "center" 
            && temps != "---------------------------------------------") { 
              nbas = atoi(temps.c_str())-1;
              infile>>temps;
              mol->mos[nmo+nbas*mol->nmo] = atof(temps.c_str());
              //if (nmo==0) cout<<nbas<<" "<<mol->mos[nmo+nbas*mol->nmo]<<" "<<temps<<endl;
              for (int j=0; j<3; j++) infile>>temps;
              if (temps == "d") {
                infile>>temps;
              }
              infile>>temps;
            }
      }//end vector coeffs
      
      infile.seekg(current);
      infile.getline(tempc,1000);
      infile.getline(tempc,1000);
      infile.getline(tempc,1000);
      }
    }//end MOs

    
/** TDDFT excited states **/
    if (temps.compare(0,10,tddft_str,0,10) == 0) {
      cout<<"getting CI vectors"<<endl;
      mol->allocateMemTddft();
      
      getnlines(infile,tempc,4,1000);

      for (int root=0; root<mol->nroots; root++) {

        infile.getline(tempc,1000);
        temps = strtok(tempc," ");
        for (int i=0; i<4; i++) temps = strtok(NULL," ");
        mol->excenergy[root+1] = atof(temps.c_str());

        getnlines(infile,tempc,2,1000);
      
        for (int k=0; k<3; k++) {
          temps = strtok(tempc," ");
          temps = strtok(NULL," ");
      
          for (int j=0; j<3; j++) {
            for (int i=0; i<2; i++) temps = strtok(NULL," ");
              mol->transmoment[k+j*3+root*9] = atof(temps.c_str());
          }

          infile.getline(tempc,1000);
        }//end trans moment

        //get oscillator strength
        temps = strtok(tempc," ");
        for (int i=0; i<3; i++) temps=strtok(NULL," ");
        mol->oscstrength[root] = atof(temps.c_str());
        infile.getline(tempc,1000);
      
        infile.getline(tempc,1000);
        temps = tempc;

        
        //get CI coeffs
        //nrpa is for RPA
        int nrpa = 1;     
        double ycoeff;
        double sumci = 0.;
        double sumci2 = 0.;
        while (temps.compare(0,10,tddft_state_stop_str,0,10) != 0 &&
              temps.compare(0,10,"Target root",0,10) !=0 ) {
          temps = strtok(tempc," ");
          temps = strtok(NULL," ");
          int row = atoi(temps.c_str())-1;
          for (int i=0; i<4; i++) temps = strtok(NULL," .");
          int col = atoi(temps.c_str())-1;
          for (int i=0; i<2; i++) temps = strtok(NULL," ");
          
          //if (nrpa%2 == 0 || mol->excMethod == "cis" || mol->excMethod == "CIS") {
            mol->ci[row + col*mol->nmo + root*mol->nmo*mol->nmo] += atof(temps.c_str());
            sumci += atof(temps.c_str())*atof(temps.c_str());
          //}
          //CTC check this and make compatible with CIS 9-12-14
          /*if (nrpa%2 == 1 ) {
            mol->ci[row + col*mol->nmo + root*mol->nmo*mol->nmo] += atof(temps.c_str());
            sumci2 += atof(temps.c_str())*atof(temps.c_str());
          }
         */ nrpa++;
          infile>>ws;
          infile.getline(tempc,1000);
          temps=tempc;
        } //end CI coeffs
//        cout<<"sum ci "<<sumci<<" "<<sumci2<<endl;
//        for (int i=0; i<mol->nmo; i++)
//          for (int j=0; j<mol->nmo; j++)
//            mol->ci[i+j*mol->nmo+root*mol->nmo*mol->nmo] /= (sumci-sumci2);
      } //end nroots
    }//end tddft
  } //end reading logfile
    return true;
}


/** Get TDDFT information from logfile **/
bool getTDDFT(string str, Molecule *mol) {

    ifstream infile;
    infile.precision(9);
    infile.open(str.c_str());

    char tempc[1000];
    int natoms=0;
    int nbasis=0;
    bool tddftstack = true;

    while(infile.getline(tempc, 1000)) {
      string temps(tempc);
     
    /** Get Stuff from TDDFT directive **/
      if ((temps.compare(0,5,tddft_task_str1,0,5) == 0 ||
        temps.compare(0,5,tddft_task_str2,0,5) == 0 ||
        temps.compare(0,5,tddft_task_str3,0,5) == 0) &&
        tddftstack) {
        
        while (temps.compare(0,3,"end",0,3) !=0 ) {
          //get nroots
          if (temps.compare(0,6,"nroots")==0 ||
            temps.compare(0,6,"NROOTS")==0 ||
            temps.compare(0,6,"Nroots")==0 ||
            temps.compare(0,6,"NRoots")==0) {
          
            temps = strtok(tempc," ");
            temps = strtok(NULL," ");
            mol->setnroots(atoi(temps.c_str()));
            tddftstack = false;
          } //end nroots

          //get spin state
          if ((temps.compare(0,9,"notriplet")==0) || 
              (temps.compare(0,9,"Notriplet")==0) ||
              (temps.compare(0,9,"NOTRIPLET")==0) ||
              (temps.compare(0,9,"NoTriplet")==0)) {
            mol->spinstate = 0; 
            cout<<"Spin state = "<<mol->spinstate<<endl;
          }

          if ((temps.compare(0,9,"nosinglet")==0) || 
              (temps.compare(0,9,"Nosinglet")==0) ||
              (temps.compare(0,9,"NOSINGLET")==0) ||
              (temps.compare(0,9,"NoSinglet")==0)) {
            mol->spinstate = 1;
            cout<<"Spin state = "<<mol->spinstate<<endl;
          }

          infile>>ws;
          infile.getline(tempc,1000);
          temps = tempc;
        }
      } //end tddft roots
     

      //Grab COM while we're here
      if (temps.compare(0,15,com_str,0,15) == 0) {
        getnlines(infile,tempc,2,1000);
        temps = strtok(tempc,"=");
        temps = strtok(NULL," ");
        mol->com[0] = atof(temps.c_str());
        mol->icom[0] = mol->com[0];
        for (int i=0; i<3; i++) temps = strtok(NULL," ");
        mol->com[1] = atof(temps.c_str());
        mol->icom[1] = mol->com[1];
        for (int i=0; i<3; i++) temps = strtok(NULL," ");
        mol->com[2] = atof(temps.c_str());
        mol->icom[2] = mol->com[2];
      }

      //Get atomic masses while we're here
      if (temps.compare(0,17,mass_str,0,17) == 0) {
        
        infile>>ws;
        getnlines(infile,tempc,3,1000);

        while(temps.compare(0,10,"Effective nuclear",0,10) != 0) {
          string temps2 = strtok(tempc," ");
          temps = strtok(NULL," ");
          mol->setAtomicMasses(temps2,temps);
          infile>>ws;
          infile.getline(tempc,1000);
          temps = tempc;
        }
      }
      
      //Get TDDFT information
      if (temps.compare(0,10,tddft_str,0,10) == 0) {

        mol->allocateMemTddft();
        
        //Get the ground state energy
        getnlines(infile,tempc,2,1000);
        temps = strtok(tempc," ");
        for (int i=0; i<3; i++) temps = strtok(NULL," ");

        mol->groundenergy = atof(temps.c_str());
        getnlines(infile,tempc,2,1000);

        for (int root=0; root<mol->nroots; root++) {

          infile.getline(tempc,1000);
          temps = strtok(tempc," ");
          for (int i=0; i<4; i++) temps = strtok(NULL," ");
        
          //get energy of root
          mol->excenergy[root+1] = atof(temps.c_str());

          getnlines(infile,tempc,2,1000);
          
          //skip for triplet
        if (mol->spinstate == 0) {
          for (int k=0; k<3; k++) {
            temps = strtok(tempc," ");
            temps = strtok(NULL," ");
      
            for (int j=0; j<3; j++) {
              for (int i=0; i<2; i++) temps = strtok(NULL," ");
                mol->transmoment[k+j*3+root*9] = atof(temps.c_str());
            }

            infile.getline(tempc,1000);
          }//end trans moment

          //get oscillator strength
          temps = strtok(tempc," ");
          for (int i=0; i<3; i++) temps=strtok(NULL," ");
          mol->oscstrength[root] = atof(temps.c_str());
          infile.getline(tempc,1000);
        } else {
          infile.getline(tempc,1000);
          infile.getline(tempc,1000);
        }

        infile.getline(tempc,1000);
        temps = tempc;
        
        //get CI coeffs
        //nrpa is for RPA
        int nrpa = 1;     
        double ycoeff;
        double sumci = 0.;
        double sumci2 = 0.;
        while (temps.compare(0,10,tddft_state_stop_str,0,10) != 0 &&
              temps.compare(0,10,"Target root",0,10) !=0 ) {
          temps = strtok(tempc," ");
          temps = strtok(NULL," ");
          int row = atoi(temps.c_str())-1;
          for (int i=0; i<4; i++) temps = strtok(NULL," .");
          int col = atoi(temps.c_str())-1;
          for (int i=0; i<2; i++) temps = strtok(NULL," ");
          
          nrpa++;
          infile>>ws;
          infile.getline(tempc,1000);
          temps=tempc;
        } //end CI coeffs
   
        } //end nroots
    }//end tddft
  //return true;
  }
}

/*** Read in MO vectors from movec file ***/
/* str             - input, file name to read
 * mol            - input, Molecule object
 * mol->occupation - output, occupation numbers
 * mol->moeigenv   - output, MO eigenvalues
 * mol->mos        - output, MO in orthonormal basis
 */

bool readMOs(string str, Molecule *mol) {
  cout.precision(10);
  char tempc[1000];
  ifstream infile;
  infile.open(str.c_str());

  if (!infile.is_open()) {
    cout<<"Error opening MO file"<<endl;
    return false;
  } else {
    cout<<"Opening file: "<<str<<endl;
  }
  
  /* Skip right to the occupation numbers */
  getnlines(infile,tempc,14,500);
  mol->nocc = 0;
  mol->nuocc = 0;
  for (int i=0; i<mol->nmo; i++) {
    infile>>tempc;
    mol->occupation[i] = (int)atof(tempc);
    if (mol->occupation[i] == 0) {
      mol->nuocc++;
    } else {
      mol->nocc++;
    }
  }
cout<<"HOMO = "<<mol->nocc<<", # virt. orb. = "<<mol->nuocc<<endl;
  /* Read in eigenvalues */
  for (int i=0; i<mol->nmo; i++) {
    infile>>tempc;
    mol->moeigenv[i] = atof(tempc);
  }

  /* Read in MO vectors */
 /* for (int i=0; i<mol->nmo; i++)
    for (int j=0; j<mol->nbasis; j++) {
      infile>>tempc;
      mol->mos[i+j*mol->nmo] = atof(tempc);
    }
*/
  return true;
}

Molecule *parseComfile(ifstream &comfile) {
  
  //parse input file
  char tempc[1000];
  string temps;
  comfile.getline(tempc,1000);
  temps = strtok(tempc,":");
  temps = strtok(NULL,": ");
  
  //Get number of molecules
  int nmol = atoi(temps.c_str());
  Molecule *mol = new Molecule[nmol];
  for (int i=0; i<nmol; i++)
    mol[i].nmol = nmol;

  //Get order of interaction
  //limited to {1,2}
  //1 = FRET
  //2 = 2nd order perturbation theory
  comfile.getline(tempc,1000);
  temps = strtok(tempc,":");
  temps = strtok(NULL,": ");
  for (int i=0; i<nmol; i++)
    mol[i].interaction = atoi(temps.c_str());

  //Number of electronic states to consider
  //on each molecule
  for (int i=0; i<nmol; i++) {
    comfile.getline(tempc,1000);
    temps = strtok(tempc,":");
    temps = strtok(NULL,": ");
    mol[i].nstates = atoi(temps.c_str());
    cout<<"There are "<<mol[i].nstates<<" excited states for molecule "<<i+1<<endl;
    
    //get natoms and charge densities of each molecule, 
    //  each one requires a separate file generated
    //  using either denstq or aodens
    //this should change if inter-excited transitions are to 
    //be included
    for (int icharges=0; icharges<2*mol[i].nstates-1; icharges++) {
      comfile.getline(tempc,1000);
      temps = strtok(tempc,":");
      temps = strtok(NULL,": ");

      //Natoms
      if (icharges == 0) {
        mol[i].natoms = getNatoms(temps,nmol,mol);
        mol[i].atoms = new Atom[mol[i].natoms];
        cout<<"Molecule "<<i+1<<" has "<<mol[i].natoms<<" atoms"<<endl;
        //Allocate atoms and all of their densities
        //We use lower triangular form for the couplings
        for (int j=0; j<mol[i].natoms; j++) {
          //this should change if inter-excited transitions are to 
          //be included
          mol[i].atoms[j].allocateCharges(mol[i].nstates*mol[i].nstates);
        }
      }
      
      //Get the transition charges.
      //The diagonal components (j==k) are the
      //state densities, while the off diagonal
      //components are the transition densities.
      //Since we are only considering transitions to/from
      //the ground state only the first column become 
      //populated with transition densities
      //We only need the lower triangle since the 

      getCharges(temps,nmol,&mol[i],icharges,mol[i].nstates);
    }

    //symmetrize charges
    for (int iatom = 0; iatom<mol[i].natoms; iatom++) {
      for (int icharge = 0; icharge<mol[i].nstates; icharge++) {
        for (int icharge2 = icharge; icharge2<mol[i].nstates; icharge2++) {
          mol[i].atoms[iatom].charges[icharge2 + icharge*mol[i].nstates] = 
              mol[i].atoms[iatom].charges[icharge + icharge2*mol[i].nstates];
        }
      }
    }
  }
  
  //get energies and dipoles of each state of each molecule
  //energies and dipoles are gotten from a tddft calculation.
  //Energies are all relative to the ground state, and we 
  //take that to be zero for all molecules.
  
  for (int i=0; i<nmol; i++) {
    comfile.getline(tempc,1000);
    temps = strtok(tempc,":");
    temps = strtok(NULL,": ");
    
    //retrieves nroots and respective energies and
    //transition dipoles for all 
    //excited states in given log file

    getTDDFT(temps,&mol[i]);

    cout<<"The transition dipole vector for S1 of molecule "<<i+1<<
    " is ("<<mol[i].transmoment[0+3*0]<<" x, "<<mol[i].transmoment[0+3*1]<<" y "<<
    mol[i].transmoment[0+3*2]<<" z) "<<endl;
    cout<<"The ground state energy of molecule "<<i+1<<" is "<<mol[i].groundenergy<<" a.u."<<endl;
  }

  //get output file name to which coupling will be written
  comfile.getline(tempc,1000);
  temps = strtok(tempc,":");
  temps = strtok(NULL,": ");
  mol[0].outputfilename = temps;
  cout<<"Couplings are written to : "<<mol[0].outputfilename<<endl;

  /** Set Molecular Mass **/
  for (int i=0; i<nmol; i++) {
    double mass = 0.;
    for (int j=0; j<mol[i].natoms; j++) {
      mass += mol[i].atoms[j].mass;
    }
    mol[i].setMass(mass);
  }
  return mol;
}

#endif // PARSE_H

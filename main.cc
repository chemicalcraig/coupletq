/*************************************************************
 * Programs to calculate the three-body
 * coupling elements and energy transfer
 * rates in a perturbative framework using
 * transition charges.
 *
 * Input: transition charges of each molecule involved,
 *        calculated using either AO's or spatial densities
 * 
 * CTC 10-10-14
 * **********************************************************/

#include "main.h"

int main(int argc, char **argv) {
  
  /** echo input file **/
  ifstream efile;
  efile.open(argv[1]);
  char tmpc[1000];
  
  cout<<"*** Start Echo Input file ***"<<endl;
  
  while (!efile.eof()) {
    efile.getline(tmpc,1000);
    cout<<tmpc<<endl;
  }
  efile.close();
  
  cout<<"*** End Echo Input File ***"<<endl;
  
  /** Read input com file **/
  cout<<"*** Parsing input file ***"<<endl;
  Reader read(argv[1]);

  /** Initialize molecules **/
  cout<<"*** Initializing molecules ***"<<endl;
  mol = initialize(read);

  /** Set COM before initial translation**/
  if (read.calc.itype != 5) {
    for (int i=0; i<read.calc.molecules; i++) {
      mol[i].setCom();
    }
  }
 
 /** set up printer and output files **/
  cout<<"*** Setting up output files ***"<<endl;
  Print print(mol);
  ofstream outfile2,pdafile;
  outfile2.open(mol[0].outputfilename.c_str());
  outfile2.precision(16);

  /** Set up the spatial grid data for translating monomers **/
  cout<<"*** Setting up translation grid ***"<<endl;
  int griddim = 3;
  for (int i=0; i<read.calc.molecules; i++) {
    mol[i].grid = new Grid[griddim];
    mol[i].griddim = griddim;
    for (int j=0; j<griddim; j++) {
      mol[i].grid[j].size = griddim;
    }
  }
  
  /** Set up initial conformations **/
  for (int i=0; i<read.calc.molecules; i++) {
    mol[i].setInit(read,i);
    /** calculate molecular diple moment from transition charges
     * unless pda (calc type=4)
     */
    if (read.calc.itype != 4)
      calcdip(mol[i]);
  }
  
/*****************  Setting up Molecular distribution ******************/
  //set initial positions
  arrangeMol(mol);
 
  /** Calculate transition dipole from charges **/
  /** This also sets the transition vector elements **/
  /** This also resets the current conformation to be the initial **/
  for (int i=0; i<read.calc.molecules; i++) {
    mol[i].setCom();
    mol[i].setPostoInit();
  }

  /** Print initial positions in xyz format **/
  ofstream posout;
  posout.open("initpos.xyz");
  int totalAtoms = 0;
  for (int i=0; i<read.calc.molecules; i++) {
    totalAtoms += mol[i].natoms;
  }
  if (read.calc.itype != 5) {
    cout<<"*** Printing initial configuration to: initpos.xyz ***"<<endl;
    posout<<totalAtoms<<endl;
    posout<<"Initial configuration"<<endl;
    for (int i=0; i<read.calc.molecules; i++) {
      for (int j=0; j<mol[i].natoms; j++) {
        posout<<mol[i].atoms[j].type<<" ";
        for (int k=0; k<3; k++) {
          posout<<mol[i].atoms[j].pos[k]<<" ";
        }
        posout<<endl;
      }
    }
  }
  posout.close();

  /** Set indicies matrix (for SSSF) **/
  
  int nindex=1;
  for (int i=0; i<read.calc.molecules; i++) {
    nindex *= mol[i].nstates;
  }
  int3 = new double[nindex*nindex];
  intham = new double[nindex*nindex];
  mol[0].nindices = nindex;
  if (read.calc.itype != 5) {
    setIndices(mol,mol[0].nmol,nindex);
    cout<<"** State indices **"<<endl;
  }
/*  for (int i=0; i<nindex; i++) 
    for (int j=0; j<mol[0].nmol; j++) 
      cout<<"ket "<<i<<" mol "<<j<<" = "<<mol[0].indices[j+i*mol[0].nmol]<<endl;
*/


  /** Create Coulomb Matrix **/
  double temp[nindex*nindex],temp2[nindex*nindex];
  Coulomb coul(nindex);

  /** Create state energy matrix **/
  double energies[nindex];
  double vec[nindex],vec2[nindex];

  /** only do this if not fret, pda, or projection **/
  if ((read.calc.itype != 1) && (read.calc.itype != 4)
      && (read.calc.itype != 5)) {
    for (int i=0; i<nindex; i++) {
      energies[i] = 0.;
      for (int m=0; m<mol[0].nmol; m++) {
        energies[i] += mol[m].excenergy[mol[0].indices[m+i*mol[0].nmol]];
      }
    energies[0] = 0.;
    cout<<"energies "<<i<<" "<<energies[i]<<endl;
  }

  /** Create unfiltered Coulomb matrix **/
  if (read.calc.itype != 5) {
    createCoulomb3(mol,coul);
  /** Filter Coulomb Matrix for energy conservation **/
  for (int i=0; i<nindex; i++) {
    //coul.int3[i+i*nindex] += energies[i];
    for (int j=0; j<nindex; j++) {
      int3[i+j*nindex] = coul.int3[i+j*nindex];
      int3[i+j*nindex] *= window(energies[i],energies[j],read.calc.ewindow,0);
    }
    int3[i+i*nindex] += energies[i];
  }
  
  /** Get eigensystem of Coulomb matrix **/
  gsl_matrix_view m = gsl_matrix_view_array(int3,nindex,nindex);
  gsl_vector *eval = gsl_vector_alloc(nindex);
  gsl_matrix *evec = gsl_matrix_alloc(nindex,nindex);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(nindex);
  gsl_eigen_symmv(&m.matrix,eval,evec,w);
  gsl_eigen_symmv_free(w);
  

  for (int i=0; i<nindex; i++) {vec[i] = gsl_matrix_get(evec,i,0);}
  double sum = 0.;
  for (int i=0; i<nindex; i++) {
    sum = 0.;
    if (read.calc.itype != 4)
      cout<<"evals "<<i<<" "<<gsl_vector_get(eval,i)<<endl;
    coul.evals3[i] = gsl_vector_get(eval,i);
    //cout<<"evals "<<i<<" "<<coul.evals3[i]<<endl;
    for (int j=0; j<nindex; j++) {
      sum += coul.int3[i+j*nindex]*vec[j];
      coul.evecs3[i+j*nindex] = gsl_matrix_get(evec,i,j);
      if (read.calc.itype != 4)
        cout<<"evecs "<<i<<" "<<j<<" "<<coul.evecs3[i+j*nindex]<<endl;
    }
  }
  }
  } //end (not projection calculation)
/*******************  Done Setting up molecules *****************************/

  /************************************
   * Do the appropriate calulation
   */

  double angle = 0.;
  double slip;
  double coupling2;

  switch (read.calc.itype) {

    /**************************************
     *    Transition charge projection
     **************************************/
    case 5:
    {
      double sum = 0.;
      cout<<"Performing a transition charge projection now."<<endl;
      cout<<"Projection written to "<<mol[0].outputfilename<<endl;
      
      /** Scale tq's **/
      scaletq(mol);

      /** Absolute difference in tq by atom **/
      double *p = projecttq(mol);
      for (int i=0; i<mol[0].natoms; i++) {
        sum += p[i];
      }
      sum /= mol[0].natoms;
      cout<<"total sum = "<<sum<<", "<<1.-sum*mol[0].natoms<<
        " "<<sum*mol[0].natoms<<endl;
    
      /** Distance-independent Coulomb interaction **/
      double interaction = getCoulombNoDist(mol,0,1,0,1,0,1);

      cout<<"interaction = "<<interaction*27211<<" meV"<<endl;

      /** Write the projection charges to file **/
      for (int i=0; i<mol[0].natoms; i++) {
        print.appendData2d(outfile2,i,p[i]);
      }

    }

    break;
    
    
    /**************************************
     *    FRET Calculation with PDA
     **************************************/
    
    case 4:
      
      cout<<"Performing a FRET calculation using PDA now."<<endl;
      cout<<"Coupling output written to "<<mol[0].outputfilename<<endl;
      
      for (int r1=0; r1<read.mol[1].mv[0].steps; r1++) {
        pdaCalc(mol,coupling);
        if (r1==0)
          cout<<"First coupling = "<<coupling*27211.396<<" meV"<<endl;
          
        /** write the coupling to file **/
        print.appendData2d(outfile2,
          mol[1].grid[read.mol[1].mv[0].iaxis].min+r1*mol[1].grid[read.mol[1].mv[0].iaxis].dgrid,
          coupling*27211.396);
        
        //translate acceptor
        mol[1].translate(read.mol[1].mv[0].iaxis,mol[1].grid[read.mol[1].mv[0].iaxis].dgrid);
        mol[1].setCom();
      }
      break;
    
    
    /**************************************
     *    FRET Calculation using tq's
     **************************************/
    
    case 1:

      cout<<"Performing a FRET calculation now."<<endl;
      cout<<"Coupling output written to "<<mol[0].outputfilename<<endl;
      for (int r1=0; r1<read.mol[1].mv[0].steps; r1++) {
        fretCalc(mol,coupling);
        if (r1==0)
          cout<<"First coupling = "<<coupling*27211.396<<" meV"<<endl;
        /** write the coupling to file **/
        print.appendData2d(outfile2,
            mol[1].grid[read.mol[1].mv[0].iaxis].min+r1*mol[1].grid[read.mol[1].mv[0].iaxis].dgrid ,//* 4.5,
            coupling*27211.396);
        
        //translate acceptor
        mol[1].translate(read.mol[1].mv[0].iaxis,mol[1].grid[read.mol[1].mv[0].iaxis].dgrid);
        mol[1].setCom();
      }
      break;
    
    
    /**************************************
     *    Dynamics Calculation
     **************************************/
    case 2:
      
      double dum;
      for (int zi=0; zi<1; zi++) {
        pertCalcEigen(mol,coul,energies,int3,intham);
        propagateTime(mol,coul,energies,read.dyn.tstart,
                      read.dyn.tfinish,read.dyn.increment,intham,read);

      }
      break;

      
    /**************************************
     *    SSSF Calculation
     *
     * NB This is set up for displacement 
     * in only one direction at the moment
     **************************************/
    case 3:

      bool closest = true;
      ofstream cfile,cmfile,crossfile1;
      ofstream crossfile2,crossfile3;
      ofstream crossfile4,crossfile5,crossfile6;
      cfile.open("coupling-3.dat");
      if (read.calc.C2_) {
        crossfile1.open("rda2-c2.dat");
        crossfile4.open("ra1a2-c2.dat");
      } else if (read.calc.C1_) {
        crossfile3.open("ra1a2-c1.dat");
        crossfile4.open("rda1da2-c1.dat");
      }
      crossfile2.open("rda1-c.dat");
      
      double r12,r23;
      double minsep = read.mol[2].mv[0].min - read.mol[1].mv[0].min;
      
      /** NB for C3 the number of steps must be the same for all
       * directions. The directions must also be listed in the same
       * order for molecules 1 and 2 in the input file (maybe, haven't
       * tested it but better safe than sorry **/
      //output matrix
      double sssfCoupling[read.mol[1].mv[0].steps*read.mol[2].mv[0].steps];
      double pos1[read.mol[1].mv[0].steps];
      double pos2[read.mol[2].mv[0].steps];
      double posMol1,posMol2,m2start;

  Molecule *molc;// = new Molecule[3];
  Coulomb coulc(nindex);
//#pragma omp parallel private(int3,posMol1,posMol2,m2start,molc,coulc)
{
  
  int3 = new double[nindex*nindex];
  molc = new Molecule[3];
  coulc.reinitialize(nindex);
  double *vecc = new double[nindex];

  for (int i=0; i<3; i++) {
    molc[i] = mol[i];
  }


//#pragma omp for ordered
    for (int r1=0; r1<read.mol[1].mv[0].steps; r1++) {

        posMol1 = molc[1].icom[read.mol[1].mv[0].iaxis] + r1*molc[1].grid[read.mol[1].mv[0].iaxis].dgrid;
        molc[1].moveTo(read.mol[1].mv[0].iaxis,posMol1);
        molc[1].setCom();


      if (!read.calc.C3_) {
        //molc[1].translate(read.mol[1].mv[0].iaxis,mol[1].grid[read.mol[1].mv[0].iaxis].dgrid);
        //molc[1].setCom();
        molc[2].resetall();
        if (read.mol[2].mv[0].min >= 0) {
          molc[2].moveTo(read.mol[2].mv[0].iaxis,mol[1].com[read.mol[2].mv[0].iaxis]+minsep);
        } else {
          molc[2].resetall();
        }
      }
      m2start=molc[2].com[read.mol[2].mv[0].iaxis];
      for (int r2=0; r2<read.mol[2].mv[0].steps; r2++) {
        if (read.calc.C1_) {
          posMol2 = molc[1].icom[read.mol[1].mv[0].iaxis] 
                  + r1*molc[1].grid[read.mol[1].mv[0].iaxis].dgrid
                  + minsep;
        
        
          posMol2 += r2*molc[2].grid[read.mol[2].mv[0].iaxis].dgrid;
        } else {
          posMol2 = molc[2].icom[read.mol[2].mv[0].iaxis]+r2*
                      molc[2].grid[read.mol[2].mv[0].iaxis].dgrid;
        }
        molc[2].moveTo(read.mol[2].mv[0].iaxis,posMol2);
        molc[2].setCom();
        createCoulomb3(molc,coulc);

      /** Filter Coulomb Matrix for energy conservation **/
        for (int i=0; i<nindex; i++) {
          for (int j=0; j<nindex; j++) {

            int3[i+j*nindex] = coulc.int3[i+j*nindex];


            int3[i+j*nindex] *= window(energies[i],energies[j],read.calc.ewindow,0);
          }
          int3[i+i*nindex] += energies[i];
        }


      /** Get eigensystem of Coulomb matrix + bare energies (H+W1 = int3) **/
        gsl_matrix_view m = gsl_matrix_view_array(int3,nindex,nindex);
        gsl_vector *eval = gsl_vector_alloc(nindex);
        gsl_matrix *evec = gsl_matrix_alloc(nindex,nindex);
        gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(nindex);
        gsl_eigen_symmv(&m.matrix,eval,evec,w);
        gsl_eigen_symmv_free(w);
  
//      for (int i=0; i<nindex; i++) {vecc[i] = gsl_matrix_get(evec,i,0);}
      for (int i=0; i<nindex; i++) {
        coulc.evals3[i] = gsl_vector_get(eval,i);
        for (int j=0; j<nindex; j++) {
          coulc.evecs3[i+j*nindex] = gsl_matrix_get(evec,i,j);
        }
      }

      /** Get perturbative correction **/
      pertCalcEigen(molc,coulc,energies,int3,intham);
      //pertCalcElements(molc,coulc,int3,energies);
      /** Write the coupling to file **/
//CTCs Change printing conditions for different configurations
//C1 - Prints DA1, DA1A2, J
      if (read.calc.C1_) {
        //appendData3d(cfile,
        //        mol[1].com[read.mol[1].mv[0].iaxis],
        //        mol[2].com[read.mol[2].mv[0].iaxis]-mol[1].com[read.mol[1].mv[0].iaxis],
        //        intham[read.calc.istate + read.calc.fstate*mol[0].nindices]);
        pos1[r1] = posMol1;//molc[1].com[read.mol[1].mv[0].iaxis];
        pos2[r2] = /*molc[2].com[read.mol[2].mv[0].iaxis]*/posMol2-molc[1].com[read.mol[1].mv[0].iaxis];
        sssfCoupling[r1+r2*read.mol[1].mv[0].steps] = 
      pertCalcElements(molc,coulc,int3,energies);//intham[read.calc.istate + read.calc.fstate*mol[0].nindices];
      } else if (read.calc.C2_) {
//C2 - Prints DA1, DA2, J
        pos1[r1] = posMol1;
        pos2[r2] = posMol2;
        sssfCoupling[r1+r2*read.mol[1].mv[0].steps] =
      pertCalcElements(molc,coulc,int3,energies);// intham[read.calc.istate + read.calc.fstate*mol[0].nindices];
        /*print.appendData3d(cfile,
                mol[1].com[read.mol[1].mv[0].iaxis],
                mol[2].com[read.mol[2].mv[0].iaxis],
                intham[read.calc.istate + read.calc.fstate*mol[0].nindices]);
      */
      } else if (read.calc.C3_) {
//C3 - Prints DA1_i, A1A2_j
        print.appendData3d(cfile,
                mol[1].com[read.mol[1].mv[0].iaxis],
                -mol[2].com[read.mol[2].mv[1].iaxis]+mol[1].com[read.mol[1].mv[1].iaxis],
                intham[read.calc.istate + read.calc.fstate*mol[0].nindices]);
      }
      
      /** Write the closest approach cross section in each direction **/
      if (read.calc.C2_ && (mol[1].com[read.mol[1].mv[0].iaxis]==mol[1].icom[read.mol[1].mv[0].iaxis])) {
        /** For use with C2, this prints DA2 **/
        print.appendData2d(crossfile1,mol[2].com[read.mol[2].mv[0].iaxis],
                intham[read.calc.istate + read.calc.fstate*mol[0].nindices]);
      }
      if (closest) {
        /** For use with C1 and C2, this prints DA1 **/
        print.appendData2d(crossfile2,mol[1].com[read.mol[1].mv[0].iaxis],
                intham[read.calc.istate + read.calc.fstate*mol[0].nindices]);
        closest = false;
      }
      if (read.calc.C1_ ) {
        if (mol[1].com[read.mol[1].mv[0].iaxis]==mol[1].icom[read.mol[1].mv[0].iaxis]) {
          /** For use with C1, this prints A1A2 for closest DA1 **/
          print.appendData2d(crossfile3,mol[2].com[read.mol[2].mv[0].iaxis]-mol[1].com[read.mol[1].mv[0].iaxis],
                intham[read.calc.istate + read.calc.fstate*mol[0].nindices]);
        }
        if (r1 == r2) {
         /** For use with C1, this prints A1=A2 **/
          print.appendData2d(crossfile4,molc[1].com[read.mol[1].mv[0].iaxis],
                sssfCoupling[r1+r2*read.mol[1].mv[0].steps]);
        }
      }
      if (read.calc.C2_ && (r1==r2)) {
        /** For use with C2, this prints A1=A2 **/
        print.appendData2d(crossfile4,molc[1].com[read.mol[1].mv[0].iaxis],
                sssfCoupling[r1+r2*read.mol[1].mv[0].steps]);
      }
      /** Move molecule 2 **/
      //molc[2].translate(read.mol[2].mv[0].iaxis,molc[2].grid[read.mol[2].mv[0].iaxis].dgrid);
      /** If configuration 3, then move molecule one the same amount **/
      if (read.calc.C3_) {
        mol[1].translate(read.mol[1].mv[0].iaxis,mol[1].grid[read.mol[2].mv[0].iaxis].dgrid);
      }
     
      /** If molecule 2 is out of bounds, skip ahead **/
      if (molc[2].com[read.mol[2].mv[0].iaxis]*molc[2].com[read.mol[2].mv[0].iaxis] > read.mol[2].mv[0].max*read.mol[2].mv[0].max) {
        continue;
      }
      
      molc[2].setCom();
      molc[1].setCom();
      //print.appendData2d(cmfile,mol[2].com[0],mol[1].com[0]);
    }//end move 2

    if (read.calc.C3_) {
      mol[1].resetExcept(read.mol[1].mv[1].iaxis);
      mol[2].resetExcept(read.mol[2].mv[1].iaxis);
      
      mol[1].moveTo(read.mol[1].mv[1].iaxis,
              mol[1].com[read.mol[1].mv[1].iaxis]
              +mol[1].grid[read.mol[1].mv[1].iaxis].dgrid);
      mol[2].moveTo(read.mol[2].mv[1].iaxis,
              mol[2].com[read.mol[2].mv[1].iaxis]
              +mol[2].grid[read.mol[2].mv[1].iaxis].dgrid);
    
    }
    closest=true;
    } //end move 1
} //end omp
    //print stuff to file  
    cout<<"Printing"<<endl;
      for (int r1=0; r1<read.mol[1].mv[0].steps; r1++) {
        for (int r2=0; r2<read.mol[2].mv[0].steps; r2++) {
          print.appendData3d(cfile,pos1[r1],
                              pos2[r2],
                              sssfCoupling[r1+r2*read.mol[1].mv[0].steps]);
        }
      }
    break; //end case 3
     

  
  } //end switch

  return 0;
}


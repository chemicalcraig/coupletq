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
  cout<<"*** Start Input file ***"<<endl;
  while (!efile.eof()) {
    efile.getline(tmpc,1000);
    cout<<tmpc<<endl;
  }
  efile.close();
  cout<<"*** End Input File ***"<<endl;
  
  /** Read input com file **/
  Reader read(argv[1]);
  if (read.calc.configuration.compare(0,2,"c1",0,2)==0) {
    C1_=true;
  } else if (read.calc.configuration.compare(0,2,"c2",0,2)==0) {
    C2_=true;
  } else if (read.calc.configuration.compare(0,2,"c3",0,2)==0) {
    C3_=true;  
  }
  mol = initialize(read);
  
  /** Set COM before initial translation**/
  for (int i=0; i<read.calc.molecules; i++) {
    mol[i].setCom();
  }

  /** set up printer and output files **/
  Print print(mol);
  ofstream outfile2,pdafile;
  remove(mol[0].outputfilename.c_str());
  outfile2.open(mol[0].outputfilename.c_str(), std::ofstream::out | std::ofstream::app);
  outfile2.precision(16);

  /** Set up the grid data **/
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


  /** Print initial positions **/
  ofstream posout;
  posout.open("initpos.xyz");
  int totalAtoms = 0;
  for (int i=0; i<read.calc.molecules; i++) {
    totalAtoms += mol[i].natoms;
  }
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
  
  /** Set indicies matrix **/
  int nindex=1;

  for (int i=0; i<read.calc.molecules; i++) {
    nindex *= mol[i].nstates;
  }
  
  int3 = new double[nindex*nindex];
  intham = new double[nindex*nindex];

  mol[0].nindices = nindex;
  setIndices(mol,mol[0].nmol,nindex);
  for (int i=0; i<nindex; i++) 
    for (int j=0; j<mol[0].nmol; j++) 
      cout<<"state "<<i<<" mol "<<j<<" "<<mol[0].indices[j+i*mol[0].nmol]<<endl;

  /** Create Coulomb Matrix **/
  double temp[nindex*nindex],temp2[nindex*nindex];
  Coulomb coul(nindex);
  
  /** Create state energy matrix **/
  double energies[nindex];
  for (int i=0; i<nindex; i++) {
    energies[i] = 0.;
    for (int m=0; m<mol[0].nmol; m++) {
        energies[i] += mol[m].excenergy[mol[0].indices[m+i*mol[0].nmol]];
    }
    energies[0] = 0.;

    cout<<"energies "<<i<<" "<<energies[i]<<endl;
  }

  /** Create unfiltered Coulomb matrix **/
  createCoulomb3(mol,coul);
  //createCoulomb3(mol,int3);
 
  /** Filter Coulomb Matrix for energy conservation **/
  for (int i=0; i<nindex; i++) {
    //coul.int3[i+i*nindex] += energies[i];
    for (int j=0; j<nindex; j++) {
      
      //coul.int3[i+j*nindex] *= window(energies[i],energies[j],read.calc.ewindow,0);
      int3[i+j*nindex] = coul.int3[i+j*nindex];
      int3[i+j*nindex] *= window(energies[i],energies[j],read.calc.ewindow,0);
      //cout<<"coul.int3 "<<i<<" "<<j<<" "<<coul.int3[i+j*nindex]<<" "<<int3[i+j*nindex]<<endl;
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
  
  double vec[nindex],vec2[nindex];
  for (int i=0; i<nindex; i++) {vec[i] = gsl_matrix_get(evec,i,0);}
  double sum = 0.;
  for (int i=0; i<nindex; i++) {
    sum = 0.;
    cout<<"evals "<<i<<" "<<gsl_vector_get(eval,i)<<endl;
    coul.evals3[i] = gsl_vector_get(eval,i);
    //cout<<"evals "<<i<<" "<<coul.evals3[i]<<endl;
    for (int j=0; j<nindex; j++) {
      sum += coul.int3[i+j*nindex]*vec[j];
      coul.evecs3[i+j*nindex] = gsl_matrix_get(evec,i,j);
      cout<<"evecs "<<i<<" "<<j<<" "<<coul.evecs3[i+j*nindex]<<endl;
    }
  }
/*
double tildeint[nindex*nindex],tempm[nindex*nindex];
double vec1[nindex];
vec1[0] = 1;
//Make C.V.C^T
  cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,mol[0].nindices,
          mol[0].nindices,mol[0].nindices,1.,coul.evecs3,
          mol[0].nindices,coul.int3,mol[0].nindices,0,tempm,mol[0].nindices);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,mol[0].nindices,
          mol[0].nindices,mol[0].nindices,1.,tempm,mol[0].nindices,
          coul.evecs3,mol[0].nindices,0,tildeint,mol[0].nindices);
 for (int i=0; i<mol[0].nindices; i++) {
    for (int j=0; j<mol[0].nindices; j++) {
      cout<<i<<" "<<j<<" "<<tempm[i+j*mol[0].nindices]<<" "
        <<tildeint[i+j*mol[0].nindices]<<endl;
    }
  }

  cblas_dgemv(CblasColMajor,CblasNoTrans,mol[0].nindices,
          mol[0].nindices,1.,coul.evecs3,mol[0].nindices,vec1,1,0.,vec2,1);

  cblas_dgemv(CblasColMajor,CblasTrans,mol[0].nindices,
          mol[0].nindices,1.,coul.evecs3,mol[0].nindices,vec2,1,0.,vec1,1);
  for (int i=0; i<mol[0].nindices; i++) {
   // cout<<i<<" transformed vec "<<vec1[i]<<endl;
  }
 


//CTC e

/*******************  Done Setting up molecules *****************************/

  /************************************
   * Do the appropriate calulation
   */

  double angle = 0.;
  double slip;
  double coupling2;

  switch (read.calc.itype) {
    case 1:
//      for (int zi=0; zi<mol[1].grid[2].ngrid; zi++) {

        for (int zi=0; 1; zi++) {
        //slip = mol[0].grid[1].min;
        //for (int thetai=0; thetai<mol[0].grid.ntheta; thetai++) {
        //for (int islip=0; islip<mol[0].grid[2].ngrid; islip++) {
          fretCalc(mol,coupling);        
          cout<<"Performing a FRET calculation now: "<<coupling*27.211396<<endl;
          exit(0);
          print.appendData2d(outfile2,mol[1].grid[2].min+zi*mol[1].grid[2].dgrid,coupling);
          //print.appendData3d(outfile2,mol[0].grid[2].min+zi*mol[0].grid[2].dgrid,slip,coupling);
          //print.appendData2d(pdafile,mol[0].grid[2].min+zi*mol[0].grid[2].dgrid,coupling2);
          
          //mol[0].rotateCom(mol[0].grid.dtheta,mol[1].com);
          //mol[0].translate(1,mol[0].grid[1].dgrid);
          //slip += mol[0].grid[1].dgrid;
        //} //end slip
        //reset x and y coordinates, keep z coordinate
        //mol[0].resetall();
        //mol[0].resetExcept(2);
        //translate vertically
        mol[1].translate(2,mol[1].grid[2].dgrid);
        //mol[0].setCom();
      }
      break;
    case 2:
      double dum;
      //gsl_matrix_view m = gsl_matrix_view_array(coul.int3,nindex,nindex);
      //gsl_vector *eval = gsl_vector_alloc(nindex);
      //gsl_matrix *evec = gsl_matrix_alloc(nindex,nindex);
      //gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(nindex);
      //gsl_eigen_symmv(&m.matrix,eval,evec,w);

      for (int zi=0; zi<1; zi++) {
        pertCalcEigen(mol,coul,energies,int3,intham);
        propagateTime(mol,coul,energies,read.dyn.tstart,
                      read.dyn.tfinish,read.dyn.increment,intham,read);

      }
      break;
    case 3:
    //NB This is set up for displacement in only one direction at the moment
      int whichaxis = 2;
      bool closest = true;
      ofstream cfile,cmfile,crossfile1,crossfile2;;
      cfile.open("coupling-3.dat");
      cmfile.open("com.dat");
      crossfile1.open("cross1.dat");
      crossfile2.open("cross2.dat");

      double r12,r23;
      double minsep = read.mol[2].mv[0].min - read.mol[1].mv[0].min;
      
      /** NB for C3 the number of steps must be the same for all
       * directions. The directions must also be listed in the same
       * order for molecules 1 and 2 in the input file (maybe, haven't
       * tested it but better safe than sorry **/
      for (int r1=0; r1<read.mol[1].mv[0].steps; r1++) {
        mol[1].setCom();
        mol[2].setCom();

        for (int r2=0; r2<read.mol[2].mv[0].steps; r2++) {
//CTCs change this while condition for C1 or C2
//C1
//#ifdef C1_
//      while (mol[2].com[read.mol[2].mv[0].iaxis] 
//            < read.mol[2].mv[0].max+mol[1].com[read.mol[1].mv[0].iaxis]) {
//C2
//#elif defined( C2_ )
//      while (mol[2].com[read.mol[2].mv[0].iaxis] 
//            > read.mol[2].mv[0].max) {
//C3
//#elif defined(C3_)
//      while (mol[2].com[read.mol[2].mv[1].iaxis] < read.mol[2].mv[1].max) {
//#endif
//CTCe
      createCoulomb3(mol,coul);
      /** Filter Coulomb Matrix for energy conservation **/
        for (int i=0; i<nindex; i++) {
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
      for (int i=0; i<nindex; i++) {
        coul.evals3[i] = gsl_vector_get(eval,i);
        for (int j=0; j<nindex; j++) {
          coul.evecs3[i+j*nindex] = gsl_matrix_get(evec,i,j);
        }
      }

      /** Get perturbative correction **/
      pertCalcEigen(mol,coul,energies,int3,intham);

      /** Write the coupling to file **/
//CTCs Change printing conditions for different configurations
//C1 - Prints DA1, DA1A2, J
      if (C1_) {

        print.appendData3d(cfile,
                mol[1].com[read.mol[1].mv[0].iaxis],
                mol[2].com[read.mol[2].mv[0].iaxis]-mol[1].com[read.mol[1].mv[0].iaxis],
                intham[read.calc.istate + read.calc.fstate*mol[0].nindices]);
      } else if (C2_) {
//C2 - Prints DA1, DA2, J
        print.appendData3d(cfile,
                mol[1].com[read.mol[1].mv[0].iaxis],
                mol[2].com[read.mol[2].mv[0].iaxis],
                intham[read.calc.istate + read.calc.fstate*mol[0].nindices]);
      } else if (C3_) {
//C3 - Prints DA1_i, A1A2_j
        print.appendData3d(cfile,
                mol[1].com[read.mol[1].mv[0].iaxis],
                -mol[2].com[read.mol[2].mv[1].iaxis]+mol[1].com[read.mol[1].mv[1].iaxis],
                intham[read.calc.istate + read.calc.fstate*mol[0].nindices]);
      }
      
      /** Write the closest approach cross section in each direction **/
      if (mol[1].com[read.mol[1].mv[0].iaxis]==mol[1].icom[read.mol[1].mv[0].iaxis]) {
        /** Cross1 **/
        print.appendData2d(crossfile1,mol[2].com[read.mol[2].mv[0].iaxis],
                intham[read.calc.istate + read.calc.fstate*mol[0].nindices]);
      }
      if (closest) {
        /** Cross2 **/
        print.appendData2d(crossfile2,mol[1].com[read.mol[1].mv[0].iaxis],
                intham[read.calc.istate + read.calc.fstate*mol[0].nindices]);
        closest = false;
      }

      /** Move molecule 2 **/
      mol[2].translate(read.mol[2].mv[0].iaxis,mol[2].grid[read.mol[2].mv[0].iaxis].dgrid);
      /** If configuration 3, then move molecule one the same amount **/
      if (C3_) {
        mol[1].translate(read.mol[1].mv[0].iaxis,mol[1].grid[read.mol[2].mv[0].iaxis].dgrid);
      }
      
      /** If molecule 2 is out of bounds, skip ahead **/
      if (mol[2].com[read.mol[2].mv[0].iaxis]*mol[2].com[read.mol[2].mv[0].iaxis] > read.mol[2].mv[0].max*read.mol[2].mv[0].max) {
        continue;
      }

      
      mol[2].setCom();
      mol[1].setCom();
      //print.appendData2d(cmfile,mol[2].com[0],mol[1].com[0]);
    }//end move 2
    
    if (!C3_) {
      mol[1].translate(read.mol[1].mv[0].iaxis,mol[1].grid[read.mol[1].mv[0].iaxis].dgrid);
      mol[1].setCom();
      mol[2].resetall();
      if (read.mol[2].mv[0].min >= 0) {
        mol[2].moveTo(read.mol[2].mv[0].iaxis,mol[1].com[read.mol[2].mv[0].iaxis]+minsep);
      }
    } else {
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
    break;

    }
  return 0;
}


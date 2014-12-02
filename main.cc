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
  
  /** Read input com file **/
  Reader read(argv[1]);
  mol = initialize(read);
  
  /** Set COM **/
  for (int i=0; i<read.calc.molecules; i++) {
    mol[i].setCom();
  }

  /** scale state charges for pmi, nelec=174 **/
  for (int i=0; i<mol[0].nmol; i++) {
    for (int j=0; j<mol[i].natoms; j++) {
      for (int k=0; k<mol[i].nstates; k++) {
        mol[i].atoms[j].charges[k+k*mol[i].nstates] *= 1;//74;
      }
    }
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
  /** Calculate transition dipole from charges **/
  /** This also sets the transition vector elements **/
  for (int i=0; i<read.calc.molecules; i++) {
    mol[i].setPostoInit();
  }

  //reset initial positions
  arrangeMol(mol);
 
  /** Print initial positions **/
  ofstream posout;
  posout.open("initpos.dat");
  for (int i=0; i<read.calc.molecules; i++) {
    for (int j=0; j<mol[i].natoms; j++) {
      for (int k=0; k<3; k++) {
        posout<<mol[i].atoms[j].pos[k]<<" ";
      }
      posout<<endl;
    }
  }

  /** Set indicies matrix **/
  int nindex=1;

  for (int i=0; i<read.calc.molecules; i++)
    nindex *= mol[i].nstates;
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
    //Convert to eV
    //energies[i] *= 27.211396;
    cout<<"energies "<<i<<" "<<energies[i]<<endl;
  }
  createCoulomb3(mol,coul);
  createCoulomb3(mol,int3);
 
  /** Filter Coulomb Matrix for energy conservation **/
  for (int i=0; i<nindex; i++) {
    for (int j=0; j<nindex; j++) {
      int3[i+j*nindex] = coul.int3[i+j*nindex];
      if (i!=j) {
        //convert from au to eV
        //int3[i+j*nindex] *= 1.;
        //coul.int3[i+j*nindex] *= 1.;
        //int3[i+j*nindex] *= window(energies[i],energies[j],read.calc.ewindow,0);
      }
      cout<<"coul.int3 "<<i<<" "<<j<<" "<<coul.int3[i+j*nindex]<<" "<<int3[i+j*nindex]<<endl;
    }
  }
  

  /** Get eigensystem of Coulomb matrix **/
  gsl_matrix_view m = gsl_matrix_view_array(int3,nindex,nindex);
  gsl_vector *eval = gsl_vector_alloc(nindex);
  gsl_matrix *evec = gsl_matrix_alloc(nindex,nindex);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(nindex);
  gsl_eigen_symmv(&m.matrix,eval,evec,w);
  gsl_eigen_symmv_free(w);
  
  createCoulomb3(mol,coul,energies,read);
  createCoulomb3(mol,int3);
 
  //coul.diagonalize(8,coul.evecs3,coul.evals3,coul.int3);
  double vec[nindex],vec2[nindex];
  for (int i=0; i<nindex; i++) {vec[i] = gsl_matrix_get(evec,i,0);}
  double sum = 0.;
  for (int i=0; i<nindex; i++) {
    sum = 0.;
    coul.evals3[i] = gsl_vector_get(eval,i);
    //cout<<"evals "<<i<<" "<<coul.evals3[i]<<endl;
    for (int j=0; j<nindex; j++) {
      sum += coul.int3[i+j*nindex]*vec[j];
      coul.evecs3[i+j*nindex] = gsl_matrix_get(evec,i,j);
      cout<<"evecs "<<i<<" "<<j<<" "<<coul.evecs3[i+j*nindex]<<endl;
    }
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
/*      for (int i=0; i<mol[0].ngriddim; i++) { //x,y,z,...
        fretCalc(mol,coupling);
        for (int j=0; j<mol[0].grid[i].ngrid; j++) { //ngrid steps
          
        }

      }
*/
//      for (int zi=0; zi<mol[1].grid[2].ngrid; zi++) {

        for (int zi=0; 1; zi++) {
//angle = 0.;
        //slip = mol[0].grid[1].min;
        //for (int thetai=0; thetai<mol[0].grid.ntheta; thetai++) {
        //for (int islip=0; islip<mol[0].grid[2].ngrid; islip++) {

          fretCalc(mol,coupling);        
          
          cout<<"Performing a FRET calculation now: "<<coupling*27.211396<<endl;
          exit(0);
          print.appendData2d(outfile2,mol[1].grid[2].min+zi*mol[1].grid[2].dgrid,coupling);
          //print.appendData3d(outfile2,mol[0].grid[2].min+zi*mol[0].grid[2].dgrid,slip,coupling);

          //pdaCalc(mol,coupling2);
          //print.appendData3d(pdafile,mol[0].grid[2].min+zi*mol[0].grid[2].dgrid,slip,coupling2);
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
      int whichaxis = 0;
      //pertCalcDegen(mol,coul,energies,int3,dum);
      //pertCalc(mol,coul,intham,energies);
      gsl_matrix_view m = gsl_matrix_view_array(coul.int3,nindex,nindex);
      gsl_vector *eval = gsl_vector_alloc(nindex);
      gsl_matrix *evec = gsl_matrix_alloc(nindex,nindex);
      gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(nindex);
      gsl_eigen_symmv(&m.matrix,eval,evec,w);

      for (int zi=0; zi<1; zi++) {
      //for (int zi=0; zi<mol[1].grid[whichaxis].ngrid; zi++) {
        //pertCalcNonDegen(mol,coul,energies,int3,&dum);
        pertCalcDegen(mol,coul,energies,int3,dum,intham,read);
        propagateTime(mol,coul,energies,0,800,0.000001,intham,read);

        exit(0);
        print.appendData2d(outfile2,mol[1].grid[whichaxis].min+zi*mol[1].grid[whichaxis].dgrid,dum);
        mol[1].translate(whichaxis,mol[1].grid[whichaxis].dgrid);
        
        createCoulomb3(mol,coul);
        createCoulomb3(mol,int3);
        for (int i=0; i<nindex; i++) {
          coul.int3[i+i*nindex] += energies[i];
          for (int j=0; j<nindex; j++) {
            if (i!=j)
 //             coul.int3[i+j*nindex] *= 10;
          }
        }
        m = gsl_matrix_view_array(coul.int3,nindex,nindex);
        //w = gsl_eigen_symmv_alloc(nindex);
        gsl_eigen_symmv(&m.matrix,eval,evec,w);

        double vec[nindex],vec2[nindex];
        for (int i=0; i<nindex; i++) {vec[i] = gsl_matrix_get(evec,i,0);}
        double sum = 0.;
        for (int i=0; i<nindex; i++) {
          sum = 0.;
          coul.evals3[i] = gsl_vector_get(eval,i);
          //coul.int3[i+i*nindex] += energies[i];
          //cout<<"evals "<<i<<" "<<coul.evals3[i]<<endl;
          for (int j=0; j<nindex; j++) {
           // if (i!=j)
             // coul.int3[i+j*nindex] *= 10;
            coul.evecs3[i+j*nindex] = gsl_matrix_get(evec,i,j);
            cout<<"evecs "<<i<<" "<<j<<" "<<coul.evecs3[i+j*nindex]<<endl;
          }
        }
      }
      gsl_eigen_symmv_free(w);
      break;
    }
  return 0;
}


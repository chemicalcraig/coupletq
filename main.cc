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
  
  //open input file
  comfile.open(argv[1]);

  /*
   * Parse input comfile and retrieve necessary data
   * this sets up the main Molecule object
   */
  mol = parseComfile(comfile);

  /** Set COM **/
  for (int i=0; i<mol[0].nmol; i++) {
    mol[i].setCom();
  }
  
  /** set up printer and output files **/
  Print print(mol);
  ofstream outfile2,pdafile;
  remove(mol[0].outputfilename.c_str());
  outfile2.open(mol[0].outputfilename.c_str(), std::ofstream::out | std::ofstream::app);
  //remove("pda.dat");
  //pdafile.open("pda.dat", std::ofstream::out | std::ofstream::app);
  outfile2.precision(16);
  pdafile.precision(16);

  /** Set up the grid data **/
  int griddim = 3;
  for (int i=0; i<mol[0].nmol; i++) {
    mol[i].grid = new Grid[griddim];
    mol[i].griddim = griddim;
    for (int j=0; j<griddim; j++) {
      mol[i].grid[j].size = griddim;
    }
  }

  /** min, max, nsteps **/
  //mol[0].rotateTheta(M_PI/2,1);
  //mol[0].grid[0].setParams(-3.,-3.,1);

  mol[1].grid[0].setParams(3., 3., 1);
  //mol[1].grid[1].setParams(1.,1.,1);
  //mol[1].grid[2].setParams(3.,3.,1);
  //mol[1].rotateTheta(-1*M_PI/2,1);
  if (mol[0].interaction > 1) {
    mol[2].grid[0].setParams(-3., -3., 1);
    //mol[2].grid[1].setParams(-1.,-1.,1);
    //mol[2].grid[2].setParams(5.,5.,1);
    //mol[2].rotateTheta(1*M_PI/2,1);
    //mol[2].rotateTheta(M_PI,0);
  }


/*****************  Setting up Molecular distribution ******************/
  /** Calculate transition dipole from charges **/
  /** This also sets the transition vector elements **/
  for (int i=0; i<mol[0].nmol; i++) {
 //   calcdip(mol[i]);
    mol[i].setPostoInit();
  }


  //reset initial positions
  arrangeMol(mol);
  for (int i=0; i<mol[0].nmol; i++) {
//    mol[i].setCom();
//    mol[i].setPostoInit();
  }
 
  /** Print initial positions **/
  ofstream posout;
  posout.open("initpos.dat");
  for (int i=0; i<mol[0].nmol; i++) {
    for (int j=0; j<mol[i].natoms; j++) {
      for (int k=0; k<3; k++) {
        posout<<mol[i].atoms[j].pos[k]<<" ";
      }
      posout<<endl;
    }
  }
    
  //CTC test start
  double temp[64],temp2[64],temp3[4],temp4[4],temp5[4];
  Coulomb coul;
  createCoulomb3(mol,coul);
  createCoulomb3(mol,int3);
  double energies[8];
  for (int i=0; i<2; i++) 
    for (int j=0; j<2; j++) 
      for (int k=0; k<2; k++) {
        int index = i+j*2+k*4;
        energies[index] = (mol[0].excenergy[i] )//- mol[0].groundenergy)
                            + (mol[1].excenergy[j])// - mol[1].groundenergy)
                            + (mol[2].excenergy[k]);// - mol[2].groundenergy);
        energies[index] *= 27.211396;
cout<<"energies "<<index<<" "<<energies[index]<<endl;
        //if (mol[0].interaction>1)
          //coul.int3[index+index*8] += energies[index];
    //    cout<<index<<" energy index "<<coul.int3[index+index*8]<<endl;
      }
  /** Filter Coulomb Matrix for energy conservation **/
  

  
  for (int i=0; i<8; i++) {
    for (int j=0; j<8; j++) {
      if (i!=j) {
      coul.int3[i+j*8] *= 10.;
      //int3[i+j*8] = coul.int3[i+j*8];
      coul.int3[i+j*8] *= window(energies[i],energies[j],0.4,0);
      //cout<<i<<" "<<j<<" "<<window(energies[i],energies[j],2,0)<<"   window"<<endl;
      }
    }
  }

  gsl_matrix_view m = gsl_matrix_view_array(int3,8,8);
  gsl_vector *eval = gsl_vector_alloc(8);
  gsl_matrix *evec = gsl_matrix_alloc(8,8);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(8);
  gsl_eigen_symmv(&m.matrix,eval,evec,w);
  gsl_eigen_symmv_free(w);

  //coul.diagonalize(8,coul.evecs3,coul.evals3,coul.int3);
  double vec[8],vec2[8];
  for (int i=0; i<8; i++) {vec[i] = gsl_matrix_get(evec,i,0);}
  double sum = 0.;
  for (int i=0; i<8; i++) {
    sum = 0.;
    coul.evals3[i] = gsl_vector_get(eval,i);
    //cout<<"evals "<<i<<" "<<coul.evals3[i]<<endl;
    for (int j=0; j<8; j++) {
      sum += coul.int3[i+j*8]*vec[j];
      coul.evecs3[i+j*8] = gsl_matrix_get(evec,i,j);
      cout<<"evecs "<<i<<" "<<j<<" "<<coul.evecs3[i+j*8]<<endl;
    }
  }
//  exit(0);
/*  
  coul.diagonalize(8,coul.evecs3,coul.evals3,coul.int3);
  for (int i=0; i<8; i++) {
    for (int j=0; j<8; j++) {
      coul.evecs3[i+j*8] = temp[i+j*8];
    }
  }
*/
  /** Check if evects are orthogonal **/
  cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,8,8,8,1.,
              coul.evecs3,8,coul.int3,8,0,temp2,8);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,8,8,8,1.,
              coul.evecs3,8,coul.evecs3,8,0,temp,8);


  for (int i=0; i<8; i++) {
    for (int j=0; j<8; j++) {
      //cout<<i<<" "<<j<<" diagonal? "<<temp[i+j*8]<<endl;
      //coul.dint3[i+j*8] = temp[i+j*8];
    }
  }
  //exit(0);
  //coul.createCoulomb2(mol);
  //coul.diagonalize(coul.n2d,coul.evecs2,coul.evals2,coul.int2);
  for (int i=0; i<8; i++) {
    for (int j=0; j<8; j++) {
      //cout<<i<<" "<<j<<" "<<coul.evecs3[i+j*8]<<endl;
    }
  }
//  exit(0);
//CTC e

/*******************  Done Setting up molecules *****************************/

  /************************************
   * Do the appropriate calulation
   */

  double angle = 0.;
  double slip;
  double coupling2;

  switch (mol[0].interaction) {
    case 1:
/*      for (int i=0; i<mol[0].ngriddim; i++) { //x,y,z,...
        fretCalc(mol,coupling);
        for (int j=0; j<mol[0].grid[i].ngrid; j++) { //ngrid steps
          
        }

      }
*/
      for (int zi=0; zi<mol[1].grid[2].ngrid; zi++) {
        //angle = 0.;
        //slip = mol[0].grid[1].min;
        //for (int thetai=0; thetai<mol[0].grid.ntheta; thetai++) {
        //for (int islip=0; islip<mol[0].grid[2].ngrid; islip++) {
          fretCalc(mol,coupling);          
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
      gsl_matrix_view m = gsl_matrix_view_array(int3,8,8);
      gsl_vector *eval = gsl_vector_alloc(8);
      gsl_matrix *evec = gsl_matrix_alloc(8,8);
      gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(8);
      gsl_eigen_symmv(&m.matrix,eval,evec,w);


      for (int zi=0; zi<mol[1].grid[whichaxis].ngrid; zi++) {
        //pertCalcNonDegen(mol,coul,energies,int3,&dum);
        pertCalcDegen(mol,coul,energies,int3,dum,intham);
        propagateTime(mol,coul,energies,0,800,0.000001,intham);
        exit(0);
        print.appendData2d(outfile2,mol[1].grid[whichaxis].min+zi*mol[1].grid[whichaxis].dgrid,dum);
        mol[1].translate(whichaxis,mol[1].grid[whichaxis].dgrid);
        
        createCoulomb3(mol,coul);
        createCoulomb3(mol,int3);
        for (int i=0; i<8; i++) {
          coul.int3[i+i*8] += energies[i];
          for (int j=0; j<8; j++) {
            if (i!=j)
              coul.int3[i+j*8] *= 10;
          }
        }
        m = gsl_matrix_view_array(int3,8,8);
        //w = gsl_eigen_symmv_alloc(8);
        gsl_eigen_symmv(&m.matrix,eval,evec,w);

        double vec[8],vec2[8];
        for (int i=0; i<8; i++) {vec[i] = gsl_matrix_get(evec,i,0);}
        double sum = 0.;
        for (int i=0; i<8; i++) {
          sum = 0.;
          coul.evals3[i] = gsl_vector_get(eval,i);
          //coul.int3[i+i*8] += energies[i];
          //cout<<"evals "<<i<<" "<<coul.evals3[i]<<endl;
          for (int j=0; j<8; j++) {
           // if (i!=j)
             // coul.int3[i+j*8] *= 10;
            coul.evecs3[i+j*8] = gsl_matrix_get(evec,i,j);
            cout<<"evecs "<<i<<" "<<j<<" "<<coul.evecs3[i+j*8]<<endl;
          }
        }
      }
      gsl_eigen_symmv_free(w);
      break;
    }
  return 0;
}


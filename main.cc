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
  mol[0].grid[2].setParams(4., 4., 1);
  //mol[0].grid[2].setParams(10., 12., 200);
  mol[2].grid[1].setParams(-1000., -12., 200);

/*****************  Setting up Molecular distribution ******************/
  /** Calculate transition dipole from charges **/
  /** This also sets the transition vector elements **/
  for (int i=0; i<mol[0].nmol; i++) {
    calcdip(mol[i]);
    mol[i].setPostoInit();
  }

  //Rotate Molecules if necessary to align mainly in the XY plane
//  for (int i=0; i<mol[0].nmol; i++)
//    mol[i].rotateTheta(M_PI/2,1);

  //reset initial positions
  arrangeMol(mol);
  for (int i=0; i<mol[0].nmol; i++) {
    mol[i].setCom();
    mol[i].setPostoInit();
  }
  
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
      for (int zi=0; zi<mol[0].grid[2].ngrid; zi++) {
        //angle = 0.;
        slip = mol[0].grid[1].min;
        //for (int thetai=0; thetai<mol[0].grid.ntheta; thetai++) {
        //for (int islip=0; islip<mol[0].grid[2].ngrid; islip++) {
          fretCalc(mol,coupling);          
          print.appendData2d(outfile2,mol[0].grid[2].min+zi*mol[0].grid[2].dgrid,coupling);
          //print.appendData3d(outfile2,mol[0].grid[2].min+zi*mol[0].grid[2].dgrid,slip,coupling);

          //pdaCalc(mol,coupling2);
          //print.appendData3d(pdafile,mol[0].grid[2].min+zi*mol[0].grid[2].dgrid,slip,coupling2);
          print.appendData2d(pdafile,mol[0].grid[2].min+zi*mol[0].grid[2].dgrid,coupling2);
          
          //mol[0].rotateCom(mol[0].grid.dtheta,mol[1].com);
          //mol[0].translate(1,mol[0].grid[1].dgrid);
          //slip += mol[0].grid[1].dgrid;
        //} //end slip

        //reset x and y coordinates, keep z coordinate
        //mol[0].resetall();
        mol[0].resetExcept(2);

        //translate vertically
        mol[0].translate(2,mol[0].grid[2].dgrid);
        //mol[0].setCom();
      }
      break;
    case 2:
      pertCalc(mol);
      break;
    }
  return 0;
}


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
  mol[0].setCom();
  mol[1].setCom();
  
  /** set up printer and output files **/
  Print print(mol);
  ofstream outfile2,pdafile;
  remove(mol[0].outputfilename.c_str());
  outfile2.open(mol[0].outputfilename.c_str(), std::ofstream::out | std::ofstream::app);
  remove("pda.dat");
  pdafile.open("pda.dat", std::ofstream::out | std::ofstream::app);
  outfile2.precision(16);
  pdafile.precision(16);

  /** Set up the grid data for translation/rotations
   * in the Etransfer calculation **/
  mol[0].grid.setParams(100,100,200,
    0.05,0.4,0.05,
    2*M_PI,100,mol[0].grid.thetamax/mol[0].grid.ntheta);

/*****************  Setting up Molecular distribution ******************/
  /** Calculate transition dipole from charges **/
  /** This also sets the transition vector elements **/
  for (int i=0; i<mol[0].nmol; i++) {
    calcdip(mol[i]);
    mol[i].setPostoInit();
  }

  //Rotate Molecules if necessary to align mainly in the XY plane
  for (int i=0; i<mol[0].nmol; i++)
    mol[i].rotateTheta(M_PI/2,1);

  //move the first molecule a bit along the z-axis
  const double transz = 3.;
  mol[0].translate(2,transz);
//  mol[0].com[2] += transz;
  
  //move molecule along y axis to slip
  const double transy = -20;
  mol[0].translate(1,transy);
//  mol[0].com[1] += transy;

  //reset initial positions
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
  ofstream outfile3;
  double coupling2;
  outfile3.open("lastcoord");
  switch (mol[0].interaction) {
    case 1:
      for (int zi=0; zi<mol[0].grid.nz; zi++) {
        //angle = 0.;
        slip = transy;
        //for (int thetai=0; thetai<mol[0].grid.ntheta; thetai++) {
        for (int islip=0; islip<mol[0].grid.ny; islip++) {
          fretCalc(mol,coupling);          
          //print.appendData2d(outfile2,transz+zi*mol[0].grid.dz,coupling);
          print.appendData3d(outfile2,transz+zi*mol[0].grid.dz,slip,coupling);

          pdaCalc(mol,coupling2);
          print.appendData3d(pdafile,transz+zi*mol[0].grid.dz,slip,coupling2);
          //print.appendData2d(pdafile,transz+zi*mol[0].grid.dz,coupling2);
          
          //outfile2<<trans+zi*mol[0].grid.dz<<" "<<angle<<" "<<coupling<<endl;
          //mol[0].rotateCom(mol[0].grid.dtheta,mol[1].com);
          mol[0].translate(1,mol[0].grid.dy);
          slip += mol[0].grid.dy;
         // mol[0].setCom();
          //angle += mol[0].grid.dtheta;
        } //end slip

        //reset x and y coordinates, keep z coordinate
        //mol[0].resetall();
        mol[0].resetExcept(2);

        //translate vertically
        mol[0].translate(2,mol[0].grid.dz);
        //mol[0].setCom();
      }
      break;
    case 2:
      pertCalc(mol);
      break;
    }
  return 0;
}


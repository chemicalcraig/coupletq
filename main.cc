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

  //parse comfile and retrieve necessary data
  mol = parseComfile(comfile);

  //set up printer
  Print print(mol);

  //perform either a FRET or 3-body calculation
  ofstream outfile2;
  outfile2.open(mol[0].outputfilename.c_str(), std::ofstream::out | std::ofstream::app);
  outfile2.precision(16);
  
  /** Set up the grid data for translation/rotations
   * in the Etransfer calculation
   */
  mol[0].grid.setParams(100,100,200,
    0.05,0.4,0.05,
    2*M_PI,100,mol[0].grid.thetamax/mol[0].grid.ntheta);

  //Rotate Molecules if necessary to align mainly in the XY plane
  for (int i=0; i<mol[0].nmol; i++)
    mol[i].rotateTheta(M_PI/2,1);

  //move the first molecule a bit along the z-axis
  const double transz = 3.;
  mol[0].translate(2,transz);
  mol[0].com[2] += transz;
  
  //move molecule along y axis to slip
  const double transy = -20;
  mol[0].translate(1,transy);
  mol[0].com[1] += transy;

  //reset initial positions
  for (int i=0; i<mol[0].nmol; i++)
    mol[i].setPostoInit();
  
  /************************************
   * Do the appropriate calulation
   */

  double angle = 0.;
  double slip;
  ofstream outfile3;
  outfile3.open("lastcoord");
  switch (mol[0].interaction) {
    case 1:
      for (int zi=0; zi<mol[0].grid.nz; zi++) {
        //angle = 0.;
        slip = transy;
        //for (int thetai=0; thetai<mol[0].grid.ntheta; thetai++) {
        for (int islip=0; islip<mol[0].grid.ny; islip++) {
          fretCalc(mol,coupling);
          print.appendData3d(outfile2,transz+zi*mol[0].grid.dz,slip,coupling);
          //outfile2<<trans+zi*mol[0].grid.dz<<" "<<angle<<" "<<coupling<<endl;
          //mol[0].rotateCom(mol[0].grid.dtheta,mol[1].com);
          mol[0].translate(1,mol[0].grid.dy);
          slip += mol[0].grid.dy;
          //angle += mol[0].grid.dtheta;
        }

        //reset x and y coordinates, keep z coordinate
        //mol[0].resetall();
        mol[0].resetExcept(2);

        //translate vertically
        //mol[0].translate(2,trans+zi*mol[0].grid.dz);
        mol[0].translate(2,mol[0].grid.dz);
      }
      break;
    case 2:
      pertCalc(mol);
      break;
    }
  return 0;
}


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
  mol = parseComfile(comfile,interactionOrder);

  //perform either a FRET or 3-body calculation
  switch (interactionOrder) {
    case 1:
      fretCalc(mol);
    case 2:
      pertCalc(mol);
  }

  return 0;
}


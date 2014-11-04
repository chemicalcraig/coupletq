#include "grid.h"
 
Grid::Grid() {
  this->min = 0.;
  this->max = 0.;
  this->dgrid = 0.;
  this->ngrid = 0;
}

void Grid::setParams(double min, double max, double dgrid) {

  this->min = min;
  this->max = max;
  this->dgrid = dgrid;

  ngrid = (this->max - this->min)/this->dgrid;
}

void Grid::setParams(double min, double max, int ngrid) {
  
  this->min = min;
  this->max = max;
  this->ngrid = ngrid;

  dgrid = (this->max - this->min)/this->ngrid;
}

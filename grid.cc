#include "grid.h"
  
void Grid::setParams(int nx,int ny,int nz,double dx,double dy,double dz,double tmax,int ntheta,double dtheta) {

  this->nx = nx;
  this->ny = ny;
  this->nz = nz;
  this->dx = dx;
  this->dy = dy;
  this->dz = dz;
  this->thetamax = tmax;
  this->ntheta = ntheta;
  this->dtheta = dtheta;
}

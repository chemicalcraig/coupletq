#ifndef GRID_H
#define GRID_H

using namespace std;

class Grid {
  public:
    double dx,dy,dz,dtheta;
    double xmin, xmax, ymin,ymax, zmin, zmax;
    double thetamin, thetamax;
    int nx,ny,nz,ntheta;

  Grid() {};
  ~Grid() {};
};
#endif

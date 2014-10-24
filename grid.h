#ifndef GRID_H
#define GRID_H

using namespace std;

class Grid {
  public:
    double dx,dy,dz,dtheta;
    double xmin, xmax, ymin,ymax, zmin, zmax;
    double thetamin, thetamax;
    int nx,ny,nz,ntheta;

  void setParams(int nx,int ny,int nz,double dx,double dy,double dz,double tmax,int ntheta,double dtheta);
  Grid() {};
  ~Grid() {};
};
#endif

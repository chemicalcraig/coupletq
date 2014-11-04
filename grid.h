#ifndef GRID_H
#define GRID_H

using namespace std;

class Grid {
  public:
    int size;
    double min,max,dgrid;
    int ngrid;

  void setParams(double min, double max, double dgrid);
  void setParams(double min, double max, int ngrid);
  Grid();
  ~Grid() {};
};
#endif

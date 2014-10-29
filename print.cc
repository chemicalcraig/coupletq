#include "print.h"

Print::Print(Molecule *mol) {
  this->mol = mol;
}

void Print::positions(string filename,const int which) {
  ofstream outfile;
  outfile.open(filename.c_str());

  for (int i=0; i<this->mol[which].natoms; i++) {
    outfile<<i<<" "<<this->mol[which].atoms[i].type<<" "
    <<this->mol[which].atoms[i].pos[0]<<" "
    <<this->mol[which].atoms[i].pos[1]<<" "
    <<this->mol[which].atoms[i].pos[2]<<endl;
  }
  outfile.close();
}

void Print::appendData3d(ofstream &outfile, double x, double y, double z) {
  outfile<<x<<" "<<y<<" "<<z<<endl;
}

void Print::appendData2d(ofstream &outfile, double x, double y) {
  outfile<<x<<" "<<y<<endl;
}



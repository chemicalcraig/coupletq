#include "reader.h"
#include "math.h"
#include <cstdlib>
#include <list>
#include <algorithm>

using namespace std;

/************************
* List of directives
* and keywords
* *********************/
const string dirs_[] = {"calculation","Calculation","CALCULATION",
                      "molecules","Molecules","MOLECULES",
                      "dynamics","Dynamics","DYNAMICS",
                      "fret","Fret","FRET"};
const string subdirs_[] = {"charges","Charges","CHARGES",
                         "move","Move","MOVE",
                         "rotate","Rotate","ROTATE",
                         "init","Init","INIT",
                         "output","Output","OUTPUT","configuration"};
const string opts_[] = {"type","ewindow","molecules","states","charges",
                           "move","rotate","tddft","start","finish","steps",
                           "increment","init","output","populations","wincrement","file","configuration","target"};
list<string> directives_(dirs_, dirs_+sizeof(dirs_)/sizeof(string));
list<string> subdirectives_(subdirs_, subdirs_+sizeof(subdirs_)/sizeof(string));
list<string> options_(opts_, opts_+sizeof(opts_)/sizeof(string));

/** Read Directive Block **/
void Reader::readBlock(string s1, ifstream &in, int molcount) {
  char c[1000];
  string s,s2;

  /** Calculation block **/
  if ((s1.compare(0,4,"calculation",0,4) == 0) ||
      (s1.compare(0,4,"Calculation",0,4) == 0) ||
      (s1.compare(0,4,"CALCULATION",0,4) == 0)) {
    while (s2.compare(0,3,"end",0,3) != 0) {
      string which="calc";
      in>>ws;
      in.getline(c,1000);
      s2=c;
      s=strtok(c," ");
      list<string>::iterator it = find(options_.begin(),options_.end(),s);
      if (it!=options_.end()) {
        /** Check for subdirective **/
        list<string>::iterator it2 = find(subdirectives_.begin(),subdirectives_.end(),*it);
        if (it2!=subdirectives_.end()) {
          //readSubBlock(which,*it2,in);
        } else {
          s=strtok(NULL," ");
          /** calc type **/
          if (string(*it).compare(0,4,"type",0,4)==0) {
            calc.type = s;
            if (s.compare(0,4,"pert",0,4)==0) {
              calc.itype = 2;
            } else if (s.compare(0,4,"fret",0,4)==0) {
              calc.itype = 1;
            } else if (s.compare(0,8,"coupling",0,8)==0) {
              calc.itype = 3;
              s=strtok(NULL," ");
              calc.istate = atoi(s.c_str());
              s=strtok(NULL," ");
              calc.fstate = atoi(s.c_str());
            }
          } else if (string(*it).compare(0,7,"ewindow",0,7)==0) {
            calc.ewindow = atof(s.c_str());
          } else if (string(*it).compare(0,9,"molecules",0,9)==0) {
            calc.molecules = atoi(s.c_str());
            mol = new Mol[calc.molecules];
          } else if (string(*it).compare(0,13,"configuration",0,13)==0) {
            calc.configuration = s;
            //if (s.compare(0,2,"c1",0,2)==0) {
            //  #define c1 1
            //} else if (s.compare(0,2,"c2",0,2) {
              //#define c2 1
            //}
            
          }
        }
      }
    } //found the end
  } //end calculation block
  
  /** Molecules block **/
  if ((s1.compare(0,4,"molecules",0,4) == 0) ||
      (s1.compare(0,4,"Molecules",0,4) == 0) ||
      (s1.compare(0,4,"MOLECULES",0,4) == 0)) {
    mol[molcount].nmov = 0;
    mol[molcount].nrot = 0;
    while (s2.compare(0,3,"end",0,3) != 0) {
      string which="mol";
      in>>ws;
      in.getline(c,1000);
      s2=c;
      s=strtok(c," ");
      list<string>::iterator it = find(options_.begin(),options_.end(),s);
      if (it!=options_.end()) {
        /** Check for subdirective **/
        list<string>::iterator it2 = find(subdirectives_.begin(),subdirectives_.end(),*it);
        if (it2!=subdirectives_.end()) {
          readSubBlock(which,*it2,in,molcount);
        } else {
          s=strtok(NULL," ");
          /** calc type **/
          if (string(*it).compare(0,6,"states",0,6)==0) {
            mol[molcount].nstates = atoi(s.c_str());
          } else if (string(*it).compare(0,5,"tddft",0,5)==0) {
            mol[molcount].tddftfile = s;
          } else if (string(*it).compare(0,6,"target",0,6)==0) {
            mol[molcount].target = atoi(s.c_str());
          }
        }
      }
    } //found the end
  } //end molecules block

/** Dynamics block **/
  if ((s1.compare(0,4,"dynamics",0,4) == 0) ||
      (s1.compare(0,4,"Dynamics",0,4) == 0) ||
      (s1.compare(0,4,"DYNAMICS",0,4) == 0)) {
    while (s2.compare(0,3,"end",0,3) != 0) {
      string which="dyn";
      in>>ws;
      in.getline(c,1000);
      s2=c;
      s=strtok(c," ");
      list<string>::iterator it = find(options_.begin(),options_.end(),s);
      if (it!=options_.end()) {
        /** Check for subdirective **/
        list<string>::iterator it2 = find(subdirectives_.begin(),subdirectives_.end(),*it);
        if (it2!=subdirectives_.end()) {
          readSubBlock(which,*it2,in,molcount);
        } else {
          s=strtok(NULL," ");
          /** calc type **/
          if (string(*it).compare(0,6,"start",0,6)==0) {
            dyn.tstart = atof(s.c_str());
          } else if (string(*it).compare(0,6,"finish",0,6)==0) {
            dyn.tfinish = atof(s.c_str());
          } else if (string(*it).compare(0,5,"steps",0,5)==0) {
            dyn.tsteps = atoi(s.c_str());
          } else if (string(*it).compare(0,9,"increment",0,9)==0) {
            dyn.increment = atof(s.c_str());
          } else if (string(*it).compare(0,10,"wincrement",0,10)==0) {
            double dum = atof(s.c_str());
            dyn.wstep = (int)(ceil(dum/dyn.increment));
          }
        }
      }
    } //found the end
  } //end dynamics block


};

/** Read Subdirective block **/
void Reader::readSubBlock(string which,string s1, ifstream &in,int molcount) {
  char c[1000];
  string s;
  /** Molecule sub blocks **/
  if (which.compare(0,3,"mol",0,3) == 0) {
    int pos = in.tellg();
    /** charges **/
    if (s1.compare(0,7,"charges",0,7) == 0) {
      int n=-1;
      while(s.compare(0,3,"end",0,3)!=0) {
        in>>ws;
        in.getline(c,1000);
        s=c;
        n++;
      }
      in.seekg(pos);
      mol[molcount].cf = new ChargeFile[n];
      mol[molcount].ncharges = n;
      for (int i=0; i<n; i++) {
        in>>ws;
        in.getline(c,1000);
        s=strtok(c," ");
        mol[molcount].cf[i].i = atoi(s.c_str());
        s=strtok(NULL," ");
        mol[molcount].cf[i].f = atoi(s.c_str());
        s=strtok(NULL," ");
        mol[molcount].cf[i].file = s;
      }
      in.getline(c,1000);
    } //end charges
  
    /** move **/
    if (s1.compare(0,4,"move",0,4) == 0) {
      int n=-1;
      int nc=0;
      int ncount=0;
      while(s.compare(0,3,"end",0,3)!=0) {
        in>>ws;
        in.getline(c,1000);
        s=c;
        if (s.compare(0,1,"#",0,1)==0) {
          //cout<<"line commented "<<s<<endl;
          nc++;
          //cout<<"nc = "<<endl;
        }
        n++;
      }
      in.seekg(pos);
      mol[molcount].mv = new Move[n-nc];
      mol[molcount].nmov = n-nc;
      
      int imv=0;
      nc=0;
      //while(s.compare(0,3,"end",0,3)!=0) {
      for (int i=0; i<n; i++) {
        in>>ws;
        in.getline(c,1000);
        s=strtok(c," ");
        if (s.compare(0,1,"#",0,1)==0) {
          nc++;
          continue;
        }
        mol[molcount].mv[i-nc].axis = s;
        if (s.compare(0,1,"x",0,1)==0)
          mol[molcount].mv[i-nc].iaxis = 0;
        else if (s.compare(0,1,"y",0,1)==0)
          mol[molcount].mv[i-nc].iaxis = 1;
        else if (s.compare(0,1,"z",0,1)==0)
          mol[molcount].mv[i-nc].iaxis = 2;
        s=strtok(NULL," ");
        mol[molcount].mv[i-nc].min = atof(s.c_str());
        s=strtok(NULL," ");
        mol[molcount].mv[i-nc].max = atof(s.c_str());
        s=strtok(NULL," ");
        mol[molcount].mv[i-nc].steps = atoi(s.c_str());
      }
      in.getline(c,1000);
    } //end move

    /** rotate **/
    if (s1.compare(0,6,"rotate",0,6) == 0) {
      int n=-1;
      while(s.compare(0,3,"end",0,3)!=0) {
        in>>ws;
        in.getline(c,1000);
        s=c;
        n++;
      }
      in.seekg(pos);
      mol[molcount].rot = new Rot[n];
      mol[molcount].nrot = n;
      
      for (int i=0; i<n; i++) {
        in>>ws;
        in.getline(c,1000);
        s=strtok(c," ");
        mol[molcount].rot[i].axis = s;
        s=strtok(NULL," ");
        mol[molcount].rot[i].theta = atof(s.c_str());
      }
      in.getline(c,1000);
    } //end rot
  } //end mol block

/** dynamics sub blocks **/
  if (which.compare(0,3,"dyn",0,3) == 0) {
    int pos = in.tellg();
    /** initial **/
    if (s1.compare(0,4,"init",0,4) == 0) {
      int n=-1;
      while(s.compare(0,3,"end",0,3)!=0) {
        in>>ws;
        in.getline(c,1000);
        s=c;
        n++;
      }
      in.seekg(pos);
      dyn.pop = new InitPop[n];
      for (int i=0; i<n; i++) {
        in>>ws;
        in.getline(c,1000);
        s=strtok(c," ");
        dyn.pop[i].mol = atoi(s.c_str());
        s=strtok(NULL," ");
        dyn.pop[i].state = atoi(s.c_str());
        s=strtok(NULL," ");
        dyn.pop[i].population = atof(s.c_str());
      }
      in.getline(c,1000);
    } //end initial population
  
    /** output **/
    if (s1.compare(0,6,"output",0,6) == 0) {
      int n=-1;
      dyn.popSet = false;
      dyn.printAll = true;
      while(s.compare(0,3,"end",0,3)!=0) {
        in>>ws;
        int pos2 = in.tellg();
        in.getline(c,1000);
        s=c;
        pos=in.tellg(); 
        char *c2 = c;
        if (s.compare(0,10,"populations",0,10)==0){
          dyn.popSet = true;
          char *cs;
          cs=strtok(c2," ");
          while(cs!=NULL) {
            cs=strtok(NULL," ");
            n++;
          }
          dyn.noutput = n;
          dyn.out = new Output[n];
          c2 = const_cast<char*>(s.c_str());
          s=strtok(c2," ");
          //s=strtok(NULL," ");
          if (n > 0) {
            for (int i=0; i<n; i++) {
              s=strtok(NULL," (,");
              cout<<s<<endl;
              if (s.compare(0,3,"all",0,3)!=0) {
                dyn.printAll = false;
              } else {
                cout<<"printing all populations"<<endl;
                continue;
              }
              dyn.out[i].mol = atoi(s.c_str());
              s=strtok(NULL,")");
              dyn.out[i].state = atoi(s.c_str());
            }
          
            //in.getline(c,1000);
            //s=strtok(c," ");
            //dyn.out[0].file = s;
          } else { 
            dyn.out = new Output[1];
          }
        } //end populations

        if (s.compare(0,4,"file",0,4)==0) {
          if (!dyn.popSet) {
            dyn.out = new Output[1];
          }
          c2 = const_cast<char*>(s.c_str());
          s=strtok(c," ");
          s=strtok(NULL," ");
          dyn.out[0].file = s;
        } //end file name
      }
    } //end output
  } //end mol block
}

/** read dynamics block **/

/** Constructor **/
Reader::Reader(string f) {
   /** Open COM file **/
   ifstream in;
   in.open(f.c_str());
   if (!in.is_open()) {
    cout<<"Cannot open file, quiting"<<endl;
    exit(0);
   }

   char c[1000];
   string s;
   int molcount = -1;
  while (in.good()) {
    /** Start reading **/
    in>>ws;
    in.getline(c,1000);
    s=c;
    
    /** skip lines that start with '#' **/
    if (s.compare(0,1,"#",0,1) == 0) {
      in>>ws;
      in.getline(c,1000);
      s=c;
    }
    
    /** check for main directive **/
    list<string>::iterator it = find(directives_.begin(),directives_.end(),s);
    if (s.compare(0,9,"molecules",0,9) == 0) {
      molcount++;
    }
    //found directive
    if (it!=directives_.end()) {
      readBlock(*it,in,molcount);
    }
  }

};

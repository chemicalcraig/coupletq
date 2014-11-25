/****************************************************
 * The reader class parses an exciton .com file
 * and sets up the appropriate calculation
 *
 * See README for appropriate format and 
 * options
 *
 * **************************************************/

class Reader() {
public:

char* filename;

string directiveList;
/** Calculation Stack **/
string type;
double ewindow;
int molecules;

/** Molecule Stack **/
int number;
int nstates;
int 

/** Constructor and Destructor **/
Reader() {};
~Reader() {};
};

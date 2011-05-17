/*
 *  simplex.H
 *  dees
 *
 *  Created by Yann Esposito on 17/04/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "general.H" // pour l'utilisation des map, hash_map...
#include "include/lp_lib.h"

typedef enum _simplexMode {realNotBounded, realZeroBounded, ZeroOrOne} simplexMode;

// this is a class used as a wrapper to other libraries
class Simplex {
private:
	int n;			// number of constraints
	lprec *lp;		// the linear programming variable (LP_SOLVE)
	double *tmpbuf; // temporary buffer used to fill the problem
	int tmpbufsize; // size of tmpbuf
	double val;     // temporary value for equations
	
	void resize_tmpbuf(void); // update the size of tmpbuf
	
public:
		Simplex(void) {
			lp = NULL;
			tmpbuf = NULL;
			tmpbufsize=0;
			n = 0;
			val=0;
		}
	~Simplex(void) {
		if (lp != NULL)
			free_lp(&lp);
		if (tmpbuf != NULL)
			free(tmpbuf);
	}
	
	// print the linear system
	inline void affiche(void) { print_lp(lp); }
	
	// set the mode for the variables ({0,1}, or Q+ or Q)
	int set_mode(const simplexMode newMode);
	
	// set the nuber of constraints
	inline lprec *set_number_of_unknown_to(const int newNumber) { 
		if (lp != NULL) 
			free(lp); 
		n = newNumber; 
		lp=make_lp(0,n);
		if (lp)
			set_verbose(lp,0);

		return (lp);
	}
	
	int setval(double newval);
	int setparam(int variable, double val);
	int add_absolute_constraint(); // use the current values setted by setval and setparam	
	int add_equality_constraint();
 	int add_equality_constraint_with_parameters(double *constraints, double val);

	// add a constraint of the form :
	// | constraints[0]X0 + ... + constraints[n-1]Xn-1 - val| <= epsilon
	int add_absolute_constraint_with_parameters(double *constraints, double val);
	int add_absolute_constraint_with_parameters(const vector<double> constraints, double val);
	int add_absolute_constraint_with_parameters(const map<int,double> constraints, double val);
	
	// return the solution int the map solution
	bool has_solution(map<int,float> &solution, float &epsilon);
	// renvoie 0 si 
	int set_minimize_epsilon_mode (double epsilon = -1);

};


/*
 *  simplex.cpp
 *  dees
 *
 *  Created by Yann Esposito on 17/04/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "simplex.H"

#include "include/lp_lib.h"
#include <math.h>

// set the mode for the variables ({0,1}, or Q+ or Q)
int Simplex::set_mode(const simplexMode newMode) {
	switch (newMode) {
		set_lowbo(lp,1,0.0); // pour epsilon c'est toujours 0
		set_upbo(lp,1,-log(0.0));    // pour epsilon c'est toujours 1
		set_int(lp,1,FALSE);      // toujours real
		case realNotBounded:
			for (int i=2 ; i<=n ; i++) {
				set_int(lp,i,FALSE);
				set_lowbo(lp,i,log(0.0));
				set_upbo(lp,i,-log(0.0));
			}
			break;
		case realZeroBounded:
			for (int i=2 ; i<=n ; i++) {
				set_int(lp,i,FALSE);
				set_lowbo(lp,i,0);
				set_upbo(lp,i,-log(0.0));
			}
			break;
		case ZeroOrOne:
			for (int i=2 ; i<=n ; i++) {
				set_binary(lp,i,TRUE);
			}
			break;
		default:
			cerr << "Simplex::set_mode : mode inconnu : " << newMode << endl;
			throw -1;
	}
	return VAL(0);
}
	

// add a equality constraint of the form :
// tmpbuf[1]X1 + tmpbuf[2]X2 + ... + tmpbuf[n]Xn  = val
// One can view X1 as an epsilon value
int Simplex::add_equality_constraint() {
	return add_equality_constraint_with_parameters(tmpbuf,val);
}


// add a equality constraint of the form :
// constraints[1]X1 + constraints[2]X2 + ... + constraints[n]Xn  = val
// One can view X1 as an epsilon value
int Simplex::add_equality_constraint_with_parameters(double *constraints,double val) {
	set_add_rowmode(lp,TRUE);
	add_constraint(lp,constraints,EQ,val);
	set_add_rowmode(lp,FALSE);
	return 0;
}


int Simplex::add_absolute_constraint(void) {
	return add_absolute_constraint_with_parameters(tmpbuf,val);
}
	
// add a constraint of the form :
// | constraints[2]X2 + ... + constraints[n]Xn - val | <= constraints[1]*X1
// One can view X1 as an epsilon value
int Simplex::add_absolute_constraint_with_parameters(double *constraints,double val) {
	set_add_rowmode(lp,TRUE);
	if (constraints[1] != 1) {
		cerr << "WARNING !!!! add_absolute_constratins_with_parameters\n";
		cerr << "constraints[n] should be 1 and it is " << constraints[n] << endl;		
	}
	// constraint of the form
	// constraints[1]X1 + ... + constraints[n]Xn >= -val
	add_constraint(lp,constraints,GE,val);
	// constraint of the form
	// -constraints[1]X1 + ... + constraints[n]Xn <= -val
	constraints[1] = -constraints[1];
	add_constraint(lp,constraints,LE,val);
	constraints[1] = -constraints[1];
	set_add_rowmode(lp,FALSE);
	return 0;
}

int Simplex::add_absolute_constraint_with_parameters(const vector<double> constraints, double val) {
	resize_tmpbuf();
	vector<double>::const_iterator i;
	int j;
	for (i=constraints.begin(), j=2 ; i!=constraints.end() ; i++, j++) {
		tmpbuf[j]=*i;
	}
	tmpbuf[1]=1;
	return add_absolute_constraint_with_parameters(tmpbuf,val);
}

int Simplex::add_absolute_constraint_with_parameters(const map<int,double> constraints, double val) {
	resize_tmpbuf();
	map<int,double>::const_iterator i;
	for (i=constraints.begin() ; i!=constraints.end() ; i++) {
		tmpbuf[i->first]=i->second;
	}
	tmpbuf[1]=1;
	return add_absolute_constraint_with_parameters(tmpbuf,val);
}

// set minimize epsilon mode or fixed epsilon mode
// if epsilon < 0 then minimize epsilon
// else fix the epsilon value
int Simplex::set_minimize_epsilon_mode (double epsilon) {
	resize_tmpbuf();
	for (int i=1;i<=n+1;i++) {
		tmpbuf[i]=0;
	}
	// tmpbuf = 0 0 ... 0
	if (epsilon < 0) { // minimize the value of epsilon
		tmpbuf[1]=1;
		set_obj_fn(lp,tmpbuf);
		set_minim(lp); // minimize the objective function
	} else {
		set_add_rowmode(lp,TRUE);
		set_obj_fn(lp,tmpbuf); // minimize the  0 function
		set_minim(lp);
		tmpbuf[1]=1;
		add_constraint(lp,tmpbuf,EQ,epsilon); // add the constraint epsilon = given value
		set_add_rowmode(lp,FALSE);
	}
	return VAL(0);
}


// répond oui ou non, il y a une solution
bool Simplex::has_solution(map<int,float> &solution, float &epsilon) {
	int res = solve(lp);
	// print_lp(lp);
	int temp ;
	if ((res == 0) ||
		(res == 1) ||
		(res == 3) ||
		(res == 4) ||
		(res == 12)  ) {
		temp = 2 + get_Nrows(lp) ;
		if (tmpbufsize < temp + n ) {
			if (tmpbufsize > 0) {
				free(tmpbuf);
			}
			tmpbuf = new double [ temp + n ];
			tmpbufsize = temp + n;
		}
		
		get_primal_solution(lp,tmpbuf);
		epsilon = tmpbuf[0];
		for (int i=1 ; i<= n ; i++) {
			solution[i]=(float)tmpbuf[temp + i -2];
		}
		return true;
	}
	else {
		solution.clear();
		return false;
	}
}


int Simplex::setval(double newval) {
	val = newval;
	return VAL(0);
}

void Simplex::resize_tmpbuf(void) {
	if (tmpbufsize < n + 2) {
		if (tmpbuf != NULL)
			free(tmpbuf);
		tmpbuf = new double[2*n+2];
		tmpbufsize=2*n+2;
	}
}

int Simplex::setparam(int variable, double val) {
	resize_tmpbuf();
	if (PFA_SAFE) {
		if ((variable<1) || (variable>n+1)) {
			cerr << "Aie Aie Aie !!! Simplex::setparam called with variable = " << variable << " where max should be ";
			cerr << n-1 << " !!! " << endl;
		}
	}
	tmpbuf[variable]=val;
	return VAL(0);
}


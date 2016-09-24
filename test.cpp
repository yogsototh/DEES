/*
 *  test.cpp
 *  dees
 *
 *  Created by Yann Esposito on 19/12/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "main.H"
#include "test.H"
#include "include/lp_lib.h"

void test_PSe(void) {
	MA A;
	SFunc PSe;
	SFunc::const_iterator q;
	int res;
	
	// A.becomeRandom(2,2);
	// A.save("test/testPSe.ma");
	A.load("test/testPSe.ma");
	
	res = A.val_PSe(PSe);
	cout << "res = " << res << endl;
	
	for (q=PSe.begin() ; q != PSe.end() ; q++) {
		cout << "P_" << q->first << "(\\Se) = " << q-> second << endl;
	}
}

void test_genere_mot_MA(void) {
    MA A;
	Word w;
    int i;
    A.load("test/testGenMot.ma");

    for (i=0;i<1;i++){
        cout << "i=" << i << " : ";
		A.genere_mot(w);
        cout << A.affiche(w);
        cout << endl;
    }
}


void test_simplex() {
	Simplex simplex;
	
	int i,j;
	int nbinc, nbcont;
	double input;
	
	cout << "Entrez le nombre d'inconnues : ";
	cin >> nbinc;
	simplex.set_number_of_unknown_to(nbinc);
	cout << "Entrez le nombres de contraintes : " ;
	cin >> nbcont;
	cout << "Entrez les contraintes (pour 3*X1 + 2*X2 = 4 écrire 3 2 4) : " << endl;
	
	for (i = 1 ; i<= nbcont ; i++) {
		for (j=1 ; j<= nbinc ; j++) {
			cin >> input;
			simplex.setparam(j,input);
		}
		cin >> input;
		simplex.setval(input);
		simplex.add_equality_constraint();
	}

	simplex.set_mode(realNotBounded);
	simplex.set_minimize_epsilon_mode(-1);
	
	map<int,float> sol;
	float epsilon;
	if (simplex.has_solution(sol,epsilon)) {
		cout << "Solution : " << endl;
		cout << "epsilon = " << epsilon << endl;
		for (i = 1 ; i<= nbinc ; i++) {
			cout << "X" << i << "=" << sol[i] << " ; ";
		}
		cout << endl;
	}
	else {
		cout << "Pas de solution..." << endl;
	}
	
}

void test(void) {
	test_simplex();
}

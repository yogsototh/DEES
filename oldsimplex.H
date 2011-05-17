#include "general.H"
#include <math.h> // pour fabs

// ------- La classe ----------
class Simplex
{
public:
    //hash_map<int, map <int, double> > a;  // La matrice des contraintes
    //hash_map<int, double> sol; // Le tableau des solutions
    //hash_map<int, int> db;  // Un tableau
    //hash_map<int, int> hb;  // un autre tableau

    map<int, map <int, double> > a;  // La matrice des contraintes
    map<int, double> sol; // Le tableau des solutions
    map<int, int> db;  // Un tableau
    map<int, int> hb;  // un autre tableau
    map<int, double>::iterator it; // un itérateur sur les colonnes de a.

public:
	Simplex() {}
    ~Simplex(){}

private:
    int ads_sortant(int m,int n,int phase);
    int ads_entrant(int n,int l);
    void pivotage(int m,int n,int l,int k);
	int affiche(int ineq1, int ineq2, int nbeq, int nbvar);
public:
    int  simplexe_dual(int ineq1, int ineq2, int eq, int n);
	int  solve(int ineq1, int ineq2, int eq, int n, bool debug=false);
};

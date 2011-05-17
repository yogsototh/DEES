/*
 *  pprfa_bm.cpp
 *  dees
 *
 *  Created by Yann Esposito on 18/04/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "pprfa.H"

RESULT PPRFA::becomeQuasiPrefixTree (
									 const Sample &S,   // les mots
									 const Word &v,     // le mot de retour
									 const SFunc &solution) // les valeurs initiales pour le retour
{
	Word w=v;
	Lettre a;
	SFunc::const_iterator s;
	Transition t;
	float val;
	
	w=v;
	a=*(w.rbegin());
	w.erase(--w.end());
	// on a v = wa
	
	// val = w^{-1}P(a\Se)
	val=exp(plog_bar(v) - plog_bar(w));
	
	
	// on supprime la transition deterministe phi(w,a,v)
	eraseTransition(w,a,v);
	
	// on fait les retours
	t.qdep = XR[w];
	t.a = a;
	for (s=solution.begin() ; s != solution.end() ; s++ ) {
		t.qarr = s->first;
		PFA::addTransition(t,s->second*val);
	}
	
	// on emmonde
	emmonde();
	
	put_prefix_tree(S,v);
	
	return VAL(0);
}



RESULT PPRFA::setRandomSolution(SFunc &solution, const StateSet &R) const
{
    StateSet::iterator q;
    for (q = R.begin() ; q != R.end() ; ++q)
    {
        solution[*q] = random(0,1);
    }
    float sum=0;
    SFunc::iterator s;
    for (s=solution.begin() ; s != solution.end() ; ++s)
    {
        sum += s->second;
    }
    for (s=solution.begin() ; s!=solution.end() ; ++s)
    {
        s->second /= sum;
    }
    return VAL(0);
}


RESULT PPRFA::choice(list<State> &X, const StateSet R) const {
	StateSet::const_iterator q;
	Graph::const_iterator r;
	LettreSFunc::const_iterator a;
	SFunc::const_iterator s;
	
	X.clear();
	for (q=R.begin() ; q!= R.end() ; q++) {	
		r=phi.find(*q);
		if (r != phi.end()) {
			for (a=r->second.begin() ; a != r->second.end() ; a++) {
				for (s=a->second.begin() ; s != a->second.end() ; s++) {
					if (R.find(s->first) == R.end()) {
						if (Q.find(s->first) != Q.end())
							X.push_back(s->first);
						else {
							cerr << "Choice !!! automate non coherent, ";
							cerr << "Etat : "<< s->first << endl;
							cerr << "Sauvegarde dans test/uncoherent" << endl;
							save("test/uncoherent");
						}
					}
				}
			}
		}
	}
	return VAL(0);
}

RESULT PPRFA::difference(const TransitionFunction T1, 
						 const TransitionFunction T2) const  {
	TransitionFunction::const_iterator t;
	float res=0;
	for (t=T1.begin() ; t != T1.end() ; t++) {
		res+=abs(t->second - T2.find(t->first)->second);
	}
	res = res/float(T1.size());
	return VAL(0);
						 }



// compte le nombre de fois que l'on passe par l'etat x en lisant S
float SPFA::count_nb_pass(const State x, const Sample &S) const {
	TransitionFunction T;
	SFunc Iota, Tau;
	
	// on met dans T toutes les transitions qui arrivent
	// dans l'etat x
	Graph::const_iterator q;
	LettreSFunc::const_iterator a;
	SFunc::const_iterator r;
	Transition t;
	t.qarr=x;
	for (q=phi.begin() ;q != phi.end() ; q++) {
		t.qdep = q->first;
		for (a = q->second.begin() ; a != q->second.end(); a ++) {
			t.a=a->first;
			if (a->second.find(x) != a->second.end()) {
				T[t]=0;
			}
		}
	}
	
	TransitionCount(S,T,Iota,Tau);
	TransitionFunction::const_iterator ti;
	float res;
	res = 0;
	for (ti=T.begin() ;  ti != T.end() ; ti++) {
		res+=ti->second;
	}
	return res;
}

// init the state solution
RESULT PPRFA::initSolution(SFunc &solution, const State x, const StateSet &R, const Sample &S) const {
	
	// a small security test
	if (! S.isprefixiel()) {
		cerr << "WARNING ! initSolution called with a not prefixialised Sample\n";
		return setRandomSolution(solution,R);
	}
	
	Word v;
	WordSet W;
	StateWordFunction wordof;
	StateSet::const_iterator r;
	Simplex simplx;		// simplx est juste l'objet Simplex qui permet de résoudre les problèmes
	int i;
	int it;	// un iterateur
	WordSet::const_iterator wr;
	StateWordFunction::const_iterator q;
	double s1;
	double s2;
	WordSet::const_iterator wi;				
	double pvw;	// v^{-1}P(w)		
	double precision;
	
	// --- Initialisation ---
	v=word_of(x); // v est le mot initial
	for (r=R.begin() ; r != R.end() ; r++) {
		wordof[*r]=word_of(*r); 
	}
	// wordof contient les les mots associes aux etats 
	ajoute_mots_associes(W,v,S,INT_MAX); // W contient les successeurs de v dans S
	
	// -- Remplissage du Simplex ---
	// number of variables = number of core states
	//n = R.size();
	
	simplx.set_number_of_unknown_to(R.size());
	simplx.setval(1);
	simplx.setparam(1,0);
	for (q = wordof.begin(),it=1; q != wordof.end(); q++)
		simplx.setparam(it++, 1);
	simplx.add_equality_constraint();
	
	simplx.setparam(1,1);
	// inequations en <= (nb=ineq1), inequations en >= (nb=ineq2)
	// les inequations X_q >= 0 sont deja
	// prisent en comptent de manière implicite
	s2=S[v];
	precision= 2/sqrt(double(s2));
	for (wi = W.begin(), i=1; wi != W.end () ; wi++, i++)
	{
		s1 = S[(v + *wi)];
		pvw = (double) s1 / (double) s2;
		simplx.setval(pvw);
		
		for (q = wordof.begin(),it=1; q != wordof.end(); q++, it++)
		{
			s1 = S[q->second + *wi];
			s2 = S[q->second];
			if (s2 == 0) {
				// cerr << "WARNING: " << q->second << ", random initialisation" << endl;
				return setRandomSolution(solution,R);
			}
			
			simplx.setparam(it, s1/s2);
		}
		simplx.add_absolute_constraint();
	}
	
	map<int,float> preciseSol;
	float epsilon;
	bool err = simplx.has_solution(preciseSol,epsilon);
	if ((! err) || (epsilon>precision))
	{	// pas de solutions
		setRandomSolution(solution,R);
		return ERR(0);
	}
	else
	{
		solution.clear ();
		for (q = wordof.begin(),it=2; q != wordof.end() ; q++, it++)
		{
			solution[q->first] = preciseSol[it];
			// cout << q->first << " SSS=" << simplx.sol[it]<<endl;
		}
		return VAL(0);
	}
}

// fonction qui teste si le retour semble plus avise
bool PPRFA::TestQPTA(const Sample &S,const State x,const StateSet R, const double precision, const double seuilbm, bool verbose) {
	PPRFA A;
	SFunc solution;
	WordSFunc::const_iterator r;
	Word v;
	Simplex simplx;
	StateSet::const_iterator s;
	PreciseSFunc F;
	int nb_tours;
	double nprime;
	
	// on recopie le PPRFA courant dans A
	A=*this;
	// on cherche le mot auquel correspond l'etat x
	save("test/tst");
	v=word_of(x);
	if (v == "") {
		cout << "compute_XR()" << endl;
		emmonde();
		compute_XR();
		v=word_of(x);
		cout << "v=" << affiche(v) << endl;
	}
	nprime=count_nb_pass(x,S);
	
	//setRandomSolution(solution,R);
	if (!S.isprefixiel()) {
		Sample S2;
		S2 = S;
		S2.prefixialise();
		initSolution(solution,x,R,S2);
	}
	else {
		initSolution(solution,x,R,S);
	}
	
	//	for (StateSet::const_iterator ri = R.begin(); ri != R.end() ; ri++) {
	//		cout << "sol[" << *ri << "]="<< solution[*ri] << endl;
	//	}
	
	A.becomeQuasiPrefixTree(S,v,solution);
	
	/* 
		// ----
	 int toto;
	 save("test/tst");
	 A.save("test/tstA");/*
	  cout << "lettre = " << v << endl;
	  cout << "entrez un chiffre pour continuer.";
	  cin >> toto;
	  // ----
	  */
	 
	 // on fait les BaumWelch
	 TransitionFunction T1,T2;
	 SFunc Iota,Tau;
	 
	 A.allTransitions(T1);
	 T2=T1;	
	 A.allStates(Tau);
	 
	 A.BaumWelch(S,T1,Iota,Tau,1, verbose);
	 A.BaumWelch(S,T2,Iota,Tau,1, verbose);
	 nb_tours=20;
	 while ((difference(T1,T2)>seuilbm) && (nb_tours-- > 0)) {
		 T2=T1;
		 A.BaumWelch(S,T1,Iota,Tau,1, verbose);
		 if (verbose) {
			 cout << "." << flush;
		 }
	 }
	 if (verbose) {
		 cout << "\n";
	 }
	 
	 if (verbose) {
		 cout << "validation de la solution (" << affiche(v) << "): ";
	 }
	 
	 // Test pour la validation de T
	 // ----
	 // if (verbose) {
	 //	 float x1, x2, x3, x4;
	 //	 x1 = A.Likelihood(S);
	 //	 x2 = S.AutoLikelihood();
	 //	 x3 = nprime;
	 //	 x4 = S.size();
	 //	 cout << "Hyp L=" << x1 << ", ";
	 //	 cout << "PTA L="<< x2 << ", ";
	 //	 cout << "nb(v)=" << x3 << ", ";
	 //	 cout << "Ssize=" << x4 << ", ";
	 //	 cout << "seuil=" << (x2-x1)*(sqrt(double(nprime))/double(S.size()));
	 //	 cout << endl;
	 // }
	 // ----
	 
	 float valtest=(abs( A.Likelihood(S) - S.AutoLikelihood() ) * sqrt(double(nprime)))  / double(S.size());
	 if (verbose) {
		 cout << valtest << " < " << precision << endl;
	 }
	 if ( valtest > precision) {
		 return false;
	 }
	 else {
		 *this=A;
		 return true;
	 }    
}

RESULT PPRFA::DEESBM(Sample &S,     // The learning sample (data)
                     const double prec,  // The precision parameter
                     const double epsprime, // The bound under which transition are erased
                     const bool verbose,  // verbose mode (show the states during construction)
                     const unsigned int maxstates, // The maximal states number after which learning is cancelled
                     const unsigned int seuil, // La periodicitÈ minimum avant de ne plus considÈrer l'Ètat
					 const double seuilbm, // la dist à partir de laquelle on considere avoir trouver le MLM
                     const unsigned int nb_tours ) // Le nombre de tour de BM
{
	
	try {
		list<State> X; // the set of considered words
		State x;
		StateSet R; // the set of traited states
		becomePrefixTree(S); // it's clear
		R.insert(iota.begin()->first);
		choice(X,R); // we add all states corresponding to letters
		while (!X.empty()) {
			x=*(X.begin()); X.pop_front();
			if ((Q.find(x) != Q.end()) && ( !TestQPTA(S,x,R, prec, seuilbm, verbose)) ) { // si on ne fait pas de retour
				R.insert(x);
				choice(X,R);
			}
		}
		erase_transitions(0.01,0);
		rend_PFA();
		return VAL(0);
	}
	catch (string e) {
		cerr << "DEESBM() ";
		cerr << e << endl;
		return ERR(0);
	}
}


RESULT PPRFA::get_solution(const PPRFA &A, 
						   SFunc &solution, 
						   Transition &t,
						   TransitionFunction &T,
						   const StateSFunc &Trad) 
{
	LettreSFunc::const_iterator a;
	SFunc::const_iterator r;
	SFunc::iterator s;
	float sum;
	a=A.phi.find(t.qdep)->second.find(t.a);
	
	for (r = a->second.begin() ; r != a->second.end() ; ++r) {
		t.qarr=r->first;
		solution[Trad.find(r->first)->second]=T[t];
	}
	
	sum = 0;
	for (s=solution.begin() ; s != solution.end() ; ++s) {
		sum += s->second;
	}
	if (sum == 0) {
		throw 1;
	}
	
	for (s = solution.begin() ; s != solution.end() ; ++s) {
		s->second /= sum;		
	}
	
	return VAL(0);
	
	// ----
	for (s=solution.begin() ; s != solution.end() ; ++s) {
		cout << "sol[" << s->first << "] = " << s->second << " ";
	}
	cout << endl;
	// ----
	
}


/*
 *  pprfa_dees.cpp
 *  dees
 *
 *  Created by Yann Esposito on 18/04/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "pprfa.H"

// ================================================================
// ===                                                          ===
// ===                          DEES                            ===
// ===                                                          ===
// ================================================================

// Methode d'apprentissage d'un PRFA prefixe
RESULT
PPRFA::DEES (T_ModeVariables modeVariables,
			 Sample & S,
			 double prec,
			 double epsprime,
			 bool verbose,
			 T_ModeReturn moderet,
			 T_ModeEpsilon modeeps,
			 unsigned int maxstates,
			 unsigned int seuil, 
			 int maxmots, 
			 int maxsearches,
			 bool bestsearch,
			 bool stepssave)
{
    try {
	// ------------------------- affichage des informations ---------------------- //
	if (verbose)
	{
		cout << "\t\t=== DEES ===";
		cout << "\nSample size = " << S.size ();
		cout << "\nprecision = " << prec;
		cout << "\nbound = " << epsprime;
		cout << "\nmoderet = ";
		switch (moderet)
		{
            case begin:
                cout << "begin";
                break;
            case end:
                cout << "end";
		}
		cout << "\nmodeeps = ";
		switch (modeeps)
		{
            case epsfixed:
                cout << "epsfixed";
                break;
            case variable:
                cout << "variable";
                break;
            case word_variable:
                cout << "word_variable";
                break;
		}
		cout << "\nmaxstates = " << maxstates;
		cout << "\nseuil = " << seuil;
		cout << "\nmaxmots = " << maxmots << endl;
	}
	
	if (prec == 0) {
		return becomePrefixTree(S);
	}
	
	// -------------------------- INITIALISATION ----------------------
	vide ();		// on initialise le PFA.
	if (S.size () == 0)
		throw 0;		// on verifie que l'echantillon n'est pas vide.
	S.prefixialise ();	// Transformation de l'Èchantillon en echantillon prefixiel.
	Sigma = S.alphabet ();	// on met ‡† jour l'alphabet
	alph = S.dictionnaire ();	// et le dictionnaire associÈ.
	
	// Declaration
	Word v;			// the word v which represent the word associated to the state to add.
	list < Word > X;		// The set of words which are potential prime residuals.
	Sample::const_iterator u;	// current word of the sample.
	float val;		// a floating value used to calculate tau.
	Word w;			// a word
	Word ww;			// another word
	Lettre a;		// a letter (we will have v=wa)
	Alphabet::const_iterator b;	// pour enumerer toutes les lettres
	SFunc solution;	// the system solution
	SFunc solutiontemp;	// a temporary solution
	StateSet::iterator q;	// Pour enumerer les √©tats
	Simplex simplx;		// L'objet simplexe qui contient les algorithmes de resolution de systemes.
	set < Word, ordre_mot > W;	// L'ensemble de mots à tester
	Word::iterator z;	// last letter
	
	// --- init variables ---
	v.clear ();		// v = epsilon (empty word)
					// X is the set of one letter words of S when prefixialised
	val = 0;
	for (u = ++(S.begin ());
		 (u != S.end ()) && (u->first.size () < 2); ++u)
	{
		X.push_back (u->first);
		val += u->second;
	}
	// We create the first state (epsilon)
	addState (v, 1, 1 - (float(val) / float(S.size ())));
	W.clear ();
	
	// liste_mots_associes(W,v,S,maxmots);
	ajoute_mots_associes(W,v,S,maxmots);
	
	// There may be a general state
	State generalstate = -1;
	
	// Step saving
	string svgfile="etape-";
	
	// -------------------------- LEARNING LOOP  -----------------------
	if (verbose)
		cout << "Ajout des états : " << endl;
	
	// For each element in X (the temporary list)
	while (!X.empty ())
	{
		v = *(X.begin ());	// v <-- min X;
		X.pop_front ();	// X <-- X\{v} ;
						// wa=v
		w = v;
		a = *(w.rbegin ());
		z = w.end ();
		--z;
		w.erase (z);	//w.pop_back ();
		
		//~ if (verbose) {
		//~ cout << "[";
		//~ cout << affiche(v) <<"]";
		//~ cout.flush();
		//~ }
		
		if (stepssave) {
			save(svgfile+affiche(v));
		}
		
		
		/// 3 possible cases : 
		///     (1) not enought data (make a loop to the general state)
		///     (2) there is no solution (add a new state and update X)
		///     (3) there is a solution (make a return of transitions)
		
		if (S[v] < seuil) // case (1) not enought data
		{
			// cout << "CASE (1)" << endl;
			if (generalstate == -1)
			{		// if the general state was not created, we create it
				generalstate = addState (v, 0, 1 / float(Sigma.size () + 1));
				for (b = Sigma.begin (); b != Sigma.end (); ++b)
				{
					MA::addTransition (generalstate, *b, generalstate, 1 / float(Sigma.size () + 1));
				}
			}
			
			addTransition (w, a, generalstate, ((double) S[v] / (double) S[w]));
			XR[w + a] = generalstate;
			if (verbose)
			{
				cout << "S";
				cout.flush ();
			}
		}
		else
		{
			solution.clear ();	// init solutions
			// calculate the solution of the linear system
			sol (modeVariables, solution, v, S, simplx, prec / pow(float(S[v]), float(0.4)), moderet, modeeps, W, maxmots, false);	
			if (solution.empty ()) // if there is no solution (case (2) add a new state and update X
			{
				// cout << "CASE (2)" << endl;

				// if there is no solution then we add the state associated with v
				// updating X and tau (val will contain tau[v])
				val = 0;
				for (b = Sigma.begin (); b != Sigma.end (); b++)
				{		// pour toute les lettres de l'alphabet
					v += *b;
					//v.push_back (*b);     // on ajoute b a la fin de v
					if (S.find (v) != S.end ())
					{	// si vb appartient a l'echantillon, alors
						X.push_back (v);	// on l'ajoute a X
						val += S[v];	// et val = val + S(vb\Sigma^*)
					}
					z = v.end ();
					--z;
					v.erase (z);	//v.pop_back ();  // on efface la derniere lettre pour que la variable v = le mot v.
				}
				
				if (verbose)
				{
					cout << "A";
					cout.flush ();
				}
				
				addState (v,0, 1 - (val / (double) S[v]));	// adding the new state
				if (size () > maxstates)
					throw 1;
				addTransition (w, a, v, ((double) S[v] / (double) S[w]));	// updating phi : phi(w,a,wa) = S(wa)/S(w)
				
			}
			else
			{
				// cout << "CASE (3)" << endl;

				// else we return transitions
				
				if (verbose)
				{
					cout << "R";
					cout.flush ();
				}
				
				for (q = Q.begin (); q != Q.end (); q++)
				{
					val = (float) solution[*q] *
					((double) S[v] / (double) S[w]);
					if (val != 0)
					{
						addTransition (w, a, *q, val);
					}
				}
			}
		}
	}
	if (verbose) 
		cout << endl;
	
	erase_transitions (epsprime, -epsprime);	// deleting transitions less than epsprime
									//best_transition_delete(S,verbose);
	if (modeVariables == nonconstrained) {
		return VAL(0);
	}
	else {
		return renormalise ();	// renormalisation in order to disable borders effects
	}
}
catch (int erreur) {
	if (PFA_VERBOSE)
	{
		if (erreur != 1)
		{
			cerr << "PFA::DEES(Sample S)" << endl;
		}
		switch (erreur)
		{
            case 0:
                cerr << "Sample vide !!!" << endl;
                break;
            case 1:		// il y a trop d'etats
                erase_transitions (epsprime);
                return renormalise ();
                break;
            default:
                cerr << "Erreur n∞" << erreur << endl;
		}
	}
	return ERR (erreur);
}
}



RESULT PPRFA::ajoute_mots_associes (set <Word, ordre_mot> &W,
                                    const Word & v,
                                    const Sample & S,
                                    int maxmots) const
{
    list < Word > L;		// une liste de mot qui sert pour calculer W
    Alphabet::const_iterator b;	// pour enumerer toutes les lettres
    Word w;				// Le mot courant
	
	
    w.clear();			// on initialise avec w <-- epsilon
						// On ajoute tous les successeurs de v
    L.push_back (w);		// et L = { eps }
    while ((!L.empty ()) && (maxmots != 0))
    {				// tant que L n'est pas vide
        w = L.front ();
        L.pop_front ();		// w = min(L) et L=L\{w}
		if (W.find(w) == W.end()) {
			--maxmots;
			W.insert (w);
		}
        for (b = Sigma.begin (); b != Sigma.end (); b++)
        {
            w += *b;  // w <-- wb
			//W.insert(w);
            if (S.find (v+w) != S.end ()) // si p_n(w)>0
            {
                L.push_back (w);
            }
            w.erase (--w.end());	//w.pop_back ();  // wb <-- w
        }
    }
	
    return VAL (0);
}

// Return the set of word used to construct the inequation system
// (the set W described in the paper in the definition of system I)
RESULT
PPRFA::liste_mots_associes (WordSet &W,
                            const Word & v,
                            const Sample & S,
                            int maxmots) const
{
    map<Word, State>::const_iterator q;	// L'etat courant
    W.clear ();		// on initialise W
					// On ajoute tous les successeurs de v
    ajoute_mots_associes (W, v, S);
    // On ajoute tous les successeurs des elements de R
    for (q = XR.begin (); q != XR.end (); q++)
    {
        ajoute_mots_associes (W, q->first, S, maxmots);
    }
    return VAL (0);
}

// renvoie vide si le systËme lineaire I(v,R,S,precision) n'a pas de solution et
// sinon renvoie une solution de I i.e. un  ensemble de variables X_q telles que
// v^{-1}P = \sum_{q \in R} X_qP_q --  en rappelant que R[w]=q => P_q = w^{-1}P
// En prenant un maximum de max variables
RESULT PPRFA::solmax (T_ModeVariables modeVariables,
					  SFunc & solution,	// La solution
                      const Word & v,	// Le residuel a tester
                      const Sample & S,	// L'echantillon
                      Simplex & simplx,	// L'objet simplex qui permet de calculer la solution
                      double precision,	// La precision (sa signification depend du mode)
                      T_ModeReturn moderet,	// Le mode de retour : beg debut de l'arbre, end √† la fin
                      T_ModeEpsilon modeeps,	// Le mode : epsfixed avec epsilon fixe, variable avec epsilon variable.
                      WordSet &W,	// L'ensemble de mots sur lesquels on teste,
									// doit contenirs les mots des √©tapes prÈcÈdentes
                      unsigned int max) const
{
    try
    {
		
        //    cout << "\t\tcreation du systeme" << endl;
		
        unsigned int i;
        int n;		// number of variables
        int err;
        WordSFunc::const_iterator q;
		
        // number of variables = number of states
        n = XR.size ();
		
		/*
		// ---------- Affichage de la liste de mots utilises --------------
          set<Word, ordre_mot>::const_iterator wy;
          map<Word, State>::const_iterator ri;
		  int compte = 1;
          cout << "\nliste de mots associes à " << v << " et à {";
          for (ri=XR.begin() ; ri != XR.end() ; ri++) {
            if (ri->first.empty())
				cout << "eps,";
            else
				cout << ri->first << ",";
          }
          cout << "}" << endl;
          for (wy = W.begin() ;  wy != W.end() ;wy++) {
            if (wy->empty())
				cout << "eps,";
            else
				cout << *wy << ",";
			if (!(compte++ % 10))
				cout << "\n";
          }
		  cout << "\nNombre de mots : " << W.size();
          cout << endl;
          // --------------------------------------------------------------
		*/
		
        // avec l'entree 0 on traite tous les mots
        if (max == 0)
            max = W.size () + 1;
		
        // inequations en <= (nb=ineq1), inequations en >= (nb=ineq2)
        // les inequations X_q >= 0 sont deja
        // prisent en comptent de mani√®re implicite
        double s1;
        double s2;
        set< Word, ordre_mot >::const_iterator wi;
		
        double pvw = 0;	// v^{-1}P(w)
		
        s2 = S[v];
		
		// Initialisation du simplexe
		simplx.set_number_of_unknown_to(XR.size()+1);
		
		// sum of all Xi is equal to 1
		simplx.setparam(1,0); // --- la variable 1 est toujours epsilon ---
		simplx.setval(1);     //
        for (q = XR.begin (); q != XR.end (); q++)
            simplx.setparam(q->second+1,1);
		simplx.add_equality_constraint();
		// cout << "XR.size() = " << XR.size() << endl;
		
		
		// loop for each word :
		simplx.setparam(1,1); // --- la variable 1 est toujours epsilon ---
				
		i=0;
		for (wi = W.begin (); wi != W.end () && i <= max; wi++)
		{
			s1 = S[(v + *wi)];
			pvw = (double) s1 / (double) s2;
			simplx.setval(pvw);
			for (q = XR.begin (); q != XR.end (); q++)
			{
				s1 = S[q->first + *wi];
				simplx.setparam(q->second+1, (double) s1 / (double) S[q->first]);
			}
			simplx.add_absolute_constraint();
			++i;
		}
		
		
        switch (modeeps)
        {
            // ---------------------------------------------------------- EPSFIXED -
			case epsfixed:
				simplx.set_minimize_epsilon_mode(precision);
				break;
			// ---------------------------------------------------------- VARIABLE -
			case variable:
				simplx.set_minimize_epsilon_mode(-1);
				break;
			default:
				cerr << "!!! PFA::sol -- modeeps " << modeeps << " inconnu !!!" << endl;
        }
		
		switch(modeVariables) {
			case determinist:
				simplx.set_mode(ZeroOrOne);
				break;
			case positive:
				simplx.set_mode(realZeroBounded);
				break;
			case nonconstrained:
				simplx.set_mode(realNotBounded);
				break;
		}
		
		map<int,float> preciseSol;
		float epsilon;
		// simplx.affiche();
		err = simplx.has_solution(preciseSol, epsilon);
		// cout << (err?"solution trouvée":"pas de solution") << ", epsilon = " << epsilon << endl;
		if ((!err) || (epsilon>precision)) {
			solution.clear();
			return ERR(0);
		}
		else {
			solution.clear();
			float tmp;
			for (map<int,float>::const_iterator q =  preciseSol.begin() ; q != preciseSol.end() ; q++) {
				if (q->first != 1) {
					if ((tmp = (float) q->second) != 0) {
						solution[q->first - 1] = (float) q->second;
					}
				}
			}
			return VAL(0);
		}
    }
    catch (int erreur)
    {
        if (PFA_VERBOSE)
        {
            cerr << "ERREUR !!! PFA::solmax !!! n¬∞" << erreur << endl;
        }
        solution.clear ();
        return erreur;
    }
}




// renvoie vide si le syst√®me lineaire I(v,R,S,precision) n'a pas de solution et
// sinon renvoie une solution de I i.e. un  ensemble de variables X_q telles que
// v^{-1}P = \sum_{q \in R} X_qP_q --  en rappelant que R[w]=q => P_q = w^{-1}P
RESULT PPRFA::sol (T_ModeVariables modeVariables, // le mode dans lequel renvoyer les variables
				   SFunc & solution,	// La solution
                   const Word & v,	// Le residuel a tester
                   const Sample & S,	// L'echantillon
                   Simplex & simplx,	// L'objet simplex qui permet de calculer la solution
                   double precision,	// La precision (sa signification d√©pend du mode)
                   T_ModeReturn moderet,	// Le mode de retour : beg debut de l'arbre, end √† la fin
                   T_ModeEpsilon modeeps,	// Le mode : epsfixed avec epsilon fixe, variable avec epsilon variable.
                   WordSet &W,	// L'ensemble de mots sur lesquels on teste,
								// doit contenirs les mots des √©tapes pr√©c√©dentes
                   int maxmots,	// le nombre de mots qu'on ajoute au maximum
                   bool Wcalculated	// true if W is considered as calculated from outside
				   ) const
{
    // unsigned int n = 2 * Q.size () + 1;
	
    // - Pour connaitre le nombre d'inequation, on doit
    // - connaitre le nombre de mots qui sont successurs de v, et de tous
    // - les ÈlÈments de R.
    // - On met l'ensemble des successeurs de v et des √©l√©ments de R de l'echantillon dans W.
	
    if (!Wcalculated)
    {
        ajoute_mots_associes (W, v, S, maxmots);
    }
	//liste_mots_associes(W,v,S,(int)((float)S.size()/(float)XR.size()));
	return solmax (modeVariables,solution, v, S, simplx,precision, moderet, modeeps, W, INT_MAX);
	
/*	
    solret = solmax (modeVariables, solution, v, S, simplx, precision, moderet,modeeps, W, n);
    if ( !solret || (n >= W.size ()))
        return err;
    else
    {
        n *= 5;
        solret = solmax (modeVariables, solution, v, S, simplx, precision,moderet, modeeps, W, n);
        if (!solret || (n >= W.size ()))
            return err;
        else
            return solmax (modeVariables,solution, v, S, simplx,precision, moderet, modeeps, W, INT_MAX);
    }
 */
}
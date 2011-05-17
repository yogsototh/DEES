/***************************************************************************
                                 pfa.cpp 
                             -------------------
    begin                : Mon 9 Dec 2002
    copyright            : (C) 2002 by Yann Esposito
    email                : esposito@cmi.univ-mrs.fr
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/*
  Pour les notes, voir le fichier "pfa.h".
  See "pfa.h" file to view notes.
*/
#include "pfa.H"
//#include <list>
#include <assert.h>

// renvoi vrai si le MA est un PFA
RESULT
PFA::rend_PFA (void)
{
	// on emmonde
	MA::emmonde();
    // renormalisation of all parameters
    renormalise();
    return VAL (0);
}

RESULT
PFA::becomeRandom (const int nb_etats, const int nb_lettres,
                   const int num_graphe, const float densite,
                   const float prob_init, const float prob_term,
                   const float min_trans, const float max_trans)
{
    // on cree un SPFA
    SPFA::becomeRandom (nb_etats, nb_lettres,
                        num_graphe, densite,
                        prob_init, prob_term,
                        min_trans, max_trans);
    // puis on l'emmonde et on le renormalise.
    return rend_PFA ();
}

// Cree un MA aleatoire avec un nombre maximal d'aretes
RESULT PFA::becomeRandomMax (const int nb_etats, const int nb_lettres,
                             const int num_graphe, const int nb_succ,
                             const int nb_init,    const int nb_term,
                             const float min_trans,  const float max_trans)
{
    SPFA::becomeRandomMax(nb_etats, nb_lettres, num_graphe, nb_succ, nb_init, nb_term,min_trans,max_trans);
    return rend_PFA();
}


RESULT
PFA::becomeRandomPrefix (const int nb_etats, const int nb_lettres,
                         const int num_graphe, const float densite, const float prob_init,
                         const float prob_term, const float min_trans,
                         const float max_trans)
{

    // On gÈnËre un MA alÈatoire avec seulement des arÍtes positives
    SPFA::becomeRandomPrefix (nb_etats, nb_lettres,
                              num_graphe, densite,
                              prob_init, prob_term,
                              min_trans, max_trans);

    // on rend un SPFA
    return rend_PFA ();
}

// Arbre prefixe
RESULT PFA::becomePrefixTree(Sample &S)
{
    bool ispref = S.isprefixiel();
    if (!ispref)
    {
        S.prefixialise();
    }

    vide(); // on vide la structure du PFA
    Sample::const_iterator w;
    Lettre a;
    Word v;

    Sigma=S.alphabet();
    alph=S.dictionnaire();

    map< Word, State > R;
    w = S.begin();
    R[w->first]=addNewState(1,0);

    // on commence apres epsilon d'ou le ++w au debut du for
    for (++w ; w != S.end() ; ++w)
    {
        v=w->first;
        R[v]=addNewState(0,0);
        a=*(v.rbegin());
        v.erase(--v.end());
        addTransition(R[v], a , R[w->first],(float)S[w->first]/(float)S[v]);
    }
    // on met ‡ jour iota et tau
    StateSet::const_iterator e;
    Graph::const_iterator q;
    LettreSFunc::const_iterator x;
    SFunc::const_iterator r;
    float sum;
    for (e=Q.begin() ; e != Q.end() ; ++e)
    {
        sum = 0;
        q=phi.find(*e);
        if (q != phi.end())
        {
            for (x=q->second.begin() ; x != q->second.end() ; ++x)
            {
                for (r=x->second.begin() ; r != x->second.end() ; ++r)
                {
                    sum += r->second;
                }
            }
        }
        tau[*e]=1-sum;
    }

    if (!ispref)
    {
        S.deprefixialise();
    }
    return VAL(0);
}



// Lissage
RESULT PFA::lisse (int mode, float delta)
{
    if (mode == 0)
    {				// Mode sale !!!
        SFunc::iterator q;
        Alphabet::const_iterator a;
        float val;
        State s;			// le nouvel etat

        if (iota.empty ())
        {
            s = addNewState(1,delta);
        }
        else
        {
            q = iota.begin ();
            val = q->second * delta;
            q->second = q->second * (1 - delta);
            s = addNewState (val, delta);
        }
        val = 1 - delta;
        for (a = Sigma.begin (); a != Sigma.end (); ++a)
            addTransition (s, *a, s, val / (Sigma.size ()));
		
		return VAL(0);
    }
    else if (mode == 1)
    {
        cerr << "PFA::lisse(), mode " << mode << " non pris en compte"
        << endl;
        return ERR (1);
    }
    else
    {
        cerr << "PFA::lisse(), mode" << mode << " non pris en compte" << endl;
        return ERR (0);
    }
}





// Vraissemblance
float PFA::Likelihood(const Sample &S)
{	
	if (S.isprefixiel()) {
		cerr << "PFA::Likelihood(const Sample &S) : S is prefixial !" << endl;
		return -1;
	}

    Sample::const_iterator w;
    float result=0;
//	cout.precision(19);
    for (w=S.begin() ; w != S.end() ; ++w)
    {
		if (isnan(plog(w->first))) {
			cerr << "ERREUR mot : " << w->first << endl;
		}
		else {
			result += plog(w->first) * w->second;
		}
    }

    return result;
}


// genere un mot ‡ l'aide du PFA.
RESULT
PFA::genere_mot (Word & w, const bool safe)
{
    if (PFA_SAFE)
    {
        if (safe)
            assert (isPFA ());
    }

    State jeton;
    float alea;
    float sum;
    SFunc::const_iterator q;
    LettreSFunc::const_iterator a;

    // w = mot vide
    w.clear ();

    // choix de l'Ètat initial q
    // on prend un nombre de [0,1[
    alea = ((float) rand ()) / INT_MAX;
    // on prend le premier etat qui fait dÈpasser le seuil
    q = iota.begin ();
    sum = q->second;
    while (sum < alea)
    {
        sum += (++q)->second;
    }
    jeton = q->first;	// on a fait notre choix

    // Fin du choix : on commence dans l'Ètat jeton
    // on genere le mot
    // alea prend une valeur au hasard de [0,1[, puis tant que l'on ne s'arrete pas
    // i.e. tant que alea > tau(jeton),
    alea = ((float) rand ()) / INT_MAX;
    while (alea >
            (sum =
                 (tau.find (jeton) !=
                  tau.end ())? (tau.find (jeton)->second) : 0))
    {
        // on choisi une transition ‡ emprunter
        a = (phi.find (jeton)->second).begin ();	//On n'a pas besoin de faire un test sur le find !
        while (sum < alea)
        {
            q = a->second.begin ();
            while ((q != a->second.end ()) && (sum < alea))
            {
                if ((sum += q->second) < alea)
                    q++;
            }
            if (sum < alea)
                a++;
        }
        // et on l'emprunte
        jeton = q->first;
        // on rajoute la lettre gÈnÈrÈe ‡ w
        w += a->first;
        //w.push_back (a->first);

        // on est dans l'Ètat "jeton" et on a ajoute la lettre "a->first"
        //puis on recommence
        alea = ((float) rand ()) / INT_MAX;
    }
    return VAL (0);
}


 // genere un Èchantillon de mots
 RESULT
 PFA::genere_echantillon (int taille, Sample & S,
                          const int num_echantillon)
 {
	 if (isPFA()) {
		 int i;
		 Word w;
		 if (num_echantillon != 0)
		 {
			 srand (num_echantillon);
		 }
 
		 S = Sample (Sigma, alph, false);
		 for (i = 0; i < taille; i++)
		 {
			 genere_mot (w, false);
			 S.insert (w);
		 }
		 return VAL (0);
	 }
	 else {
		 return MA::genere_echantillon(taille, S, num_echantillon);
	 }
 }
 
 // genere un Èchantillon de mots : l'echantillon reste identique si lot est identique
 RESULT PFA::genere_echantillon (const int taille,
                                 Sample & S,
                                 const char *Filename,
                                 const int num_echantillon)
 {
     if (PFA_DEBUG)
     {
         if (Filename == NULL)
         {
             cerr << " PFA::echantillon(int taille, char *Filename) : Filename = NULL !!!" << endl;
         }
     }
     genere_echantillon (taille, S, num_echantillon);
     return S.save (Filename);
 }

// ---------- DISTANCES ---------

// renvoie la perplexité du PFA par rapport à l'echantillon S
double
PFA::perplexite (const Sample & S) const
{
    double res, sum;
    Sample::const_iterator w;

	if (isPFA()) {
		res = sum = 0;
		for (w = S.begin (); w != S.end (); S.next (w))
		{
			res += plog (w->first) * w->second;
			sum += w->second;
		}
	}
	else {
		res = sum = 0;
		for (w = S.begin (); w != S.end (); S.next (w))
		{
			res += log (r(w->first)) * w->second;
			sum += w->second;
		}		
	}
    return -res / sum;
}

// renvoie la perplexitÈ du PFA par rappor ‡ l'echantillon S
double
PFA::divKL (PFA &B, const Sample & S) const
{
    double res, res2, sum;
    Sample::const_iterator w;

    res = res2 = sum = 0;
    for (w = S.begin (); w != S.end (); S.next (w))
    {
		res += plog (w->first) * w->second;
		res2 += B.plog(w->first) * w->second;
		sum += w->second;
    }
    return (res2 - res) / sum;
}


// renvoie la norme 2 entre ce PFA et le PTA correspondant ‡ S
double
PFA::d_2(const Sample &S) const {
    Sample::const_iterator w;
    double res, sum;
    double pA, pB;
	double taille=double(S.size());

    res = sum = 0;
    for (w = S.begin (); w != S.end (); w++)
    {
        if (w->second >= S.seuil)
        {
            pA = p (w->first);
            pB = double(w->second)/taille;
            res += (pA - pB) * (pA - pB) * pA;
            sum += pA;
        }
    }

    return sqrt (res / sum);
}

// renvoie la norme 2 entre ce PFA et un autre
double
PFA::d_2 (const PFA & B, const Sample & S) const
{
    Sample::const_iterator w;
    double res, sum;
    double pA, pB;

    res = sum = 0;

    for (w = S.begin (); w != S.end (); w++)
    {
        if (w->second >= S.seuil)
        {
            pA = p (w->first);
            pB = B.p (w->first);
            res += (pA - pB) * (pA - pB) * pA;
            sum += pA;
        }
    }

    return sqrt (res / sum);
}

// renvoie la norme L1 entre ce PFA et un autre
double
PFA::L_1 (const PFA & B, const Sample & S) const
{
    Sample::const_iterator w;
    double res, sum;
    double pA, pB;

    res = sum = 0;
    for (w = S.begin (); w != S.end (); w++)
    {
        if (w->second >= S.seuil)
        {
            pA = p (w->first);
            pB = B.p (w->first);
            res += (pA - pB >= 0) ? (pA - pB) * pA : (pB - pA) * pA;	// |pA(w) - pB(w)|*pA(w)
            sum += pA;
        }
    }

    return sqrt (res / sum);
}


// renvoie la norme 2 entre ce PFA et un autre
double
PFA::d_2nlog (const PFA & B, const Sample & S) const
{
    Sample::const_iterator w;
    double res, sum;
    double pA, pB;
	
    res = sum = 0;
	
    for (w = S.begin (); w != S.end (); w++)
    {
        if (w->second >= S.seuil)
        {
            pA = p_directe (w->first);
            pB = B.p_directe (w->first);
            res += (pA - pB) * (pA - pB) * pA;
            sum += pA;
        }
    }
	
    return sqrt (res / sum);
}

// renvoie la norme L1 entre ce PFA et un autre
double
PFA::L_1nlog (const PFA & B, const Sample & S) const
{
    Sample::const_iterator w;
    double res, sum;
    double pA, pB;
	
    res = sum = 0;
    for (w = S.begin (); w != S.end (); w++)
    {
        if (w->second >= S.seuil)
        {
            pA = p_directe (w->first);
            pB = B.p_directe (w->first);
            res += (pA - pB >= 0) ? (pA - pB) * pA : (pB - pA) * pA;	// |pA(w) - pB(w)|*pA(w)
            sum += pA;
        }
    }
	
    return sqrt (res / sum);
}



// renvoie la moyenne des Ècarts de log en fonction d'un Èchantillon S
double
PFA::dlog (const PFA & B, const Sample & S) const
{
    Sample::const_iterator w;
    double res, sum, pa;
    res = sum = 0;
    for (w = S.begin (); w != S.end (); w++)
    {
        sum += w->second;
        pa = plog (w->first) - B.plog (w->first);
        if (pa < 0)
            pa = -pa;
        res += pa * w->second;
    }
    return res / sum;
}

// renvoie la moyenne des Ècarts de log en fonction d'un Èchantillon S
double
PFA::dlog_safe (const PFA & B, const Sample & S) const
{
    Sample::const_iterator w;
    double res, pa;
    int sum, nb_err, nb_mots;

    res = sum = nb_err = nb_mots = 0;
    for (w = S.begin (); w != S.end (); w++)
    {
        if (w->second >= S.seuil)
        {
            nb_mots++;
            pa = plog (w->first) - B.plog (w->first);
            if (!isnan (pa))
            {
                sum += w->second;
                if (pa < 0)
                    pa = -pa;
                res += pa * w->second;
            }
            else
            {
                if (PFA_VERBOSE)
                {
                    cerr << "\tErreur : ";
                    Word::const_iterator va;
                    for (va = w->first.begin ();
                            va != w->first.end (); va++)
                        cerr << alph.find (*va)->
                        second << " ";
                    cerr << endl;
                    cerr << "\tpA = " << plog (w->
                                               first) <<
                    " et pB = " << B.plog (w->
                                           first)
                    << endl;
                }
                nb_err++;
            }
        }
    }
    if (PFA_VERBOSE)
    {
        if (nb_err > 0)
            cerr << nb_err * 100 /
            nb_mots << "% d'erreur de structure." << endl;
    }
    return res / sum;
}

// renvoie la divergence de Kullback Leibler relativemant ‡ l'Èchantillon S
double
PFA::divKL (const PFA & B, const Sample & S)
{
    return B.perplexite (S) - perplexite (S);
}

// return true if it is a PRFA
bool
PFA::isPRFA (bool verbose)
{
    if (!isPFA())
    {
        return false;
    }

    // !!!!!!!!!!! VERSION NON DEFINITIVE !!!!!!!!!!!!!!
    // !!!!! reduce is not yet implemented !!!!!!!!!

    //reduce();

    map < State, bool > resStates; // resStates[q]=true iff q is reach uniquely by some word.
    set
        < StateSet> Memory; // the set of states even founded.
    StateSet QI; // l'ensemble des Ètats initiaux
    StateSet R; // Un ensemble d'Ètats
    StateSet Qm; // l'ensemble des Ètats dont on a dÈja trouvÈ un mot caractÈristique.
    set
        < StateSet >::const_iterator Ri;
    set
        < Word > W; // un ensemble de mots
    Word w;
    Word::iterator b;
    StateSet::const_iterator q;
    SFunc::iterator r;
    Alphabet::const_iterator a;

    //~ // -- //
    //~ cout << "L'ensemble d'Ètats est : {";
    //~ StateSet::const_iterator s;
    //~ for (s=Q.begin() ; s != Q.end() ; ++s)
    //~ cout << *s << ",";
    //~ cout << "}" << endl;
    //~ // -- //

    // QI devient l'ensemble des Ètats initiaux
    for (r=iota.begin();r != iota.end() ; ++r)
    {
        QI.insert(r->first);
    }

    //~ // -- //
    //~ cout << "QI={";
    //~ for (s=Q.begin() ; s != Q.end() ; ++s)
    //~ cout << *s << ",";
    //~ cout << "}" << endl;
    //~ // -- //


    // W = Sigma
    for (a=Sigma.begin() ; a != Sigma.end() ; ++a)
    {
        w.clear();
        w += *a;
        //w.push_back(*a);
        W.insert(w);
    }

    //	Si il n'y a qu'un seul Ètat initial, alors
    // son mot caractÈristique est epsilon
    if (QI.size()==1)
    {
        Qm=QI;
        if (verbose)
        {
            cout << "mot caractÈristique de " << *(Qm.begin()) << "\n";
            cout << "epsilon\n";
            cout << "longueur 0" << endl;
        }
    }

    // QI est la premiËre partie possible
    Memory.insert(QI);

    while(!W.empty() && !(Qm==Q))
    {
        w=*(W.begin());        // w=min W
        W.erase(W.begin());	 // W <- W \ {w}
        R=delta(QI,w);

        //~ // -- Affichage de delta (QI, w) -- //
        //~ cout << "delta(QI," << affiche(w) << ")={";
        //~ for (s = R.begin() ; s != R.end() ; ++s)
        //~ cout << *s << ",";
        //~ cout << "}" << endl;
        //~ // -- //

        // si R \notin Memory
        if (Memory.find(R) == Memory.end())
        {
            // on rajoute R ‡ Memory
            Memory.insert(R);
            // W = W \cup w\Sigma
            for (a=Sigma.begin() ; a != Sigma.end() ; ++a)
            {
                w += *a;
                //w.push_back(*a);
                W.insert(w);
                b=w.end();
                w.erase(--b);
            }

            // Si R contient un seul Ètat
            if (R.size() == 1)
            {
                // Si cet Ètat n'a pas dÈja eu un mot caractÈristique
                if (Qm.find(*(R.begin())) == Qm.end())
                {
                    // Alors on a trouver un nouvel Ètat dont on connait le
                    // mot caractÈristique (w)
                    Qm.insert(*(R.begin()));
                    if (verbose)
                    {
                        cout << "le mot caracteristique de l'Ètat " << *(R.begin()) << " est :\n";
                        cout << affiche(w) << "\n";
                        cout << "longueur " << w.size() << endl;
                    }
                }
            }
        }
    }
    return Qm==Q;
}

bool PFA::isPDFA()
{
    if (!isPFA())
    {
        return false;
    }
    erase_transitions(0,0); // erasing null transition
    if (iota.size() != 1)
        return false;
    Graph::const_iterator q;
    LettreSFunc::const_iterator a;
    SFunc::const_iterator r;
    for (q=phi.begin() ; q != phi.end() ; ++q)
    {
        for (a=q->second.begin() ; a!=q->second.end() ; ++a)
        {
            if (a->second.size() > 1)
                return false;
        }
    }
    return true;
}

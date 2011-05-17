/***************************************************************************
                                 pprfa.C
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
  Pour les notes, voir le fichier "pprfa.h".
  See "pprfa.h" file to view notes.
*/

#include "pprfa.H"

// Best Suppression
float PPRFA::best_transition_delete(Sample Stest, bool verbose) {

	if (phi.empty()) {
		return ERR(0);
	}
	else if (phi.begin()->second.empty()) {
		return ERR(1);
	}	
	else if (phi.begin()->second.begin()->second.empty()) {
		return ERR(2);
	}

    // Recherche des transitions à supprimer
    Graph::iterator q;
    LettreSFunc::iterator a;
    SFunc::iterator s;
    float min;
	PPRFA C=*this;

	float perplex;
	float minperplex;
	float epsprime;

	lisse();
    minperplex = perplexite(Stest);
    if (verbose) {
		cout << "\tepsprime = 0";
		cout << ", distance = " << minperplex << endl;
	}
    min=0;
	perplex=minperplex;
    do
    {
        epsprime = min;
		minperplex = perplex;
		
        min = phi.begin()->second.begin()->second.begin()->second;

        for (q = phi.begin(); q != phi.end(); ++q)
        {
            for (a = q->second.begin (); a != q->second.end (); ++a)
            {
                for (s = a->second.begin (); s != a->second.end (); ++s)
                {
                    if ((s->second < min) && (s->second > epsprime))
                    {
                        min = s->second;
                    }
                }
            }
        }
        if (min == epsprime) {
            break;
		}
        erase_transitions (min);
        rend_PFA ();
        lisse ();
        perplex = perplexite (Stest);
        *this = C;
        if (verbose) {
            cout << "\tepsprime = " << min;
			cout << ", distance = " << perplex << endl;
		}
    }
    while (perplex <= minperplex);
    if (verbose) {
        cout << "\t|| Meilleur epsprime trouve = " << epsprime;
		cout << ", mindist = " << minperplex << " ||" << endl;
	}
    *this = C;
    erase_transitions (epsprime);
    rend_PFA();

	return d_2(Stest);
}



RESULT PPRFA::eraseTransition(const Word &w, const Lettre l, const Word &v) {
	return PFA::eraseTransition(XR[w],l,XR[v]);
}

// fonction qui associe ‡ un etat le mot correspondant
Word PPRFA::word_of(const State x) const {
	WordSFunc::const_iterator r;
	for (r=XR.begin() ; r != XR.end() ; r++) {
		if (r->second == x)
			return r->first;
	}
	// on sort de la boucle sans avoir trouve !!!
	cerr << "WORD_OF !!! etat : " << x << endl;
	return "";
}

RESULT PPRFA::compute_XR(void) {
	StateSet W;
	StateSet R;
	State etat_cour;
	Word u;
	Graph::const_iterator q;
	LettreSFunc::const_iterator a;
	SFunc::const_iterator s;
	StateWordFunction WS;
    W.insert (iota.begin ()->first);
	
	XR.clear();
	XR[u]=iota.begin()->first;
	WS[iota.begin()->first]=u;
    while (!W.empty ())
    {
        etat_cour = *W.begin ();
        W.erase (W.begin ());
        q = phi.find (etat_cour);
        if (q != phi.end ())
        {
            for (a = q->second.begin (); a != q->second.end (); ++a)
            {
                for (s = a->second.begin ();  s != a->second.end (); ++s)
                {
					if (R.find(s->first) == R.end()) {
						u=WS[etat_cour] + string(1,a->first);
						XR[u]=s->first;
						WS[s->first]=u;
						W.insert (s->first);
						R.insert(s->first);
					}
                }
            }
        }
    }
	return VAL(0);
}




// addStates from q to generate w from P_q(w)
RESULT PPRFA::add_word( const State q, const Word &w) {
	if (w.empty()) {
		if ((tau.find(q) == tau.end()) || (tau.find(q)->second == 0)) {
			tau[q] = 0.01;
			return VAL(1);
		}
		else {
			return ERR(1);
		}
	}
	else if (delta(q, string(1,*w.begin())).empty()) {
		WordSFunc::const_iterator xr;
		for (xr = XR.begin() ; 
			(xr != XR.end()) && (xr->second != q) ; 
			xr++) 
		;
		
		if (xr == XR.end()) {
			return ERR(2);
		}
		
		Word u=xr->first;
		State r,s;
		Word::const_iterator a;
		r=q;
		for (a = w.begin() ; a != w.end() ; a++) {	
			u+=*a;
			s=addState(u,0,0);
			PFA::addTransition(r, *a , s, 0.01);
			r=s;
		}		
		return VAL(0);
	}
	else {
		return ERR(2);
	}
}

RESULT PPRFA::emmonde(void) {
	PFA::emmonde();
	WordSFunc::iterator r,s;
	for (r=XR.begin() ; r != XR.end() ; ) {
		if (Q.find(r->second) == Q.end()) {
			s=r;
			s++;
			XR.erase(r);
			r=s;
		}
		else {
			r++;
		}
	}
	return VAL(0);
}

// Add Prefix many Tree in order to recognize all the sample
RESULT PPRFA::put_prefix_tree(const Sample &S, const Word &v) {
	Sample::const_iterator w;
	StateSet Q1;
	StateSet::const_iterator q;
	Word::const_iterator a;
	for (w=S.begin() ;w != S.end() ; w++) {
		//~ // ----
		//~ cout << "Mot : " << w->first << endl;
		//~ // ----
	
		// we are in a PPRFA then there is exactly one initial state
		Q1.clear();
		Q1.insert(iota.begin()->first); 
		for (a=w->first.begin() ; a != w->first.end() ; a++) {
		
			//~ // ----
			//~ cout << "Q1={";
			//~ for (StateSet::const_iterator s=Q1.begin() ; s != Q1.end() ; s++) {
				//~ cout << *s << ",";
			//~ }
			//~ cout << "}" << endl;
			//~ // ----
		
			// si pour un q delta(q,a) est vide, on ajoute les etat
			// a partir de q pour reconnaitre la fin du mot
			for (q=Q1.begin() ;q != Q1.end() ;q++) {
				if (delta(*q , string(1,*a)).empty()) {
					//~ // ----
					//~ cout << "add_word(" << *q << ",'" << flush;
					//~ cout << string(a,w->first.end()) << "','" << flush;
					//~ cout << string(w->first.begin(), a) << "')"<< endl;
					//~ // ----
					add_word(*q,string(a,w->first.end()));
				}
			}
			Q1 = delta(Q1,string(1,*a));
		}

		//~ // ----
		//~ cout << "Q1={";
		//~ for (StateSet::const_iterator s=Q1.begin() ; s != Q1.end() ; s++) {
			//~ cout << *s << ",";
		//~ }
		//~ cout << "}" << endl;
		//~ // ----

		for (q = Q1.begin() ; q != Q1.end() ; q++) {
			if (tau.find(*q) == tau.end() || (tau.find(*q)->second == 0)) {
				//~ // ----
				//~ cout << "add_word(" << *q << ",'" << flush;
				//~ cout << "','" << flush;
				//~ cout << w->first << "')"<< endl;
				//~ // ----
				add_word(*q,"");
			}
		}
	}
	return VAL(0);
}

// Arbre prefixe
RESULT PPRFA::becomePrefixTree(Sample &S)
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
    
    w = S.begin();
    XR[w->first]=addNewState(1,0);

    // on commence apres epsilon d'ou le ++w au debut du for
    for (++w ; w != S.end() ; ++w)
    {
        v=w->first;
        XR[v]=addNewState(0,0);
        a=*(v.rbegin());
        v.erase(--v.end());
        PFA::addTransition(XR[v], a , XR[w->first],(float)S[w->first]/(float)S[v]);
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



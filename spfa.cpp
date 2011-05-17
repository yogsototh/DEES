/***************************************************************************
                                  spfa.cpp
                             -------------------
    begin                : Sun 7 Dec 2002
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
  Pour les notes, voir le fichier "spfa.h".
  See "ffa.h" file to view notes.
*/

#include "spfa.H"
#include <assert.h>

// Ajoute un Ètat
RESULT SPFA::addNewState (const float init, const float term)
{
    if (PFA_SAFE)
    {
        assert (init >= 0);
        assert (init <= 1);
        assert (term >= 0);
        assert (term <= 1);
    }
    return MA::addNewState (init, term);
    if (PFA_SAFE)
    {
        assert (isSPFA ());
    }
}

// Fonction qui ajoute une transition
RESULT SPFA::addTransition (const Transition & t, const float val)
{
    return MA::addTransition (t, val);
    if (PFA_SAFE)
        assert (isSPFA ());
}

// Ajout d'une transition
RESULT SPFA::addTransition (const State qdep, const Lettre a, const State qarr, const float val)
{
    return MA::addTransition (qdep, a, qarr, val);
    if (PFA_SAFE)
    {
        assert (isSPFA ());
    }
}

// Donne un SPFA aleatoire
RESULT SPFA::becomeRandom (const int nb_etats,
                           const int nb_lettres,
                           const int num_graphe,
                           const float densite,
                           const float prob_init,
                           const float prob_term,
                           const float min_trans, const float max_trans)
{

    // On gÈnËre un MA alÈatoire avec seulement des arÍtes positives
    MA::becomeRandom (nb_etats, nb_lettres,
                       num_graphe, densite,
                       prob_init, prob_term, min_trans, max_trans);

    // On rend le MA en SPFA
    return renormalise ();
}

// Cree un MA aleatoire avec un nombre maximal d'aretes
RESULT SPFA::becomeRandomMax (const int nb_etats, const int nb_lettres,
                              const int num_graphe, const int nb_succ,
                              const int nb_init,    const int nb_term,
                              const float min_trans,  const float max_trans)
{
    MA::becomeRandomMax(nb_etats, nb_lettres, num_graphe, nb_succ, nb_init, nb_term,min_trans,max_trans);
    return renormalise();
}


RESULT SPFA::becomeRandomPrefix (const int nb_etats,
                                 const int nb_lettres,
                                 const int num_graphe,
                                 const float densite,
                                 const float prob_init,
                                 const float prob_term,
                                 const float min_trans, const float max_trans)
{

    // On gÈnËre un MA alÈatoire avec seulement des arÍtes positives
    MA::becomeRandomPrefix (nb_etats, nb_lettres,
                             num_graphe, densite,
                             prob_init, prob_term, min_trans, max_trans);

    // on rend un SPFA
    return renormalise ();
}


RESULT SPFA::renormalise ()
{
    if (!Q.empty())
    {
        // --- On renormalise la fonction d'initialisation ---
        SFunc::iterator s;
        double sum, sum2;
        // on gere le cas ou il n'y a pas d'etat initial
        if (iota.empty ()) 
		{
            iota[*(Q.begin ())] = 1;
        }
        else
        {
            // on calcule 1/somme des iota
            sum = 0;
            for (s = iota.begin (); s != iota.end (); ++s)
            {
                sum += s->second;
            }
            
            // on renormalise tous les iota
            sum2 = 0;
            for (s = iota.begin (); s != iota.end (); s++)
            {
                s->second = s->second / sum;
                sum2 += s->second;
            }
            // on repare encore un peu
            iota.begin()->second += 1 - sum2;
        }

        // Pour chaque etat on regarde la somme de ce qui sort
        StateSet::const_iterator q;
        Graph::iterator g;
        LettreSFunc::iterator a;
        for (q = Q.begin (); q != Q.end (); q++)
        {	// pour tous les Ètats
            // on calcule ce qui sort de l'etat
            sum = val_outstate(*q);
            if (sum == 0)	// si rien ne sort de l'etat, on met ‡ jour tau
            {
                tau[*q] = 1;
            } 
			else
            {
                // on renormalise en prenant en compte les problËmes de virgule
				g=phi.find(*q);
                sum2 = 0;
                if (g != phi.end ())
                {
                    for (a = g->second.begin (); a != g->second.end (); ++a)
                    {
                        for (s = a->second.begin (); s != a->second.end (); ++s)
                        {
                            s->second /= sum;
                            sum2 += s->second;
                        }
                    }
                }
                // si la somme ne fait pas encore 1 on
                // met ‡ jour la valeur d'une arete ou de tau
                if (sum2 != 1)
                {
					// si tau est non nul c'est lui que l'on modifie
                    if ((tau.find (*q) != tau.end ()) || (g == phi.end()))
                    {
                        tau[*q] = 1 - sum2;
                    }
                    else
                    { // sinon on repare un peu une transition
                        bool trouve = false;
                        for (a = g->second.begin (); 
							(a != g->second.end ()) && !trouve; 
							++a)
                        {
                            s = a->second.begin ();
                            if (s != a->second.end ())
                            {
                                trouve = true;
                                s->second += 1 - sum2;
                            }
                        }
                    }
                }
            }
        }

        if (PFA_DEBUG)
        {
            if (!isSPFA ())
            {
                cerr << "ERREUR de la fonction SPFA::renormalise()" <<
                endl;
            }
        }
        return isSPFA ();
    }
    else
        return false;
}
// le logarithme de la fonction phi Ètendue aux mots
double
SPFA::philog (const State q, const Word & u, const State s,
              const Dictionnaire * dico) const
{
    Word::const_iterator a;
    PreciseSFunc V;	// Le vecteur contenant log(p(w))
    PreciseSFunc Vbis;
    StateSet::const_iterator r;
    PreciseSFunc::iterator si;
    double tmp;

    if (dico == NULL)
    {
        V[q] = 0;
        for (a = u.begin (); a != u.end (); a++)
        {


            Vbis.clear ();
            for (r = Q.begin (); r != Q.end (); r++)
            {
                // Calcul du nouveau vecteur
                tmp = 0;
                for (si = V.begin (); si != V.end (); si++)
                {
                    if (val_phi (si->first, *a, *r) != 0)
                    {
                        if (Vbis.find (*r) !=
                                Vbis.end ())
                        {
                            Vbis[*r] =
                                sumlog (Vbis[*r], si->second + log(val_phi(si->first,*a, *r)));
                        }
                        else
                        {
                            Vbis[*r] = si->second + log(val_phi(si->first,*a, *r));
                        }
                    }
                }
            }
            V = Vbis;
        }
        return V[s];
    }
    else			// les dictionnaires ne sont peut-etre pas Èquivalents
    {
        Word w;
        Dictionnaire::const_iterator tl;
        TrueLettre b;

        w.clear ();
        for (a = u.begin (); a != u.end (); ++a)
        {
            tl = dico->find (*a);
            if (tl == dico->end ())
            {
                cerr << "Erreur !!! SPFA::philog, bad dictionnary !" << endl;
            }
            else
            {
                b = tl->second;
                for (tl = alph.begin ();
                        (tl->second != b)
                        && (tl != alph.end ()); ++tl)
                    ;
                if (tl == alph.end ())
                {
                    cerr << "Erreur 2 !!! SPFA::philog, bad dictionnary" << endl;
                }
                else
                {
                    w += tl->first;
                }
            }
        }

        V[q] = 0;
        for (a = w.begin (); a != w.end (); ++a)
        {
            Vbis.clear ();
            for (r = Q.begin (); r != Q.end (); ++r)
            {
                // Calcul du nouveau vecteur
                tmp = 0;
                for (si = V.begin (); si != V.end (); ++si)
                {
                    if (val_phi (si->first, *a, *r) != 0)
                    {
                        if (Vbis.
                                find (*r) != Vbis.end ())
                        {
                            Vbis[*r] =sumlog (Vbis[*r],si->second+log(val_phi(si-> first, *a,*r)));
                        }
                        else
                        {
                            Vbis[*r] =si->second +log (val_phi(si-> first,*a,*r));
                        }
                    }
                }
            }
            V = Vbis;
        }
        return V[s];
    }

}

// la proba conditionnÈe i.e. p(q|u)
double
SPFA::plog (const State q, const Word & u, const Dictionnaire * dico) const
{
    Word::const_iterator a;
    PreciseSFunc V;	// Le vecteur contenant log(p(w))
    PreciseSFunc Vbis;
    SFunc::const_iterator e;
    StateSet::const_iterator r;
    PreciseSFunc::iterator s;
    double res, tmp;

    for (e = iota.begin (); e != iota.end (); ++e)
        V[e->first] = log (e->second);

    if ((dico == NULL) || (*dico == alph))
    {
        for (a = u.begin (); a != u.end (); ++a)
        {
            Vbis.clear ();
            for (r = Q.begin (); r != Q.end (); ++r)
            {
                // Calcul du nouveau vecteur
                tmp = 0;
                for (s = V.begin (); s != V.end (); ++s)
                {
                    if (val_phi (s->first, *a, *r) != 0)
                    {
                        if (Vbis.find (*r) !=
                                Vbis.end ())
                        {
                            Vbis[*r] = sumlog(Vbis[*r],
                                              s->second + log(val_phi(s->first,*a, *r)));
                        }
                        else
                        {
                            Vbis[*r] =
                                s->second + log(val_phi(s->first,*a, *r));
                        }
                    }
                }
            }
            V = Vbis;
        }

        //  res=exp(V[q] - sum_r V[r])
        if (V.find (q) != V.end ())
        {
            res = log ((double) 0);
            for (s = V.begin (); s != V.end (); ++s)
                res = sumlog (res, s->second);
            res = V[q] - res;
        }
        else
        {
            res = log ((double) 0);
        }

        return res;
    }
    else			// dico != NULL
    {
        Word w;
        Dictionnaire::const_iterator tl;
        TrueLettre b;

        w.clear ();
        for (a = u.begin (); a != u.end (); ++a)
        {
            tl = dico->find (*a);
            if (tl == dico->end ())
            {
                cerr << "Erreur !!! SPFA::philog, bad dictionnary !" << endl;
            }
            else
            {
                b = tl->second;
                tl = alph.begin() ;
                while	((tl != alph.end()) && (tl->second != b))
                    ++tl;
                if (tl == alph.end ())
                {
                    cerr << "Erreur 2 !!! SPFA::philog, bad dictionnary" << endl;
                }
                else
                {
                    w += tl->first;
                }
            }
        }


        for (a = w.begin (); a != w.end (); ++a)
        {
            Vbis.clear ();
            for (r = Q.begin (); r != Q.end (); ++r)
            {
                // Calcul du nouveau vecteur
                tmp = 0;
                for (s = V.begin (); s != V.end (); ++s)
                {
                    if (val_phi (s->first, *a, *r) != 0)
                    {
                        if (Vbis.find (*r) != Vbis.end ())
                        {
                            Vbis[*r] =
                                sumlog(Vbis[*r],
                                       s->second + log(val_phi(s->first,*a, *r)));
                        }
                        else
                        {
                            Vbis[*r] = s->second + log(val_phi(s->first,*a, *r));
                        }
                    }
                }
            }
            V = Vbis;
        }

        //  res=exp(V[q] - sum_r V[r])
        if (V.find (q) != V.end ())
        {
            res = log ((double) 0);
            for (s = V.begin (); s != V.end (); ++s)
                res = sumlog (res, s->second);
            res = V[q] - res;
        }
        else
        {
            res = log ((double) 0);
        }

        return res;
    }
}

// la fonction p directe
double SPFA::p_directe (const Word & u, const Dictionnaire * dico) const
{
	PreciseSFunc F;
	PreciseSFunc::const_iterator x;
	SFunc::const_iterator t;
	double res;
	if ((dico == NULL) || (*dico == alph)) {
		forward(F,u);
		res = 0.0;
		for (x=F.begin() ; x!=F.end() ; ++x) {
			t = tau.find(x->first);
			if (t != tau.end()) {
				res += x->second * double(t->second);
			}
		}
		return res;
	}
	else {
		// le dictionnaire n'est pas identique, on doit traduire le mot
		Word v;
		Word::const_iterator a;
        Dictionnaire::const_iterator tl;
        Dictionnaire::const_iterator tl2;
        map<Lettre, Lettre> traducteur;
        for (tl = dico->begin() ; tl != dico->end() ; ++tl)
        {
            tl2 = alph.begin();
            while ((tl2 != alph.end()) && (tl2->second != tl->second))
                ++tl2;
            if (tl2 == alph.end())
            {
                cerr << "Erreur !!! SPFA::philog, bad dictionnary ! lettre " << tl->second << endl;
            }
            else
            {
                traducteur[tl->first] = tl2->first;
            }
        }
		for (a=u.begin() ; a != u.end() ; ++a) {
			v += traducteur[*a];
		}
		
		forward(F,v);
		res = 0.0;
		for (x=F.begin() ; x!=F.end() ; ++x) {
			t = tau.find(x->first);
			if (t != tau.end()) {
				cout << "^^" << x->second << " * " << t->second << endl;
				res += x->second * double(t->second);
			}
		}
		return res;
	}
}


// la fonction p directe
double SPFA::p_bar_directe (const Word & u, const Dictionnaire * dico) const
{
	PreciseSFunc F;
	PreciseSFunc::const_iterator x;
	SFunc::const_iterator t;
	double res;
	if ((dico == NULL) || (*dico == alph)) {
		forward(F,u);
		res = 0.0;
		for (x=F.begin() ; x!=F.end() ; ++x) {
				res += x->second;
		}
		return res;
	}
	else {
		// le dictionnaire n'est pas identique, on doit traduire le mot
		Word v;
		Word::const_iterator a;
        Dictionnaire::const_iterator tl;
        Dictionnaire::const_iterator tl2;
        map<Lettre, Lettre> traducteur;
        for (tl = dico->begin() ; tl != dico->end() ; ++tl)
        {
            tl2 = alph.begin();
            while ((tl2 != alph.end()) && (tl2->second != tl->second))
                ++tl2;
            if (tl2 == alph.end())
            {
                cerr << "Erreur !!! SPFA::philog, bad dictionnary ! lettre " << tl->second << endl;
            }
            else
            {
                traducteur[tl->first] = tl2->first;
            }
        }
		for (a=u.begin() ; a != u.end() ; ++a) {
			v += traducteur[*a];
		}
		
		forward(F,v);
		res = 0.0;
		for (x=F.begin() ; x!=F.end() ; ++x) {
			res += x->second;
		}
		return res;
	}
}


// le logarithme de la fonction p
double SPFA::plog (const Word & u, const Dictionnaire * dico) const
{
	PreciseSFunc F;
	PreciseSFunc::const_iterator x;
	SFunc::const_iterator t;
	double res;
	double tmp;
	if ((dico == NULL) || (*dico == alph)) {
		logforward(F,u);
		res = log(0.0);
		for (x=F.begin() ; x!=F.end() ; ++x) {
			t = tau.find(x->first);
			if (t != tau.end()) {
				tmp=x->second+log(double(t->second));
				if (!isinf(tmp)) {				
					res = sumlog(res, tmp);
				}
			}
		}
		return res;
	}
	else {
		// le dictionnaire n'est pas identique, on doit traduire le mot
		Word v;
		Word::const_iterator a;
        Dictionnaire::const_iterator tl;
        Dictionnaire::const_iterator tl2;
        map<Lettre, Lettre> traducteur;
        for (tl = dico->begin() ; tl != dico->end() ; ++tl)
        {
            tl2 = alph.begin();
            while ((tl2 != alph.end()) && (tl2->second != tl->second))
                ++tl2;
            if (tl2 == alph.end())
            {
                cerr << "Erreur !!! SPFA::philog, bad dictionnary ! lettre " << tl->second << endl;
            }
            else
            {
                traducteur[tl->first] = tl2->first;
            }
        }
		for (a=u.begin() ; a != u.end() ; ++a) {
			v += traducteur[*a];
		}
		
		logforward(F,v);
		res = log(0.0);
		for (x=F.begin() ; x!=F.end() ; ++x) {
			t = tau.find(x->first);
			if (t != tau.end()) {
				res = sumlog(res,x->second + log(double(t->second)));
			}
		}
		return res;
	}
}

// ---------------- calcul des probas

double SPFA::p_bar(const Word &u, const Dictionnaire * dico) const  {
	return exp(plog_bar(u));
}

double SPFA::plog_bar(const Word &u, const Dictionnaire * dico) const {
	// PreciseSFunc F;
	// PreciseSFunc::const_iterator f;
	// double res;
	// logforward(F,u);
	// res = log(0.0);
	// for (f=F.begin() ; f != F.end() ; f++) {
	// 	res = sumlog(res, f->second);
	// }
	// return res;
	// -------------
	
	
	PreciseSFunc F;
	PreciseSFunc::const_iterator x;
	SFunc::const_iterator t;
	double res;
	if ((dico == NULL) || (*dico == alph)) {
		logforward(F,u);
		res = log(0.0);
		for (x=F.begin() ; x!=F.end() ; ++x) {
			if (!isinf(x->second)) {
				res = sumlog(res,x->second);
			}
		}
		return res;
	}
	else {
		// le dictionnaire n'est pas identique, on doit traduire le mot
		Word v;
		Word::const_iterator a;
        Dictionnaire::const_iterator tl;
        Dictionnaire::const_iterator tl2;
        map<Lettre, Lettre> traducteur;
        for (tl = dico->begin() ; tl != dico->end() ; ++tl)
        {
            tl2 = alph.begin();
            while ((tl2 != alph.end()) && (tl2->second != tl->second))
                ++tl2;
            if (tl2 == alph.end())
            {
                cerr << "Erreur !!! SPFA::philog, bad dictionnary ! lettre " << tl->second << endl;
            }
            else
            {
                traducteur[tl->first] = tl2->first;
            }
        }
		for (a=u.begin() ; a != u.end() ; ++a) {
			v += traducteur[*a];
		}
		
		logforward(F,v);
		res = log(0.0);
		for (x=F.begin() ; x!=F.end() ; ++x) {
			if (!isinf(x->second)) {
				res = sumlog(res,x->second);
			}
		}
		return res;
	}
	
}

RESULT SPFA::forward(PreciseSFunc &F, const Word &u) const {
	PreciseSFunc init;
	SFunc::const_iterator q;
	for (q=iota.begin() ; q!= iota.end() ; q++)
	{
		init[q->first] = q->second;
	}
	return forwardprecalc(F,u,init);
}

RESULT SPFA::logforward(
    PreciseSFunc &F,
    const Word &u) const
{
    PreciseSFunc init;
    SFunc::const_iterator q;
    for (q=iota.begin() ; q!=iota.end() ; ++q)
    {
        init[q->first] = q->second;
    }
    return logforwardprecalc(F,u,init);
}


RESULT SPFA::forwardprecalc(
							   PreciseSFunc &F,
							   const Word &u,
							   const PreciseSFunc &init_vector) const
{
    F.clear();
	
    PreciseSFunc V;
    PreciseSFunc::const_iterator e;
    Word::const_iterator a;
    PreciseSFunc Vbis;
    PreciseSFunc::iterator s;
	Graph::const_iterator q;
	LettreSFunc::const_iterator b;
	SFunc::const_iterator r;
	
	
	for (e = init_vector.begin (); e != init_vector.end (); ++e) {
		if (e->second != 0) {
			V[e->first] = e->second;
		}
	}
	
    for (a = u.begin (); a != u.end (); ++a)
    {
        Vbis.clear ();
		for (s=V.begin() ; s != V.end() ; ++s) {
			q=phi.find(s->first);
			if (q != phi.end() ) {
				b=q->second.find(*a);
				if (b != q->second.end()) {
					for (r=b->second.begin() ; r != b->second.end() ; ++r) {
						if (Vbis.find(r->first) != Vbis.end()) {
							Vbis[r->first] = Vbis[r->first] + r->second*s->second;
						}
						else {
							Vbis[r->first] = r->second*s->second;
						}
					}
				}
			}
		}
		V = Vbis;
    }
	F=V;
    return VAL(0);
}

RESULT SPFA::logforwardprecalc(
    PreciseSFunc &F,
    const Word &u,
    const PreciseSFunc &init_vector) const
{
    F.clear();

    PreciseSFunc V;
    PreciseSFunc::const_iterator e;

    for (e = init_vector.begin (); e != init_vector.end (); ++e) {
		if (e->second != 0) {
			V[e->first] = log (e->second);
		}
	}

    Word::const_iterator a;
    PreciseSFunc Vbis;
    PreciseSFunc::iterator s;

	Graph::const_iterator q;
	LettreSFunc::const_iterator b;
	SFunc::const_iterator r;

    for (a = u.begin (); a != u.end (); ++a)
    {
        Vbis.clear ();
		for (s=V.begin() ; s != V.end() ; ++s) {
			q=phi.find(s->first);
			if (q != phi.end() ) {
				b=q->second.find(*a);
				if (b != q->second.end()) {
					for (r=b->second.begin() ; r != b->second.end() ; ++r) {
						if (Vbis.find(r->first) != Vbis.end()) {
							Vbis[r->first] = sumlog(Vbis[r->first],
												log(r->second) + s->second);
						}
						else {
							Vbis[r->first] = log(r->second) + s->second;
						}
					}
				}
			}
		}
		V = Vbis;
	
		//~ for (e=Vbis.begin() ; e!=Vbis.end() ; ++e) {
			//~ if (!isinf(e->second)) {
				//~ V[e->first]=e->second;
			//~ }
		//~ }
    }
	F=V;
    return VAL(0);
}


// ---------------- LA PARTIE BAULM WELCH

RESULT SPFA::logforward(
    list < PreciseSFunc > &F,
    const Word &u) const
{
    PreciseSFunc init;
    SFunc::const_iterator q;
    for (q=iota.begin() ; q!=iota.end() ; ++q)
    {
        init[q->first] = q->second;
    }
    return logforwardprecalc(F,u,init);
}

RESULT SPFA::logforwardprecalc(
    list < PreciseSFunc > &F,
    const Word &u,
    const PreciseSFunc &init_vector) const
{
    F.clear();

    PreciseSFunc V;
    PreciseSFunc::const_iterator e;
    for (e = init_vector.begin (); e != init_vector.end (); ++e) {
		if (e->second != 0) {
			V[e->first] = log (e->second);
		}
	}
    F.push_back(V);

    Word::const_iterator a;
    PreciseSFunc Vbis;
    PreciseSFunc::iterator s;

	Graph::const_iterator q;
	LettreSFunc::const_iterator b;
	SFunc::const_iterator r;

    for (a = u.begin (); a != u.end (); ++a)
    {
        Vbis.clear ();
		for (s=V.begin() ; s != V.end() ; ++s) {
			q=phi.find(s->first);
			if (q != phi.end() ) {
				b=q->second.find(*a);
				if (b != q->second.end()) {
					for (r=b->second.begin() ; r != b->second.end() ; ++r) {
						if (Vbis.find(r->first) != Vbis.end()) {
							Vbis[r->first] = sumlog(Vbis[r->first],
												log(r->second) + s->second);
						}
						else {
							Vbis[r->first] = log(r->second) + s->second;
						}
					}
				}
			}
		}
		V = Vbis;
		//~ for (e=Vbis.begin() ; e!=Vbis.end() ; ++e) {
			//~ if (!isinf(e->second)) {
				//~ V[e->first]=e->second;
			//~ }
		//~ }
        F.push_back(V);
    }
    return VAL(0);
}

RESULT SPFA::logbackward(
    list < PreciseSFunc > &B,
    const Word &u) const
{
    // copie de tau dans term
    PreciseSFunc term;
    SFunc::const_iterator q;
    for (q=tau.begin() ; q != tau.end() ; ++q)
    {
        term[q->first]=q->second;
    }
    return logbackwardprecalc(B,u,term);
}

RESULT SPFA::logbackwardprecalc(
    list < PreciseSFunc > &B,
    const Word &u,
    const PreciseSFunc &term) const
{

	B.clear();

    PreciseSFunc::const_iterator e;
    PreciseSFunc V;
    for (e = term.begin (); e != term.end (); ++e) {
		if (e->second != 0) {
			V[e->first] = log (e->second);
		}
	}

    B.push_front(V);

	Word::const_reverse_iterator a;
	Graph::const_iterator q;
	LettreSFunc::const_iterator b;
	SFunc::const_iterator s;
    PreciseSFunc Vbis;


    for (a = u.rbegin (); a != u.rend (); ++a)
    {
        Vbis.clear ();
		for (q=phi.begin() ; q != phi.end() ; ++q) {
			b=q->second.find(*a);
			if (b != q->second.end()) {
				for (s=b->second.begin() ; s != b->second.end() ; ++s) {
					if (V.find(s->first) != V.end()) {
						if (Vbis.find(q->first) == Vbis.end()) {
							Vbis[q->first] = log(s->second) + V[s->first];
						}
						else {
							Vbis[q->first] = sumlog( Vbis[q->first],
											log(s->second) + V[s->first]);
						}
					}
				}
			}
		}
	
		V=Vbis;
		B.push_front(V);
    }
    return VAL(0);
}

// Compte le nombre de fois que l'on passe par chaque transition de
// facon pondÈree.
RESULT SPFA::TransitionCount(const Sample &S,
                              TransitionFunction &T,
							  SFunc &Iota,
							  SFunc &Tau) const
{
	TransitionFunction::iterator t;
	SFunc::iterator s;
	Sample::const_iterator u;
	// on initialise le compte d'utilisation des transitions
    for (t=T.begin() ; t != T.end() ; ++t)
    {
		t->second = 0;
	}
	for (s=Iota.begin() ; s != Iota.end() ; ++s) {
		s->second = 0;
	}
	for (s=Tau.begin() ; s != Tau.end() ; ++s) {
		s->second=0;
	}
	
    for (u=S.begin() ; u != S.end() ; ++u)
    {
		TransitionCount(u->first, T, Iota, Tau, u->second);
		//~ // ----
        //~ cout << "\t: " << affiche(u->first) << "," << u->second << " : ";
        //~ for (t=T.begin() ; t != T.end() ; ++t) {
			//~ cout << "(";
			//~ cout << t->first.qdep << ",";
			//~ cout << t->first.a << ",";
			//~ cout << t->first.qarr << ")=";
			//~ cout << t->second << ",";
        //~ }			
		//~ for (s=Iota.begin() ; s != Iota.end() ; ++s) {
			//~ cout << " i(" << s->first << ")=" << s->second;
		//~ }
		//~ for (s=Tau.begin() ; s!=Tau.end() ; ++s) {
			//~ cout << " t(" << s->first << ")=" << s->second;
		//~ }
        //~ cout << endl;
		//~ // ----
    }
	return VAL(0);
}


// Compte le nombre de fois que l'on passe par chaque transition de
// facon pondÈree.
RESULT SPFA::TransitionCount(const Word &u,
                              TransitionFunction &T,
							  SFunc &Iota,
							  SFunc &Tau,
                              const int nb_occurences) const
{

    // Algorithme renvoie sum_{t} \eta_t(i, a)
    // \eta_t(i,a) = \alph_t(i) x \phi(i, w_t=a ,j) x \beta_t(j) / P(w)
    // on calcule alpha et beta pour tout i

    list < PreciseSFunc > F; // forward
    list < PreciseSFunc > B; // backward
    list < PreciseSFunc >::const_iterator f;
	list < PreciseSFunc >::const_iterator b;
    logforward(F,u);
    logbackward(B,u);
	

	Word::const_iterator a;	
	PreciseSFunc::const_iterator x;
	PreciseSFunc::const_iterator y;
	Transition t;
	double plogu=log(0.);// =plog(u);
	double tmp;

	//~ // AFFICHAGE DE F et B
	//~ cout << "Forward" << endl;
	//~ for (f=F.begin() ; f != F.end() ; ++f) {
		//~ for (x = f->begin() ; x != f->end() ; ++x) {
			//~ cout << "(" << x->first << "," << x->second << ") ";
		//~ }
		//~ cout << endl;
	//~ }

	//~ cout << "Backward" << endl;
	//~ for (b=B.begin() ; b != B.end() ; ++b) {
		//~ for (y = b->begin() ; y!=b->end() ; ++y) {
			//~ cout << "(" << y->first << "," << y->second << ") ";
		//~ }
		//~ cout << endl;
	//~ }

	f=F.begin();
	b=B.begin();
	for (x=f->begin() ; x != f->end() ; ++x) {
		y=b->find(x->first);
		plogu = sumlog(plogu, x->second + y->second);
	}

	//~ cout << "plogu = " << plogu << endl;
	//~ cout << "plogu calcule = " << plog(u) << endl;


	f=F.begin();
	b=B.begin();
	for (x=f->begin() ; x != f->end() ; ++x) {
		y=b->find(x->first);		
		if ((Iota.find(x->first) != Iota.end()) && (y != b->end())) {
			tmp = exp(x->second + y->second + log(float(nb_occurences)) - plogu);
			if (! isnan(tmp)) {
				Iota[x->first] += tmp;
			}
			else {
				cout << "IOTA nan : " << x->first;
				cout << "," << y->first;
				cout << ":" << x->second;
				cout << " " << y->second;
				cout << " " << log(float(nb_occurences));
				cout << " " << plogu;
				cout << endl;
			}
		}		
	}

	f=--F.end();
	b=--B.end();
	for (x=f->begin() ; x != f->end() ; ++x) {
		y = b->find(x->first);
		if ((Tau.find(x->first) != Tau.end()) && (y != b->end())) {
			tmp = exp(x->second + y->second + log(float(nb_occurences)) - plogu);
			if (!isnan(tmp)) {
				Tau[x->first] += tmp;
			}
			else {
				cout << "TAU nan : " << x->first;
				cout << "," << y->first;
				cout << ":" << x->second;
				cout << " " << y->second;
				cout << " " << log(float(nb_occurences));
				cout << " " << plogu;
				cout << endl;
			}
		}
	}

    for (a = u.begin(), 
		f=F.begin() , 
		b=B.begin() ,
		b++ ; // on avance de ce qu'il faut
	
		a != u.end() ; 
	
		++a, ++f, ++b) 
	{
		for (x=f->begin() ; x != f->end() ; ++x) 
		{
			for (y=b->begin() ; y != b->end() ; ++y) 
			{
				t.qdep = x->first;
				t.a = *a;
				t.qarr = y->first;
				if (T.find(t) != T.end())
				{
					if (val_phi(x->first, *a, y->first) != 0) {
						tmp= exp(x->second + 
									y->second + 
									log(val_phi(x->first, *a, y->first)) -
									plogu + log(float(nb_occurences)));	
					
						if (!isnan(tmp)) {
							T[t] +=tmp;
						}
						else {
							cout << "PHI nan : " << x->first;
							cout << "," << *a;
							cout << "," << y->first;
							cout << ":" << x->second;
							cout << " " << y->second;
							cout << " " << log(val_phi(x->first, *a, y->first));
							cout << " " << plogu;
							cout << " " << log(float(nb_occurences));
							cout << endl;
						}
					}
				}
			}
		}		
	}
	return VAL(0);
}

RESULT SPFA::BaumWelch(const Sample &S,
                        TransitionFunction &T,
						SFunc &Iota,
						SFunc &Tau,
                        int nb_tours,
                        bool verbose)
{

	if (S.isprefixiel()) {
		cerr << "BaumWelch():: Sample prefixiel !" << endl;
		return ERR(1);
	}

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// Il faut Iota = tous les Ètats initiaux ou rien
	// si (q,a,r) est dans T ou q dans Tau, il faut tous les (q,x,s) + tau(q).
	// Sinon la renormalisation devient fausse
    Sample::const_iterator u;
    TransitionFunction::iterator t;

	// on initialise la probabilite de sortir de chaque etats
	SFunc count_out;
	Graph::const_iterator q;
	LettreSFunc::const_iterator b;
	SFunc::iterator s;
//	int i=0;
	double sum;

	while (nb_tours-- != 0)
    {	
        //~ if (verbose)
        //~ {
            //~ cout << "[BM " << nb_tours << "] ";
            //~ cout.flush();
        //~ }	
	
		TransitionCount(S,T,Iota,Tau);
	
		for (t=T.begin() ; t != T.end() ; ++t) {
			count_out[t->first.qdep] += t->second;
		}
		for (s=Tau.begin() ; s!=Tau.end() ; ++s) {
			count_out[s->first] += s->second;
		}	
		sum=0;
		for (s=Iota.begin() ; s!= Iota.end() ; ++s) {
			sum +=s->second;
		}
	
        for (t=T.begin() ; t != T.end() ; ++t)
        {
			if (t->second != 0) {
				SPFA::addTransition(t->first, t->second / count_out[t->first.qdep]);
			}
			else {
				SPFA::addTransition(t->first, 0);
			}
        }
		for (s=Tau.begin() ; s != Tau.end() ; ++s) {
			if (s->second != 0) {
				tau[s->first]= s->second / count_out[s->first];
			}
			else {
				tau[s->first]=0;
			}
		}
		for (s=Iota.begin() ; s != Iota.end() ; ++s) {
			iota[s->first] = s->second / sum;
		}
	
		erase_transitions();
	
        //~ if (verbose)
        //~ {
			//~ cout << endl;
            //~ for (t=T.begin() ; t != T.end() ; ++t)
            //~ {
				//~ cout << "(";
                //~ cout << t->first.qdep;
				//~ cout << ",";
				//~ cout << t->first.a;
				//~ cout << ",";
				//~ cout << t->first.qarr;
				//~ cout << ")=";
				//~ cout << t->second / count_out[t->first.qdep];
				//~ if (++i%5 == 0) cout << "\n";
				//~ else cout << " ";
            //~ }
			//~ for (s = Tau.begin() ; s != Tau.end() ; ++s) {
				//~ cout << "t(" << s->first << ")=" << s->second / count_out[s->first];
				//~ if (++i%5 == 0) cout << "\n";
				//~ else cout << " ";
			//~ }
			//~ for (s = Iota.begin() ; s != Iota.end() ; ++s) {
				//~ cout << "i(" << s->first << ")=" << s->second / sum;
				//~ if (++i%5 == 0) cout << "\n";
				//~ else cout << " ";
			//~ }
			//~ cout << endl;
        //~ }
		
		count_out.clear();
    }

    return VAL(0);
}

// fonction qui calcule les probas qui sortes d'un Ètat
// transitions + arret
float SPFA::val_outstate(const State q) const {
	float sum;
	Graph::const_iterator r;
	LettreSFunc::const_iterator a;
	SFunc::const_iterator s;

	sum=val_tau(q);
	r=phi.find(q);
	if (r != phi.end() ) {
		for (a=r->second.begin() ; a != r->second.end() ; ++a) {
			for (s=a->second.begin() ; s != a->second.end() ; ++s) {
				sum += s->second;
			}
		}
	}
	return sum;
}

RESULT SPFA::erase_bad_states(void) {

	StateSet::iterator q;
	LettreSFunc::iterator a;
	SFunc::iterator s;
	for (q=Q.begin() ; q != Q.end() ; ++q) {
		if (val_outstate(*q) == 0) {
			erase(*q);
		}
	}
	return VAL(0);
}

/***************************************************************************
                       spfa.h  -  Semi-Probabilistic Finite Automaton
                             -------------------
    begin                : 7 Dec 2002
    copyright            : (C) 2002 by Yann Esposito
    email                : esposito@cmi.univ-mrs.fr
    ***************************************************************************
 
    ***************************************************************************
    *                                                                         *
    *   This program is free software; you can redistribute it and/or modify  *
    *   it under the terms of the GNU General Public License as published by  *
    *   the Free Software Foundation; either version 2 of the License, or     *
    *   (at your option) any later version.                                   *
    *                                                                         *
    ***************************************************************************
 
  ____________ HISTORIQUE et DÈscription plus prÈcise _____________________
 
Classe qui hÈrite de la classe ffa, la seule diffÈrence est que les
valeurs de iota, phi et tau doivent Ítre positives et
sum_{q} iota[q]=1.
 
pour tout etat q : tau[q] + sum_{a,r} phi[q][a][r] = 1
  _________________________________________________________________________
*/

#ifndef SPFA_H
#define SPFA_H

#include "ma.H"
#include "sample.H"
#include <list>


typedef map<Word, State, ordre_mot> WordSFunc;
typedef map<State,State> StateSFunc;
typedef map<State,Word> StateWordFunction;
typedef set<Word,ordre_mot> WordSet;


class SPFA : public MA
{
public :

    // Ajoute un Ètat
    //RESULT addState(const float init=0, const float term=1);
    RESULT addNewState (const float init, const float term);
    // Fonction qui ajoute une transition
    RESULT addTransition(const Transition &t, const float val);
    // Ajout d'une transition
    RESULT addTransition(const State qdep, const Lettre a, const State qarr, const float val);

    RESULT becomeRandom(const int num_etat,        // le nombre d'Ètats
                        const int num_lettre,      // le nombre de lettres
                        const int num_graphe = 0,  // le numero du graphe
                        const float densite=0.3,  // densitÈ du ffa (0 ffa vide, 1 ffa complet)
                        const float prob_init=0.3, // la probabilitÈ pour un Ètat d'etre initial
                        const float prob_term=0.3, // probabilitÈ pour un Ètat d'Ítre terminal
                        const float min_trans = 0, // la valeur minimale des transitions
                        const float max_trans = 1); // la valeur maximale des transitions

    RESULT becomeRandomPrefix(const int nb_etats,        // le nombre d'Ètats
                              const int nb_lettres,      // le nombre de lettres
                              const int num_graphe = 0,  // le numero du graphe
                              const float densite=0.3, // densitÈ du ffa (0 ffa vide, 1 ffa complet)
                              const float prob_init=0.3, // la probabilitÈ pour un Ètat d'etre initial
                              const float prob_term=0.3, // probabilitÈ pour un Ètat d'Ítre terminal
                              const float min_trans = 0, // la valeur minimale des transitions
                              const float max_trans = 1); // la valeur maximale des transitions

    // Cree un MA aleatoire avec un nombre maximal d'aretes
    RESULT becomeRandomMax (const int nb_etats, const int nb_lettres,
                            const int num_graphe, const int nb_succ,
                            const int nb_init,    const int nb_term,
                            const float min_trans,  const float max_trans);

    // --- Output methods ---
protected:
    // log (exp(p1) + exp(p2)) = p1 + log(1 + exp(p2-p1)) pour p1 > p2
    // good approximation if p1=log(x1) and p2=log(x2) sumlog(p1,p2) return a good approximation
    // of log(x1 + x2)
    inline double sumlog(const double p1, const double p2) const
    {
		return (p1>p2)?(p1 + log(1+exp(p2-p1))):(p2 + log(1+exp(p1-p2)));
    }
public:
    // extension of the phi function to words
    inline double phiext(const State q, const Word &u, const State s, const Dictionnaire *dico=NULL) const
    {
        return exp(philog(q,u,s,dico));
    }
    // logarith of the extended phi
    double philog(const State q, const Word &u, const State s, const Dictionnaire *dico=NULL) const;
    // la fonction p renvoie la valeur d'un mot (probabilitÈ dans les PFAs)
    inline double p(const Word &u, const Dictionnaire *dico=NULL) const
    {
        return exp(plog(u,dico));
    }
    // le logarithme de la fonction p
    double plog(const Word &u, const Dictionnaire *dico=NULL) const;

	// la fonction p calculée de façon directe !!! ATTENTION AUX ARRONDIS !!
    double p_directe(const Word &u, const Dictionnaire *dico=NULL) const;

	
    // probability to be in state q having read the word u
    inline double p(const State q, const Word &u, const Dictionnaire *dico=NULL)
    {
        return exp(plog(q,u,dico));
    }
    // logarithm of the probability to be in state q having read the word u
    double plog(const State q, const Word &u, const Dictionnaire *dico=NULL) const;

	// return probability to begin by u
	double p_bar(const Word &u, const Dictionnaire * dico=NULL) const;
	// return probability to begin by u without log trick
	double p_bar_directe(const Word &u, const Dictionnaire * dico=NULL) const;
	// return the log probability to begin by u
	double plog_bar(const Word &u, const Dictionnaire * dico=NULL) const;

	// return the forward vector
    RESULT logforward(PreciseSFunc &F, const Word &u) const;
    RESULT logforwardprecalc(PreciseSFunc &F,
                             const Word &u, const PreciseSFunc &init) const;

	RESULT forward(PreciseSFunc &F, const Word &u) const;
    RESULT forwardprecalc(PreciseSFunc &F,
                             const Word &u, const PreciseSFunc &init) const;
	
	
    // return the forward list vectors
    RESULT logforward(list < PreciseSFunc > &F, const Word &u) const;
    RESULT logforwardprecalc(list < PreciseSFunc > &F,
                             const Word &u, const PreciseSFunc &init) const;
    RESULT logbackward(list < PreciseSFunc > &B, const Word &u) const;
    RESULT logbackwardprecalc(list < PreciseSFunc > &B, // la liste des vecteurs
                              const Word &u, // le mot
							  const PreciseSFunc & term) const; // le vecteur terminal
							  

	// Apprentissage
	RESULT BaumWelch(const Sample &S, 
					TransitionFunction &T, 
					SFunc &Iota,
					SFunc &Tau,
					int nb_tours, 
					bool verbose=false);
	RESULT TransitionCount(	const Sample &S,
                            TransitionFunction &T,
							SFunc &Iota,
							SFunc &Tau) const;

    RESULT TransitionCount(	const Word &u, 
							TransitionFunction &T, 
							SFunc &Iota, 
							SFunc &Tau, 
							const int nb_apparitions) const;
							
	float count_nb_pass(const State x, const Sample &S) const;

public:
    RESULT renormalise(void);
	RESULT erase_bad_states(void);
private:
	float val_outstate(const State q) const;
};
#endif

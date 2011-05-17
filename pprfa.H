/***************************************************************************
                          pfa.h  -  Probabilistic Finite Automaton
                             -------------------
    begin                : 7 Dec 2002
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

/*________________________  DÈscription plus precise ______________________*\
  Classe qui herite de la classe spfa, la seule diffÈrence est que l'on
  peut generer des mots avec une distribution de probabilite.
  \***************************************************************************/

#ifndef PPRFA_H
#define PPRFA_H

#include "pfa.H"
#include "simplex.H"


typedef enum {determinist, positive, nonconstrained} T_ModeVariables;
typedef enum {begin, end} T_ModeReturn;
typedef enum {epsfixed, variable, word_variable} T_ModeEpsilon;


class PPRFA : public PFA
{
protected:
    WordSFunc XR;

public:
    // become the PTA representing S but successor of v are returned
    RESULT becomeQuasiPrefixTree (	const Sample &S,   const Word &v, const SFunc &solution);
	// erase the transition w,l,v
	RESULT eraseTransition(const Word &w, const Lettre l, const Word &v);
	// emmonde
	RESULT emmonde(void);

protected:
    // fonction necessaire ‡ sol
    RESULT liste_mots_associes(WordSet &W,
                               const Word &v,
                               const Sample &S,
                               int maxmots = 10) const;
    // fonction necessaire ‡ liste_mots_associes
    RESULT ajoute_mots_associes(WordSet &W,
                                const Word &v,
                                const Sample &S,
                                int maxmots = INT_MAX) const;
    // renvoie vide si le systËme lineaire I(v,R,S,precision) n'a pas de solution et
    // sinon renvoie une solution de I i.e. un  ensemble de variables X_q telles que
    // v^{-1}P = \sum_{q \in R} X_qP_q --  en rappelant que R[w]=q => P_q = w^{-1}P
    // on ne prend en compte que au plus max mots.
    RESULT solmax(T_ModeVariables modeVariables,
				  SFunc &solution,  // La solution
                  const Word &v,
                  const Sample &S,
                  Simplex &simplx,
                  const double precision,
                  const T_ModeReturn moderet,
                  const T_ModeEpsilon modeeps,
                  set
                      <Word, ordre_mot> &W,
                      unsigned int max) const;
    // renvoie vide si le systËme lineaire I(v,R,S,precision) n'a pas de solution et
    // sinon renvoie une solution de I i.e. un  ensemble de variables X_q telles que
    // v^{-1}P = \sum_{q \in R} X_qP_q --  en rappelant que R[w]=q => P_q = w^{-1}P
    RESULT sol(T_ModeVariables modeVariables,
			   SFunc &solution,  // La solution
               const Word &v,
               const Sample &S,
               Simplex &simplx,
               const double precision,
               const T_ModeReturn moderet,
               const T_ModeEpsilon modeeps,
               set
                   <Word, ordre_mot> &W,
                   int maxmots=10,
                   bool Wcalculated = false) const;

    inline RESULT addTransition(const Word &wdep, const Lettre a, const Word &warr, const float val)
    {
        return PFA::addTransition(XR[wdep], a, XR[warr], val);
    }

    inline RESULT addTransition(const Word &wdep, const Lettre a, const State qarr, const float val)
    {
        return PFA::addTransition(XR[wdep], a, qarr, val);
    }

    inline void vide()
    {
        XR.clear();
        PFA::vide();
    }


    inline RESULT addState(const Word &w, const float init=0, const float term=0)
    {
        WordSFunc::const_iterator q;
        q=XR.find(w);
        if (q != XR.end()) {
            return MA::addState( q->second, init, term);
		}
        else
        {
            return XR[w]=addNewState(init, term);
        }
    }

	float best_transition_delete(Sample Stest, bool verbose = false);


public:
	// become PTA
	RESULT becomePrefixTree(Sample &S);
    // Apprentissage de PRFA, renvoie un PRFA prefixe
    // prec prend la precision, et mode est le mode d'apprentissage
    RESULT DEES(T_ModeVariables modeVariables, // the learning space for variables
				Sample &S,     // The learning sample (data)
                const double prec,  // The precision parameter
                const double epsprime, // The bound under which transition are erased
                const bool verbose = false,  // verbose mode (show the states during construction)
                const T_ModeReturn moderet = end, // return mode (end of the tree or begin of the tree)
                const T_ModeEpsilon modeeps = variable, // epsilon mode (epsfixed or variable)
                unsigned int maxstates = INT_MAX, // The maximal states number after which learning is cancelled
                unsigned int seuil=10,  // La periodicitÈ minimum avant de ne plus considÈrer l'Ètat
                int maxmots=10, // le nombre maximal de mots que l'on reconsidËre ‡ chaque Ètape
                int maxsearches=0, // le nombre de recherche dichotomique pour trouver le meilleur epsilon
                bool bestsearch=true, // vrai si on recherche  le meilleur epsilon
				bool stepssave=false); // vrai sauvegarde chaque etape de l'apprentissage

    // Apprentissage de PRFA, renvoie un PRFA prefixe
    // prend la precision, et mode est le mode d'apprentissage
    RESULT DEESBM(Sample &S,     // The learning sample (data)
                  const double prec,  // The precision parameter
                  const double epsprime, // The bound under which transition are erased
                  const bool verbose = false,  // verbose mode (show the states during construction)
                  const unsigned int maxstates = INT_MAX, // The maximal states number after which learning is cancelled
                  const unsigned int seuil=10,  // La periodicitÈ minimum avant de ne plus considÈrer l'Ètat
				  const double seuilbm=0.0001, // Le seuil en dessous duquel on considère avoir trouver le MLM
                  const unsigned int nb_tours=10); // Le nombre de tours pour BM

	Word word_of(const State x) const;


protected:
		
	// ##### Pour la version Baum Welch #####
    RESULT setRandomSolution(SFunc &solution, const StateSet &R) const ;
    RESULT BMsol (SFunc &solution, const Word &v,const Sample &S, const Sample &Sdep, const float prec,int nb_tours, bool verbose = false);
	// recupere la solution en fonction des valeurs recuperes apres Baum Welch
	RESULT get_solution(const PPRFA &A,SFunc &solution, Transition &t,TransitionFunction &T,const StateSFunc &Trad);
	// Ajoute les etats necessaires pour que tous les mots vS soient generables
	RESULT put_prefix_tree(const Sample &S, const Word &v);
	// Ajoute les etats necessaires pour que le mot w deviennent generable
	RESULT add_word(const State q, const Word &w);
	// Liste les successeurs d'elements de Q qui ne sont pas dans Q
	RESULT choice(list<State> &X, const StateSet Q) const;
	RESULT difference(const TransitionFunction T1, const TransitionFunction T2) const;
	// Initialise la solution en utilisant une resolution de systeme d'equation
	RESULT initSolution(SFunc &solution, const State x, const StateSet &R, const Sample &S) const;
	// Construit un quasi pta et teste sa qualite par rapport a l'echantillon
	bool TestQPTA(const Sample &S,const State x, const StateSet R,const double precision, const double seuilbm, bool verbose);
	// remet a jour la liste des mots associes aux etats
	RESULT compute_XR(void);

	};

#endif

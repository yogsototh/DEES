/***************************************************************************
                          ma.h  -  Multiplicity Automata
							-------------------
    begin                : Tue Jul 16 2002
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

/*_________________________________________________________________________*\
  Un automate fini fonctionnel se défini comme un quintuplet :
  Sigma : l'alphabet, un ensemble fini de lettre. 
  Ici, les lettre seront des mots.
  Q    : L'ensemble des états. Pour l'instant on se limite aux entiers.
  iota : une fonction de Q dans IR
  tau  : une fonction de Q dans IR
  phi  : une fonction de Q x Sigma x Q dans IR
 
  \***************************************************************************/

#ifndef MA_H
#define MA_H

#include "general.H"
#include "sample.H"
#include <math.h>

typedef set < State > StateSet;	// type Ensemble d'états
typedef map < State, float >SFunc;	// Type fonction State --> IR
typedef map < Lettre, SFunc > LettreSFunc;	// Type fonction Lettre x State --> IR
// typedef hmext::hash_map < State, LettreSFunc > Graph;	// Le type Graph
typedef map < State, LettreSFunc > Graph;	// Le type Graph
struct ordre_transition
{
    bool operator()(const Transition t1, const Transition t2) const
    {
        return
            (t1.qarr < t2.qarr) ||
            (
                (t1.qarr == t2.qarr) &&
                (
                    (t1.a < t2.a) ||
                    (
                        (t1.a == t2.a) && (t1.qdep < t2.qdep)
                    )
                )
            );
    }
};
typedef map<Transition, float, ordre_transition> TransitionFunction;


typedef enum { paffa, paalergia, pamdi} FormatSauvegarde;

class MA
{
    // ----------------  Structure de donnée ------------------------
public:
    // Le nom du MA
    string name;
protected:
    // Le dictionnaire
    Dictionnaire alph;	// permet d'associer des chaines de caractères aux lettre
    // l'alphabet
    Alphabet Sigma;
    // L'ensemble d'états
    StateSet Q;
    // La fonction iota
    SFunc iota;
    // La fonction tau
    SFunc tau;
    // La fonction phi;
    Graph phi;

    // ------------------------- Méthodes ---------------------------
protected:		// méthodes internes
    // renvoie un flottant aleatoire compris entre min et max
    float random (float min = 0, float max = 1) const;
    // renvoie vrai si x est à peu près égal à 1;
    inline bool simeqone (double x) const
    {
        return ((x > 0.9999) && (x < 1.0001));
    }

    // --- Constructeurs et destructeurs --
public:
    MA (int nb_lettres = 0, int nb_etats = 0);
    ~MA ();

    // --- Affichage ---
    // Renvoie la chaine de caractère d'affichage
    string message (void) const;
    // Affichage sous la sortie standart
    RESULT affiche (void) const;

    // --- Input output ---
    // Sauvegarde le MA dans Filename
	RESULT save (ofstream &fp, const FormatSauvegarde c_format = paffa) const;
	RESULT save (const char *Filename, const FormatSauvegarde format = paffa) const;
    inline RESULT save (const string Filename) const {return save (Filename.c_str ());}
    // Charge le MA de Filename
    RESULT load (const char *Filename);
    RESULT load (string Filename) {return load (Filename.c_str ());}

    // save in dot format
    RESULT save_dot (const char *Filename,
                     const int precision = 2,
                     const bool noms_etats = false,
                     const bool multiple_trans = false) const;
    RESULT save_dot (const string Filename,
                     int precision = 2,
                     bool noms_etats = false,
                     bool multiple_trans = false) const;

    // --- Accesseurs ---
	// true if the MA is empty
	inline bool empty() const {return Q.empty();}
	// return the size of the MA
    inline unsigned long size () const {return Q.size ();}
	// return the number of letter of the MA
    inline unsigned long nb_lettres () const {return Sigma.size ();}
	// return safely iota[q] (even if iota[q]=0)
    inline float val_iota (const State q) const {
        SFunc::const_iterator i = iota.find (q);
        return i == iota.end () ? 0 : i->second;
    }
	// return safely tau[q] (even if tau[q]=0)
    inline float val_tau (const State q) const{
        SFunc::const_iterator i = tau.find (q);
        return i == tau.end () ? 0 : i->second;
    }
	// return safely phi[t] (even if phi(t)=0 )
    inline float val_phi (const Transition &t) const {
        return val_phi(t.qdep, t.a, t.qarr);
    }
	// return safely phi[t] (even if phi(t)=0 )
    inline float val_phi (const State q, const Lettre a, const State s) const {
        Graph::const_iterator qi = phi.find (q);
        if (qi == phi.end ()) return 0;
        LettreSFunc::const_iterator ai = qi->second.find (a);
        if (ai == qi->second.end ()) return 0;
        SFunc::const_iterator si = ai->second.find (s);
        if (si == ai->second.end ()) return 0;
        return si->second;
    }

	// The transitionFunction T is the set of all transitions of the MA
	void allTransitions(TransitionFunction &T) const;
	// The stateFucntion S is the set of all state of the MA
	void allStates(SFunc &S) const;

	// the delta function
    StateSet delta(const State q, const Word &w) const;
    StateSet delta(const StateSet &R, const Word &w) const;
	
	// direct read only accessor
	inline Dictionnaire dictionnaire(void) const {return alph;};

    // ---  Modifieurs ---
    // Ajout une lettre par defaut à l'alphabet
    Lettre addLettre (void);
    // Add the letter buf
    Lettre addLettre (const TrueLettre & buf);
    // Ajoute un état
    //State addState (const float init = 0, const float term = 1);
    State addNewState (const float init, const float term);

    // Ajout l'état q s'il n'existe pas déjà
    State addState (const State q, const float init = 0, const float term = 1);

    // Fonction qui ajoute une transition
    RESULT addTransition (const Transition & t, const float val);
    // Ajout d'une transition
    RESULT addTransition (const State qdep, const Lettre a, const State qarr, const float val);
    // erase badstate
    RESULT erase (State badstate);
    // clean transitions with parameter between minval and maxval
    RESULT erase_transitions (const double maxval = 0, const double minval = 0);
	// erase cleanly the transition phi(qbeg, l ,qend)
	RESULT eraseTransition(const State qbeg, const Lettre l, const State qend);
	
    // Become a Random MA
    RESULT becomeRandom (const int nb_etats,	// le nombre d'états
                         const int nb_lettres,	// le nombre de lettres
                         const int num_graphe = 0,	// le numero du graphe
                         const float densite = 0.3,	// densité du ffa (0 ffa vide, 1 ffa complet)
                         const float prob_init = 0.3,	// la probabilité pour un état d'etre initial
                         const float prob_term = 0.3,	// probabilité pour un état d'être terminal
                         const float min_trans = 0,	// la valeur minimale des transitions
                         const float max_trans = 1);	// la valeur maximale des transitions
						 
    // Cree un MA prefixe aleatoire
    RESULT becomeRandomPrefix (const int nb_etats,	// le nombre d'états
                               const int nb_lettres,	// le nombre de lettres
                               const int num_graphe = 0,	// le numero du graphe
                               const float densite = 0.3, // densité du ffa (0 ffa vide, 1 ffa complet)
                               const float prob_init = 0.3,	// la probabilité pour un état d'etre initial
                               const float prob_term = 0.3,	// probabilité pour un état d'être terminal
                               const float min_trans = 0,	// la valeur minimale des transitions
                               const float max_trans = 1);	// la valeur maximale des transitions
							   
    // Cree un MA aleatoire avec un nombre d'aretes maximal fixé
    RESULT becomeRandomMax (const int nb_etats, const int nb_lettres,
                            const int num_graphe=0, const int nb_succ=2,
                            const int nb_init=1,    const int nb_term=INT_MAX,
                            const float min_trans=0,  const float max_trans=1);
							
	// Cree un MA aleatoire avec le maximum de controle
	// nb_etats : number of states of the MA
	// nb_lettres : number of letters of the alphabet of the MA
	// num_graphe : if different of 0, the same number gives always the same MA
	// max_nb_init : the maximal number of initial states
	// max_nb_succ : the maximal number of successors by states-letter
	// max_nb_term : the maximal number of terminal states
	// min_trans : minimal value of transitions
	// max_trans : maximal value of transitions
	// min_iota : minimal value for iota(q)
	// max_iota : maximal value for iota(q)
	// min_tau : minimal value for tau(q)
	// max_tau : maximal value for tau(q)
	// nomalize : if different of 0, normalize sum of iota and sum of tau + successor transitions
	// prob_init : probability for a state to be initial
	// prob_trans : probability to create each transition
	// prob_term : probability for a state to be terminal
	RESULT becomeRandomControl (const int nb_etats, const int nb_lettres, const int num_graphe, const int max_nb_succ, const int max_nb_init, const int max_nb_term, const float min_trans,  const float max_trans, const float min_iota, const float max_iota, const float min_tau, const float max_tau, const float normalization, const float prob_init, const float prob_trans, const float prob_term);
		
	
	RESULT becomeComplete (const int nb_etats, const int nb_lettres);

    // Fonction qui rend le MA vide (pas d'états pas de lettre)
    void vide ();
	// exactly as vide, clear the MA
	inline void clear() { vide(); }
	// emmonde the MA
	// delete all reachable states and all states from which no state can be reached
	RESULT emmonde(void);


	// --- Fonctions dŽpendant de l'objet ---
    // the next word using the alphabet Sigma
    void next (Word & w) const;
    // renvoie la chaine correspondant au mot
    string affiche (const Word & w) const;
    void normaliseVecteur(SFunc &V) const;
	RESULT val_PSe(SFunc &PSe) const;
	// generate a word from the model using the P_r function where r
	// is the rationnal serie of the MA
    // genere un mot en corelation avec la distribution P_r associe
    // au MA
    RESULT genere_mot(Word &w, const bool safe = true);
    // genere un échantillon de mots
    RESULT genere_echantillon(const int taille, Sample &S, const int num_echantillon=0);
    // genere un échantillon de mots
    RESULT genere_echantillon(const int taille, Sample &S, const char *Filename, const int num_echantillon=0);

    // --- Test ---
    // Teste la cohérence de la structure de donnée
    bool coherent (void) const;
    // Teste si le MA est coherent et est un SPFA. valeur de iota, tau et phi dans [0,1]
    bool isSPFA (void) const ;
    // renvoi vrai si le MA est un PFA
    bool isPFA (void) const;
	
	private:
	RESULT renormalize_state(const State q, const float renormalization, const float min_tau, const float max_tau, const float min, const float max);

	
};

#endif

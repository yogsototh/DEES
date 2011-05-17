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

/*________________________  DÈscription plus prÈcise ______________________*\
  Classe qui hÈrite de la classe spfa, la seule diffÈrence est que l'on
  peut gÈnÈrer des mots avec une distribution de probabilitÈ.
  \***************************************************************************/

#ifndef PFA_H
#define PFA_H

#include "spfa.H"
#include "sample.H"
//#include "simplex.H"


// ------------------------ La Classe ----------------------------------

class PFA : public SPFA
{
public:
    RESULT rend_PFA(void);

    // cree un PFA aleatoire
    RESULT becomeRandom(int num_etat, int num_lettre, // le nombre d'etats et le nombre de lettres
                        int num_graphe = 0,           // le numero du graphe (donner deux fois le mÍme numero avec
                        // les memes paramËtres donne le mÍme graphe
                        float densite=0.3,  // densitÈ du ffa (0 ffa vide, 1 ffa complet)
                        float prob_init=0.3, // probabilitÈ d'Ítre initial
                        float prob_term=0.3, // probabilitÈ d'Ítre terminal
                        float min_trans=0,
                        float max_trans=1);   // valeur maximale des transitions

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
                            const int num_graphe=0, const int nb_succ=2,
                            const int nb_init=1,    const int nb_term=INT_MAX,
                            const float min_trans=0,  const float max_trans=1);

    // genere un mot en corelation avec la distribution
    RESULT genere_mot(Word &w, const bool safe = true);
    // -- dans MA maintenant -- // genere un Èchantillon de mots
    RESULT genere_echantillon(const int taille, Sample &S, const int num_echantillon=0);
    // -- dans MA maintenant -- // genere un Èchantillon de mots
    RESULT genere_echantillon(const int taille, Sample &S, const char *Filename, const int num_echantillon=0);

    // transforme un MA avec valeurs positives en SPFA
    RESULT rend_spfa(void);

    // Lissage
    RESULT lisse(int mode=0, float delta=0.01);

    // become the PTA representing S
    RESULT becomePrefixTree(Sample &S);

    // ----------- DISTANCES ENTRE PFAs -----------

    // renvoie la perplexitÈ du PFA par rappor ‡ l'echantillon S
    double perplexite(const Sample &S) const;
    // renvoie la divergence de Kullback Leibler relativemant ‡ l'Èchantillon S
    double divKL(const PFA &B, const Sample &S);
	// renvoie la divergence de Kullback Leibler relativemant ‡ l'Èchantillon S en utilisant B
	double divKL (PFA &B, const Sample & S) const;
	// renvoie la norme 2 entre ce PFA et le PTA correspondant ‡ S
	double d_2(const Sample &S) const;
    // renvoie la norme 2 entre ce PFA et un autre
    double d_2(const PFA &B, const Sample &S) const;
    // revoie la norme L1 entre ce PFA et un autre en fonction de l'Èchantillon
    double L_1(const PFA &B, const Sample &S) const;
    // renvoie la norme 2 entre ce PFA et un autre calculée sans passer par les logs
    double d_2nlog(const PFA &B, const Sample &S) const;
    // revoie la norme L1 entre ce PFA et un autre en fonction de l'Èchantillon calculé sans passer par les logs
    double L_1nlog(const PFA &B, const Sample &S) const;
    // renvoie la moyenne des Ècarts de log en fonction d'un Èchantillon S
    double dlog(const PFA &B, const Sample &S) const;
    // renvoie dlog sans prendre en compte les mots non gÈnÈrables par l'un des deux PFAs
    double dlog_safe(const PFA &B, const Sample &S) const;
    // renvoie la vraissemblance
    float Likelihood(const Sample &S);


    // return true if the MA is a PRFA
    bool isPRFA(bool verbose=false);

    // return true if the PFA is a PDFA
    bool isPDFA();

};

#endif

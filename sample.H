/***************************************************************************
                       echantillon.h - Sample de mots
                             -------------------
    begin                : 20 Jan 2003
    copyright            : (C) 2003 by Yann Esposito
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
  Classe échantillon, échantillon de mots.
  \***************************************************************************/

#ifndef ECHANTILLON_H
#define ECHANTILLON_H

// le type de sortie
typedef enum {ffa, alergia, mdi} T_Format;

#include "general.H"
#include <list>
// ----------------------- Declaration ---------------------------------
// ordre naturel sur les mots
struct ordre_mot
{
    bool operator()(const Word u, const Word v) const
    {
        if (u.size() < v.size())
            return true;
        else if (u.size() > v.size())
            return false;
        else
        {
            less<Word> infeq;
            return infeq(u,v);
        }
    }
};

// Type echantillon qui suit l'ordre naturel sur les mots
typedef map<Word, int, ordre_mot> EnsWord;  // Le type ensemble de mots

// ------------------------ La Classe ----------------------------------
class Sample
{
protected:
    Alphabet Sigma;      // Ensemble de lettre de l'echantillon
    Dictionnaire alph;   // permet d'associer des éléments plus complexes aux lettres.
    EnsWord S;            // l'ensemble de mots
    unsigned int taille; // la taille de l'echantillon
    bool prefixiel;      // vrai si l'echantillon est prefixiel.
public:
    int seuil;  // Le seuil à partir duquel on permet de visualiser un élément de l'échantillon

public:

    // --- Constructeur et Destructeurs ---
    Sample() : taille(0), seuil(0)
    {}
    Sample(const Alphabet A, const Dictionnaire D, const bool pref=false, int n_seuil=0) :
            Sigma(A), alph(D), taille(0), prefixiel(pref), seuil(n_seuil)
    {}
    ~Sample(void)
    {}

    // --- Les accesseurs ---
    // Renvoie le nombre de fois qu'est present le mot w
    inline unsigned int count(const Word w) const
    {
		EnsWord::const_iterator s;
		s = S.find(w);
        return (s!=S.end())?s->second:0;
    }
    inline unsigned int operator[](const Word w) const
    {
        return count(w);
    }
    // Renvoie la taille de l'echantillon
    inline unsigned int size() const
    {
        return taille;
    }
    // Renvoie l'alphabet de l'echantillon
    inline Alphabet alphabet() const
    {
        return Sigma;
    }
    // Renvoie le dictionnaire
    inline Dictionnaire dictionnaire() const
    {
        return alph;
    }
    // Renvoie l'ensemble de mot
    inline EnsWord ensWord() const
    {
        return S;
    }
    // Renvoie le nombre de mots contenus dans S (pour la fct° precision)
    inline unsigned long nb_mots() const
    {
        return S.size();
    }
    // Renvoie alph
    inline Dictionnaire get_alph() const
    {
        return alph;
    }

    // Revoie vrai si l'echantillon est prefixiel
    inline bool isprefixiel() const
    {
        return prefixiel;
    }

    // iterateurs
    typedef EnsWord::const_iterator const_iterator;
    inline const_iterator begin() const
    {
        return S.begin();
    }
    inline const_iterator end() const
    {
        return S.end();
    }
    inline const_iterator find(Word w) const
    {
        EnsWord::const_iterator u=S.find(w);
        return ((u!=S.end()) && (u->second >= seuil))?u:S.end();
    }

public:
    // renvoie le prochain iterateur pointant vers un mot apparaissant
    // un nombre supérieur de fois au seuil.
    inline void next(const_iterator &w) const
    {
        do
            w++;
        while ((w->second < seuil) && (w != S.end()));
    }

    // --- Les méthodes ---
    // Ajout une lettre par defaut à l'alphabet
    Lettre addLettre(const Lettre l);
    // Ajout de la lettre specifiee à l'alphabet
    Lettre addLettre(const string l);
	// Ajout de la lettre specifie a l'alphabet c etant la lettre associe pour les mots
	Lettre addLettre (const Lettre c, const string l);
    // Insere le mot w, si safe est vrai, on met à jour l'alphabet et le dictionnaire.
    RESULT insert(const Word w, const unsigned int num=1, const bool safe=true);
    // Rend l'echantillon vide
    RESULT vide();
    inline RESULT clear()
    {
        return vide();
    };
    // Sauvegarde l'echantillon dans le fichier Filename
    RESULT save(const char *Filename, T_Format format=ffa) const;
    inline RESULT save(const string Filename, T_Format format=ffa) const
    {
        return Sample::save(Filename.c_str(), format);
    }
    // Charge l'echantillon à partir d'un fichier
    RESULT load(const char *Filename);
    inline RESULT load(const string Filename)
    {
        return load((char *)Filename.c_str());
    }
    // Prefixialise l'echantillon
    RESULT prefixialise();
    // Deprefixialise l'echantillon
    RESULT deprefixialise();
    // Suppression des valeurs inferieures au seuil
    inline RESULT seuille(const int n_seuil)
    {
        return seuil = n_seuil;
    }
    // Sépare un échantillon en un échantillon principal et un échantillon test
    int separe(Sample &Stest, const float proportion);
    // Ajoute un échantillon à un autre
    Sample &operator+=(const Sample &Stest);

    // return the likelihood of the MLM
    float AutoLikelihood() const;
	
	// return the list of words successor of v (in the prefixial sample)
	Sample prefix_successor_of(const Word &v) ;
	// return words begining by v
	Sample begin_by(const Word &v);
	// become S inter S2
	void inter(const Sample &S2);

};

#endif

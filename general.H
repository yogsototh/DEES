/***************************************************************************
             general.h  -  This file must be included in all Project
                             -------------------
    begin                : Thu Jul 18 2002
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

#ifndef GENERAL_H
#define GENERAL_H

using namespace std;
#include <string>
#include <set>
#include <map>
//#include <vector>
//#include <hash_set>
#include <iostream>		// pour cout <<
#include <fstream>		// pour fic <<
#include <cmath>

// ------------- pour le portage --------------

// pour compiler sur le ppc

#ifdef __ppc__
#define isnan __isnan
#define isinf __isinf
#endif

// pour utiliser les hash_map
#ifdef __GNUC__
      #if __GNUC__ < 3
        #include <hash_map.h>

namespace hmext
{
using ::hash_map;
}
; // inherit globals
#else
        #include <ext/hash_map>
        #if __GNUC_MINOR__ == 0

namespace hmext = std;               // GCC 3.0
#else

namespace hmext = ::__gnu_cxx;       // GCC 3.1 and later
#endif
      #endif
      #else      // ...  there are other compilers, right?

namespace hmext = std;
#endif

// ------------------------------------------------

// Make DEBUG false to disable the debug mode (faster but you'll lose all error messages).
const bool PFA_DEBUG = true;

// Make VERBOSE false to disable the verbose mode (less messages)
const bool PFA_VERBOSE = true;

// Make SAFE false to disable the SAFE mode (faster but DANGEROUS !!!!)
const bool PFA_SAFE = true;

// Definition of the RESULT type.
typedef int RESULT;

// Test if the result of a function is OK.
inline RESULT
OK (RESULT x)
{
    return x >= 0;
}
inline RESULT
VAL (RESULT x)
{
    return x;
}
inline RESULT
ERR (RESULT x)
{
    if (x >= 0)
        return -1 - x;
    else
        return -1 + x;
}

// ===============  Définition des types necessaires aux structures de données ==============
// ___   Située hors des classes pour éviter l'utilisation de l'opérateur de résolution   ___

typedef int State;		// Le type état
typedef char Lettre;		// Le type lettre
typedef basic_string< Lettre > Word; // Le type Word
//~ class Word : public basic_string < Lettre >	// Le type mot
//~ {
//~ public:
//~ inline Lettre pop_back() {
//~ iterator e=end();
//~ --e;
//~ Lettre a=*e;
//~ erase(e);
//~ return a;
//~ }
//~ };

typedef set< Lettre > Alphabet;	// type Alphabet
typedef string TrueLettre;	// real letters.
//typedef hash_map < Lettre, TrueLettre > Dictionnaire;	// type Dictionnaire
typedef map < Lettre, TrueLettre > Dictionnaire;	// type Dictionnaire (association lettre et vraie lettre)

//typedef hmext::hash_map < State, double >PreciseSFunc;	// Type fonction State --> IR, plus précise
typedef map < State, double >PreciseSFunc;	// Type fonction State --> IR, plus précise

// Le type transition n'est pas forcement necessaire mais il est joli
// Type transition est un triplet State x Lettre x State
typedef struct str_Transition
{
    State qdep;
    Lettre a;
    State qarr;
}
Transition;

#endif

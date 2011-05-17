/***************************************************************************
                               echantillon.C
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

/*
  Pour les notes, voir le fichier "echantillon.h".
  See "echantillon.h" file to view notes.
*/
#include "sample.H"
#include <algorithm>

// Ajout une lettre par defaut à l'alphabet
Lettre Sample::addLettre (const Lettre l)
{
    Lettre a;
    if (!Sigma.empty ())
    {
        a = *(--(Sigma.end ())) + 1;
    }
    else
    {
        a = 'a';
    }
    Sigma.insert (a);
    alph[a] = l;
    return a;
}

// Ajout de la lettre spécifiée à l'alphabet
Lettre Sample::addLettre (const string l)
{
    Lettre a;
    if (!Sigma.empty ())
    {
        a = *(--(Sigma.end ())) + 1;
    }
    else
    {
        a = 'a';
    }
    Sigma.insert (a);
    alph[a] = l;
    return a;
}

// Ajout de la lettre specifiee a l'alphabet
Lettre Sample::addLettre (const Lettre c, const string l)
{
    Sigma.insert (c);
    alph[c] = l;
    return c;
}

// Insere num fois le mot w
RESULT Sample::insert (const Word w,
                            const unsigned int num,
                            const bool safe)
{
    if (safe)
    {
        // on verifie si le mot contient des lettres en plus
        Word::const_iterator a;
        for (a = w.begin (); a != w.end (); a++)
        {
            if (Sigma.find (*a) == Sigma.end ())
            {
                addLettre (*a);
            }
        }
    }

    if (prefixiel)
    {
        Word	u;
        Word::const_iterator a;
        u.clear ();
        for (a = w.begin (); a != w.end (); a++)
        {
            u += *a;
            S[u] += num;
        }
    }
    else
    {
        S[w] += num;
    }

    taille += num;

    return VAL (0);
}

// Rend l'echantillon vide
RESULT Sample::vide ()
{
    S.clear ();
    Sigma.clear ();
    alph.clear ();
    taille = 0;
    prefixiel = false;
    return VAL (0);
}

// Sauvegarde un echantillon dans le fichier Filename
RESULT Sample::save (const char *Filename, T_Format format) const
{				
    try
    {
        if (Filename == NULL)
            throw -1;
        ofstream
        fp (Filename);
        if (!fp.good ())
            throw -2;

        EnsWord::const_iterator w;
        Word::const_iterator l;
        Dictionnaire::const_iterator a;
        string msg;
        int i;

        switch (format)
        {
        case ffa:	// ========================= format ffa =============================
            // Affichage de l'entete
            fp << "Sample" << endl;	// Affichage de Sample
            // Affichage de la taille
            fp << "Taille " << size () << endl;
            // ---- ne pas modifier l'ordre Nb_Lettre puis Alphabet pour le "load" ---
            // Affichage du nombre de lettres
            fp << "Nb_Lettres " << alph.size () << endl;
            // Affichage de l'alphabet
            fp << "Alphabet" << endl;
            for (a = alph.begin (); a != alph.end (); a++)
            {
                fp << a->first << " " << a->second << "\n";
            }
            fp << endl;
            // Affichage du drapeau prefixiel
            fp << "Prefixiel " << prefixiel << endl;
            // Affichage de l'echantillon proprement dit
            for (w = S.begin (); w != S.end (); w++)
            {
                for (l = w->first.begin ();
                        l != w->first.end (); l++)
                {
                    fp << alph.find (*l)->
                    second << " ";
                }
                fp << ": " << w->second;
                fp << endl;
            }
            break;

        case alergia:	// format alergia
            for (w = S.begin (); w != S.end (); w++)
            {
                msg="";
                for (l = w->first.begin (); l != w->first.end (); l++)
                {
                    msg += alph.find (*l)->second;
                }
                for (i = 0; i < w->second; i++)
                    fp << msg << endl;
            }
            break;

        case mdi:	// format MDI
            for (w = S.begin (); w != S.end (); w++) {
                msg.clear ();
                for (l = w->first.begin ();  l != w->first.end (); l++) {
                    msg += alph.find(*l)->second;
                    msg += " ";
                }
                for (i = 0; i < w->second; i++) {
                    fp << msg << endl;
				}
            }
            break;

        default:
            cerr << "Sample::save()" << endl;
            cerr << "Format " << format << "inconnu !!!" <<
            endl;
            throw - 3;
        }
        return VAL (0);
    }
    catch (int e)
    {
        if (PFA_VERBOSE)
        {
            cerr << " PFA::echantillon(int taille, char *Filename)" << endl;
            switch (e)
            {
            case -1:
                cerr << " : Filename = NULL !!!" << endl;
                break;
            case -2:
                cerr << "Impossible d'ecrire dans le fichier '" << Filename << "' !!!" << endl;
            default:
                cerr << "Unknown error" << endl;
            }
        }
        return e;
    }
    catch (...)
    {
        cerr << " PFA::echantillon(int taille, char *Filename)"
        << endl;
        cerr << "Unknown error !!!" << endl;
        return -1;
    }
}

// Charge un echantillon
RESULT Sample::load (const char *Filename)
{
    // les échantillons sont sous la forme d'un liste de lignes
    // mot : nb_apparitions
    try
    {
        // on verifie qu'il n'y ai pas une erreur du passage du nom
        if (Filename == NULL)
            throw - 1;
        ifstream
        fp (Filename);	// on declare le descripteur de fichier et on l'ouvre
        // on verifie qu'il n'y ait pas d'erreur d'ouverture
        if (!fp.good ())
            throw - 2;


        Word w;	// le mot
        map < string, Lettre > inv_alph;	// le dictionnaire inversé
        int i;
        int nb_lettres;
        bool	safe;	// des drapeaux fin de fichier et mode safe
        char	c;
        string buf;	// le buffer

        // on vide la structure
        vide ();

        fp >> buf;	// on lit le premier mot
        if (buf == "Sample")
        {		// ========================== TYPE de fichier echantillon interne =

            safe = true;
            // Tant que l'on est pas à la fin du fichier
            while (!fp.eof ())
            {	// boucle principale
                fp >> buf;
                if (buf == "%")
                {	// -------------------- Commentaire
                    while ((fp.good ()) && (c != '\n'))
                    {
                        fp >> c;
                    }
                }
                else if (buf == "Taille")
                {	// ------- Taille de l'echantillon
                    fp >> taille;	// on lit la taille
                }
                else if (buf == "Nb_Lettres")
                {	// --- On lit l'alphabet
                    fp >> nb_lettres;	// on lit le nombre de lettres
                    fp >> buf;
                    if (buf != "Alphabet")
                        throw - 3;
					Lettre c;
                    for (i = 0; (i < nb_lettres) && (fp.good ()); i++)
                    {
						fp >> c;
                        fp >> buf;	// on lit la lettre
                        addLettre (c,buf);
                    }
                    Dictionnaire::iterator x;
                    for (x = alph.begin ();
                            x != alph.end (); x++)
                    {
                        inv_alph[x->second] =
                            x->first;
                    }
                    safe = false;	// pas besoin d'etre en mode safe pour le chargement de l'echantillon
                }
                else if (buf == "Prefixiel")
                {	// ------ On regarde si l'echantillon est prefixiel
                    fp >> prefixiel;	// on lit la valeur du drapeau "prefixiel"
                }
                else
                {	// ------------------- On charge les données
                    while (fp.good ())
                    {
                        // on lit le mot
                        w.clear ();
                        while ((fp.good ())
                                && (buf != ":"))
                        {
                            if (safe)
                            {	// si on ne connait pas l'alphabet, on verifie l'arrivee de nouvelles lettres.
                                if (inv_alph.
                                        find (buf)
                                        ==
                                        inv_alph.
                                        end ())
                                {
                                    inv_alph[buf] = addLettre (buf);
                                }
                            }
                            w += inv_alph.
                                 find (buf)->
                                 second;
                            //            w.push_back(inv_alph.find(buf)->second);
                            fp >> buf;	// on lit la prochaine lettre.
                        }
                        // on a fini de lire le mot
                        if (fp.good ())
                            fp >> S[w];	// on lit le nombre de fois qu'est apparu ce mot.
                        fp >> buf;
                    }
                }
            }
            return VAL (0);
        }
        else if (buf == "Text_without_space")
        {		// ================================ Mode mots (alergia et mdi)-
            w.clear ();	// on initialise le mot courant.
            prefixiel = false;
            fp.get(c); // on lit le caractere '\n'
            while (!fp.eof())
            {
                fp.get(c);
                // on ajoute si nécessaire les nouvelles lettres
                if ((c != '\n') && (c != '\t') && (c != ' ')) // on lit une nouvelle lettre
                {
                    if (Sigma.find (c) == Sigma.end ()) // si c \notin Sigma
                    {
                        Sigma.insert (c);
                        alph[c] = c;
                    }
                    w += c;	// on ralonge le mot
                } else // on est a la fin de la phrase
                {
                    ++S[w];	// on ajoute le mot un coup
                    w.clear();	// et on le réinitialise
                }
            }

            // on supprime l'artefact pour epsilon
            w.clear();
            --S[w];
            if (S.begin()->second == 0)
                S.erase(S.begin());

            EnsWord::iterator wi;
            taille=0;
            for (wi=S.begin() ; wi != S.end() ; ++wi)
            {
                taille += wi->second;
            }

            return VAL (0);
        }
        else if (buf == "Text_with_spaces")
        {		// ========================== mode phrases ---
            //caractere de separation : .
            Word curmot;
            Lettre curlettre = 'a';

            curmot.clear ();
            prefixiel = false;
            Dictionnaire::const_iterator d;
            while (fp.good ())
            {
				bool findephrase=false;
                fp >> buf;	// on lit la premiere lettre (un mot)
                //      cout << "lettre : " << buf << flush;
                c = *(buf.rbegin ());	// c is the last letter of the string buf
                if ((c == '.') || (c == '!')  || (c='?')) {
					buf.erase(--buf.end());
					findephrase=true;
				}
                // on cherche a savoir si on a deja insere w
                for (d = alph.begin ();(d != alph.end ()) && (d->second != buf); d++)
                    ;
                if (d == alph.end ())
                {
                    //      cout << ": jamais vue." << endl;
                    Sigma.insert (curlettre);
                    alph[curlettre] = buf;
                    curmot += curlettre;
                    ++curlettre;
                }
                else
                {
                    //      cout << ": deja vue." << endl;
                    curmot += d->first;
                }
                c = *(buf.rbegin ());	// c is the last letter of the string buf
                if (findephrase)
                {
                    // si on arrive en fin de phrase
                    cout << "fin de la phrase : " << endl;
                    for (Word::const_iterator mi = curmot.begin() ;  mi != curmot.end() ; ++mi)
						cout << alph[*mi] << " ";
                    cout << endl;
                    ++S[curmot];
                    ++taille;
                    curmot.clear ();
					findephrase=false;
                }
            }
            return VAL (0);
        }
        else if (buf == "Amaury")
        {		// ========================== mode phrases pour Amaury ---
            //caractere de separation : .
            Word curmot;
            // Lettre curlettre = 'a';

            curmot.clear ();
            prefixiel = false;
            Dictionnaire::const_iterator d;
            while (fp.good ())
            {
                ++S[curmot];
                ++taille;
                curmot.clear ();
            }
            return VAL (0);
        }
        else
            throw - 4;
    }
    catch (int e)
    {
        if (PFA_VERBOSE)
        {
            cerr << " Sample::load(char *Filename)" << endl;
            switch (e)
            {
            case -1:
                cerr << "Filename = NULL !!!" << endl;
                break;

            case -2:
                cerr << Filename << "is unreadable !!!" <<
                endl;
                break;

            case -3:
                cerr << "Erreur de coherence lors du chargement du fichier :" << Filename << endl;
                cerr << "\"Nb_Lettres\" n'est pas suivit d'un chiffre puis de \"Alphabet\"" << endl;
                break;

            case -4:
                cerr << "Format of the file " << Filename <<
                " is not recognized" << endl;
                cerr << "Add 'Text_without_space' or 'Text_with_spaces' at the begining of the file" << endl;
                cerr << "to make the file in a valid format."
                << endl;
                break;

            default:
                cerr << "Erreur non repertoriée" << endl;
            }
        }
        return e;
    }
    catch (...)
    {
        cerr << "Sample::load() Erreur non repertoriee et non traitee" << endl;
        return -1;
    }
}

// Genere l'echantillon prefixe associe à l'échantillon non prefixe S.
RESULT Sample::prefixialise ()
{
    if (!prefixiel)
    {
        Word  u;
        Word::const_iterator a;	// la lettre courante
        EnsWord::const_iterator w;	// le mot courant
        EnsWord Spref;	// l'echantillon prefixe

        u.clear ();
        Spref[u] = taille;
        for (w = S.begin (); w != S.end (); ++w)
        {
            u.clear ();
            for (a = w->first.begin (); a != w->first.end (); ++a)
            {
                u += *a;
                Spref[u] += w->second;
            }
        }
        S = Spref;
        prefixiel = true;

        return VAL (0);
    }
    else
    {
        return ERR (0);
    }
}

// Deprefixialise
RESULT Sample::deprefixialise ()
{
    if (!prefixiel)
    {
        return VAL (1);
    }
    else
    {
        EnsWord R;
        EnsWord::const_iterator w, u;
        Alphabet::const_iterator a;
        int sum;
        for (w = S.begin (); w != S.end (); w++)
        {
            sum = 0;
            for (a = Sigma.begin (); a != Sigma.end (); a++)
            {
                u = S.find (w->first + *a);
                if (u != S.end ())
                    sum += u->second;
            }
            if (w->second != sum) {
                R[w->first] = w->second - sum;				
			}
        }
        S = R;
        prefixiel = false;
        return VAL (0);
    }
}

// Sépare un échantillon en un échantillon principal et un échantillon test
int
Sample::separe (Sample & Stest, const float proportion)
{
	deprefixialise();
	Stest.clear();
	Stest.alph=alph;
	Stest.Sigma=Sigma;
	int nb_suppr = int(float(taille) * proportion);
	int num;
	EnsWord::iterator w;

	while (nb_suppr-- != 0) {
		w=S.begin();
		num = int(float(rand ())*float(taille) / float(INT_MAX));
		num -= w->second;
		while (num > 0) {			
			w++;
			num -= w->second;	
		}
		Stest.insert( w->first );
		if (--(w->second) == 0) {
			S.erase(w);			
		}
		--taille;
	}

	return Stest.size();
	


    //~ deprefixialise ();
    //~ EnsWord::iterator w;
    //~ int i;
    //~ int nb_suppr;
    //~ Stest.clear ();		// first we clear the test sample
    //~ Stest.alph = alph;
    //~ Stest.Sigma = Sigma;

    //~ for (w = S.begin (); w != S.end (); ++w)
    //~ {
        //~ nb_suppr = 0;
        //~ for (i = 0; i < w->second; ++i)
        //~ {
            //~ float val;
            //~ val = (float) ((float) random () / (float) INT_MAX);
            //~ if (val < proportion)
            //~ {
                //~ Stest.insert (w->first);
                //~ ++nb_suppr;
            //~ }
        //~ }
        //~ w->second -= nb_suppr;
        //~ taille -= nb_suppr;
    //~ }
    //~ // Suppression des mots effaces de l'echantillons
    //~ EnsWord::iterator wtmp;
    //~ for (w = S.begin (); w != S.end (); ++w)
    //~ {
        //~ if (w->second == 0)
        //~ {
            //~ wtmp = w;
            //~ ++wtmp;
            //~ S.erase (w);
            //~ w = wtmp;
        //~ }
    //~ }
    //~ return Stest.size ();
}

// Ajoute un échantillon à un autre
Sample & Sample::operator+= (const Sample & Stest)
{
    //  if (alph != Stest.alph) throw 1;
    Sample::const_iterator w;
    for (w = Stest.begin (); w != Stest.end (); ++w)
        insert (w->first, w->second);
    return *this;
}

float Sample::AutoLikelihood() const
{
    if (prefixiel)
    {
		cerr << "Sample::AutoLikelihood() : echantillon prefixiel !" << endl;
		return -1;
    }

    Sample::const_iterator w;
    unsigned int si = size();
    float result=0;
    for (w=begin() ; w != end() ; ++w)
    {
        result += (log ( float(w->second) / float(si) )) * w->second;
	}
    return result;
}

Sample Sample::begin_by(const Word &v) {
	bool waspref=prefixiel;
	if (!prefixiel) {		
		prefixialise();
	}

	list<Word> W;
	Sample X;
	Word u;

	X.Sigma=Sigma;
	X.alph=alph;
	X.prefixiel=true;

	EnsWord::const_iterator w;
	w=S.find(v);
	if (w != S.end()) {
		W.push_back(w->first);
		X.S[w->first]=w->second;
	}
	
	Alphabet::const_iterator a;
	while (!W.empty()) {
		u=*(W.begin()); 
		W.pop_front();		
		for (a=Sigma.begin() ; a!=Sigma.end() ;++a) {
			u+=*a;
			w=S.find(u);
			if (w != S.end()) {
				W.push_back(w->first);
				X.S[w->first]=w->second;
			}
			u.erase(--u.end());
		}
	}

	X.deprefixialise();
	
	for (w=X.S.begin() ; w != X.S.end() ; ++w) {
		X.taille += w->second;	
	}

	if (!waspref) {
		deprefixialise();
	}

	return X;
}

Sample Sample::prefix_successor_of(const Word &v) {
	bool waspref=prefixiel;
	if (!prefixiel) {		
		prefixialise();
	}

	list<Word> W;
	Sample X;
	Word u;
	if (S.find(v) != S.end()) {
		W.push_back(u); X.insert(u);
	}
	
	Alphabet::const_iterator a;
	while (!W.empty()) {
		u=*(W.begin()); W.pop_front();		
		for (a=Sigma.begin() ; a!=Sigma.end() ;++a) {
			u+=*a;
			if (S.find(v+u) != S.end()) {
				W.push_back(u);
				X.insert(u);
			}
			u.erase(--u.end());
		}
	}

	if (!waspref) {
		deprefixialise();
	}

	return X;
}


void Sample::inter(const Sample &S2) {
	Sample::const_iterator u;
	EnsWord::iterator w;
	for (w = S.begin() ; w != S.end() ; ++w) {
		w->second = 0;
	}
	for (w = S.begin() ; w!=S.end() ; ++w) {
		u=S2.find(w->first);
		if (u != S2.end()) {
			w->second = u->second;
		}
	}

    // erase for iota
	EnsWord::iterator stmp;
    for (w=S.begin() ; w != S.end() ; )
    {
        if (w->second == 0)
        {
            stmp = w;
            ++stmp;
            S.erase (w);
            w = stmp;
        }
        else
            ++w;
    }
}

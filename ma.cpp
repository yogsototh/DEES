/***************************************************************************
                                 ma.C
                             -------------------
    begin                : Tue Jul 16 2002
    copyright            : (C) 2002 by Yann Esposito
    email                : esposito@cmi.univ-mrs.fr
***************************************************************************/

/*
  Pour les notes, voir le fichier "ffa.h".
  See "ffa.h" file to view notes.
*/

#include "ma.H"
#include <stdio.h>
#include <assert.h>
#include <sstream>
#include <algorithm>
#include "simplex.H"

RESULT MA::save_dot (const string Filename, int precision, bool noms_etats,
                      bool multiple_trans) const
{
    return save_dot (Filename.c_str (), precision, noms_etats, multiple_trans);
}


/// ------- MÈthodes internes ---------
float
MA::random (float min, float max) const
{
    if (PFA_SAFE)
    {
        assert (min <= max);
    }
    return ((((float) rand () + 1) / INT_MAX) * (max - min)) +
           min;
}

/// --------------------- Constructeurs et Destructeurs ----------------------
// Constructeur par dÈfaut
MA::MA (int nb_lettres, int nb_etats)
{
    // Initialisation de l'alphabet
    int i;
    for (i = 0; i < nb_lettres; ++i)
    {
        Sigma.insert ('a' + i);
    }

    // Initialisation de l'ensemble d'Ètats
    for (i = 0; i < nb_etats; ++i)
        Q.insert (i);

    // Initialisation de iota
    StateSet::iterator q;
    for (q = Q.begin (); q != Q.end (); ++q)
        iota[*q] = 0;
    if (nb_etats != 0)
        iota[*(Q.begin ())] = 1;

    // Initialisation de tau
    for (q = Q.begin (); q != Q.end (); ++q)
        tau[*q] = 1;

    // Initialisation de phi
    // On laise phi complËtement vide, l'automate n'a pas d'arÍte ‡ l'initialisation.
}

// Destructeur par dÈfaut
MA::~MA ()
{}


/// -------------------------- Affichage ----------------------------------

// Fonction d'affichage d'automate
string MA::message (void) const
{
    string
    mesg;	// La variable qui va contenir le message a afficher.
    char
    flottant[3 * (sizeof (float)) + 1];	// un tableau suffisamment grand pour Ècrire un flottant en dÈcimal.
    mesg = "--- AFFICHAGE de " + name + " ---\n";

    if (PFA_VERBOSE)
    {
        // Affichage du dictionnaire
        for (Dictionnaire::const_iterator d = alph.begin ();
                d != alph.end (); ++d)
        {
            mesg += "alph[";
            mesg += d->first;
            mesg += "]=" + d->second + " ; ";
        }
        mesg += "\n";
    }

    // Affichage de l'alphabet
    Alphabet::const_iterator l;
    mesg += "Sigma = {";
    if (Sigma.empty ())
    {
        mesg += "}\n";
    }
    else
    {
        for (l = Sigma.begin (); l != Sigma.end ();)
        {
            if (alph.find (*l) == alph.end ())
            {
                mesg += *l;
            }
            else
            {
                mesg += (alph.find (*l))->second;
            }

            if (++l != Sigma.end ())
            {
                mesg += ", ";
            }
            else
            {
                mesg += "}\n";
            }
        }
    }

    // Affichage des Ètats
    StateSet::const_iterator e;
    mesg += "Q = {";
    e = Q.begin ();
    if (e == Q.end ())
    {
        mesg += "}\n";
    }
    while (e != Q.end ())
    {
        sprintf (flottant, "%d", *e);
        mesg += flottant;
        e++;
        if (e != Q.end ())
        {
            mesg += ", ";
        }
        else
        {
            mesg += "}\n";
        }
    }

    // Affichage de iota
    SFunc::const_iterator s;
    stringstream convert;
    string nombre;
    mesg += "iota : ";
    for (s = iota.begin (); s != iota.end (); ++s)
    {
        mesg += "(";
        sprintf (flottant, "%d", s->first);
        mesg += flottant;
        mesg += ", ";
        //convert << s->second;
        //convert >> nombre;
        //mesg += nombre;
        sprintf (flottant, "%f", s->second);
        mesg += flottant;
        mesg += ") ";
    }
    mesg += "\n";

    // Affichage de tau
    mesg += "tau : ";
	int i=0;
    for (s = tau.begin (); s != tau.end (); ++s)
    {
        mesg += "(";
        sprintf (flottant, "%d", s->first);
        mesg += flottant;
        mesg += ", ";
        sprintf (flottant, "%f", s->second);
        mesg += flottant;
        mesg += ") ";
		if (i++%5 == 0) {
			mesg += "\n";
		}
    }
    mesg += "\n";

    // Affichage de phi
    mesg += "phi : ";
    if (phi.empty ())
        mesg += "vide";
    else
    {
        LettreSFunc::const_iterator a;
        Graph::const_iterator q;
        for (q = phi.begin (); q != phi.end (); q++)
        {
            for (a = q->second.begin ();
                    a != q->second.end (); a++)
            {
                for (s = a->second.begin ();
                        s != a->second.end (); s++)
                {
                    mesg += "(";
                    sprintf (flottant, "%d",
                             q->first);
                    mesg += flottant;
                    mesg += ",";

                    // on Ècrit la bonne lettre
                    if (alph.find (a->first) ==
                            alph.end ())
                        mesg += a->first;
                    else
                        mesg += alph.find (a->first)->second;

                    mesg += ",";
                    sprintf (flottant, "%d", s->first);
                    mesg += flottant;
                    mesg += ") = ";
                    sprintf (flottant, "%f", s->second);
                    mesg += flottant;
                    mesg += " ; ";
                }
            }
			mesg+="\n";
        }
    }
    mesg += "\n";

    return mesg;
}

// Chargement d'un fichier dot pour l'affichage de graphe
RESULT
MA::save_dot (const char *Filename,
               const int precision,
               const bool noms_etats,
               const bool multiples_trans) const
{
    try
    {
        if (Filename == NULL)
        {
            throw 0;
        }

        ofstream fp (Filename);
        if (!fp)
        {
            throw 1;
        }
        fp.precision (precision);	//fp.flags(ios::scientific);

        if (!name.empty ())
            fp << "digraph \"" << name << "\" {" << endl;
        else
            fp << "digraph G {" << endl;

        // On met le graphe en mode horizontal (Left to Right)
        fp << "graph [rankdir=LR]" << endl;

        Graph::const_iterator q;
        StateSet::const_iterator e;
        LettreSFunc::const_iterator a, b;
        SFunc::const_iterator s, r;
        StateSet W;
		Dictionnaire::const_iterator tl;

        for (e = Q.begin (); e != Q.end (); ++e)
        {
            // On crÈe la flÍche d'initialisation si il y a lieu
            if ((iota.find (*e) != iota.end ()) && (iota.find (*e)->second != 0))
            {
                fp << "iota" << *e << "[shape=plaintext, label = ";
                fp << (iota.find (*e))->
                second << " psfrag = \"$";
                fp << (iota.find (*e))->
                second << "$\", height=.15]" <<
                endl;
                fp << "iota" << *e << " -> " << *e <<
                "[minlen=1]" << endl;
            }

            // on cree l'etat avec ses attributs
            fp << *e;
            fp << " [ label = \"";
            if (noms_etats)
            {
                fp << *e << "\\n";
            }
            //      if ((iota.find(*e) != iota.end()) &&
            //          (iota.find(*e)->second != 0)) {
            //fp << "i=" << (iota.find(*e))->second << "\\n";
            //
            if ((tau.find (*e) != tau.end ())
                    && (tau.find (*e)->second != 0))
            {
                fp << (tau.find (*e))->second;
                fp << "\",peripheries=2";
            }
            else
            {
                fp << "\"";
            }
            fp << "]" << endl;

            // on ecrit ses successeurs
            q = phi.find (*e);
            if (q != phi.end ())
            {
                if (multiples_trans)
                {
                    for (a = q->second.begin ();a != q->second.end (); a++)
                    {
                        for (s =(a->second).begin ();s !=(a->second).end ();s++)
                        {
                            if (s->second !=0)
                                fp << "\t" << q->first;
								fp << " -> " << s->first;
								fp << "[label=\"";
								tl = alph.find(a->first);
								if (tl == alph.end()) {
									throw 3;
								}
								fp << tl->second;
								fp << "," << s->second;
								fp << "\"]" << endl;
                        }
                    }
                }
                else
                {	
					// le W sert a ne pas s'occuper plusieurs fois de l'etat d'arrive
                    W.clear();
                    for (a = q->second.begin (); a != q->second.end (); a++)
                    {
                        for (s =(a->second).begin ();s !=(a->second).end ();s++)
                        {
                            if (s->second !=0)
                            {
                                if (W.find(s->first)==W.end())
                                {
                                    W.insert (s->first);
                                    fp << "\t" << q->first << " -> " << s->first;
                                    fp << "[label=\"";
									tl=alph.find(a->first);
									if (tl == alph.end())
										throw 3;
									fp << tl->second << "," << s->second;
                                    for (b = a; ++b != q->second.end ();)
                                    {
                                        r = b->second.find (s->first);
                                        if (r != b->second.end ())
                                        {
                                            fp << "\\n";
											tl = alph.find (b->first);
											if (tl == alph.end())
												throw 3;
											fp << tl->second << "," << r->second;
                                        }
                                    }
                                    fp << "\"]" << endl;
                                }
                            }
                        }
                    }
                }
            }
        }

        fp << "}" << endl;
        return VAL (0);
    }
    catch (int erreur)
    {
        if (PFA_VERBOSE)
        {
            cerr << "ERREUR MA::save_dot" << endl;
            switch (erreur)
            {
            case 0:
                cerr << "\tFilename == NULL !!!" << endl;
                break;
            case 1:
                cerr << "\tErreur d'ouverture du fichier " << Filename << endl;
                break;
			case 3:
				cerr << "\tErreur du dictionnaire toutes les lettres ne sont pas definies" << endl;
            default:
                cerr << "Erreur n∞" << erreur << endl;
                break;
            }
        }
        return ERR (erreur);
    }
}

// Fonction d'affichage d'automate sur la sortie standart
RESULT MA::affiche (void) const
{
    cout << message ();
    return VAL (0);
}

/// --------------------- Entree sortie su fichier ------------------------------

RESULT MA::save (ofstream &fp, const FormatSauvegarde c_format) const 
{
    try
    {
        // les itÈrateur pour la structure
        Graph::const_iterator q;
        LettreSFunc::const_iterator a;
        SFunc::const_iterator s;
        StateSet::const_iterator e;
		Dictionnaire::const_iterator l;
	
		// la precision de sortie est maximale
        fp.precision (20);
        // first case format ffa or mdi
        switch (c_format)
        {
		case paffa:
            // calculation of the number of edges
            // Entete writing
            fp << "% MA file generated" << "\n";
			fp << "name \t" << name << "\n";
					
			fp << "alph\n";
			for (l=alph.begin() ; l != alph.end() ; ++l) {
				fp << l->first << " \t" << l->second << "\n";
			}
			fp << ".\n";
		
			fp << "iota\n";
			for (s=iota.begin() ;s != iota.end() ; ++s) {
                fp << s->first << " \t" << s->second << "\n";
			}
			fp << ".\n";
		
			fp << "tau\n";
			for (s=tau.begin() ; s != tau.end() ;++s) {
				fp << s->first << " \t" << s->second << "\n";
			}
			fp << ".\n";
		
			fp << "phi\n";
			for (q=phi.begin() ; q != phi.end() ; ++q) {
				for (a = q->second.begin() ; a != q->second.end() ; ++a) {
					for (s=a->second.begin() ; s != a->second.end() ; ++s) {
						fp << q->first << " \t";
						fp << a->first << " \t";
						fp << s->first << " \t";
						fp << s->second << "\n";
					}
				}
			}
			fp << ".\n";
		
            return 0;
            break;

        case paalergia:
            // =========== Format alergia !!! format avec perte d'information !!! ====
            // =========== il faut travailler avec des PDFAs !!!!!!!!!! ==============
            if (phi.find (0) == phi.end ())
                throw - 3;

            for (q = phi.begin (); q != phi.end (); ++q)
            {
                for (a = q->second.begin ();
                        a != q->second.end (); ++a)
                {
                    for (s = a->second.begin ();
                            s != a->second.end (); ++s)
                    {
                        fp << q->first << " ";
						fp << alph.find (a->first)->second;
						fp << " " << s->first;
						fp << " " << s->second << endl;
                    }
                }
            }
            return 0;
            break;
        default:
            throw - 4;
        }
    }
    catch (int e)
    {
        if (PFA_VERBOSE)
        {
            cerr << "ERROR in MA::Save !!!!" << endl;
            cerr << " \t";
            switch (e)
            {
            case -3:
                cerr << "There is no initial state" <<
                endl;
            case -4:
                cerr << "unrecognized format : " <<
                c_format << endl;
            default:
                cerr << "Erreur non prise en compte !" <<
                endl;
            }
        }
        return e;
    }
    catch (...)
    {
        if (PFA_VERBOSE)
        {
            cerr << "ERROR in MA::Save !!!" << endl;
            cerr << "Unexpected error *!#@???" << endl;
        }
        return ERR (0);
    }
}

// Sauvegarde le MA dans Filename
RESULT
MA::save (const char *Filename, const FormatSauvegarde c_format) const
{
	try {
		// On ouvre le fichier
		if (Filename == NULL)
			throw -1;
		ofstream fp (Filename);
		if (!fp.good ())
			throw -2;
	
		return MA::save(fp, c_format);
	}
	catch (int e) {
        if (PFA_VERBOSE)
        {
            cerr << "ERROR in MA::Save !!!!" << endl;
            cerr << " \t";
            switch (e)
            {
            case -1:
                cerr << "Filename == NULL !!!" << endl;
            case -2:
                cerr << "I can't open the file" <<
                Filename << endl;
            default:
                cerr << "Erreur non prise en compte !" <<
                endl;
            }
        }
        return e;
	}	
}

// erase cleanly all transitions with parameter between minval and maxval
RESULT MA::erase_transitions (const double maxval, const double minval)
{
    Graph::iterator q, qtmp;
    LettreSFunc::iterator a, atmp;
    SFunc::iterator s, stmp;

    // erase for iota
    for (s=iota.begin() ; s != iota.end() ; ) {
        if ((s->second >= minval) && (s->second <= maxval)) {
            stmp = s;
            ++stmp;
            iota.erase (s);
            s = stmp;
        }
        else {
            ++s;
		}
    }

    // erase for tau
    for (s=tau.begin() ; s != tau.end() ; ) {
        if ((s->second >= minval)&& (s->second <= maxval)) {
            stmp = s;
            ++stmp;
            tau.erase (s);
            s = stmp;
        }
        else {
            ++s;
		}
    }

    // erase for transitions
    for (q = phi.begin (); q != phi.end ();) {
        for (a = q->second.begin (); a != q->second.end ();) {
            for (s = a->second.begin (); s != a->second.end ();) {
                if ((s->second >= minval) && (s->second <= maxval)) {
                    stmp = s;
                    ++stmp;
                    a->second.erase (s);
                    s = stmp;
                }
                else {
                    ++s;
				}
            }
            if (a->second.empty ()) {
                atmp = a;
                ++atmp;
                q->second.erase (a);
                a = atmp;
            }
            else {
                ++a;
			}
        }
        if (q->second.empty ()) {
            qtmp = q;
            ++qtmp;
            phi.erase (q);
            q = qtmp;
        }
        else {
			++q;
		}
    }

    return 0;
}


// erase cleanly the transition phi(qbeg,l,qend)
RESULT MA::eraseTransition(const State qbeg, const Lettre l, const State qend) {
	Graph::iterator q;
	LettreSFunc::iterator a;
	SFunc::iterator s;

    q = phi.find(qbeg); 
	if (q != phi.end ())
    {
        a = q->second.find(l); 
		if (a != q->second.end ())
        {
            s = a->second.find(qend); 
			if (s != a->second.end ())
            {
                a->second.erase (s);                
            }
		
            if (a->second.empty ())
            {
                q->second.erase (a);
            }
        }
        if (q->second.empty ())
        {
            phi.erase (q);
        }
    }
	return VAL(0);
}


// emmonde l'automate
RESULT MA::emmonde(void) {
    // --- Explication de l'algorithme ------
    // On fait un parcour en profondeur en notant tous les Ètats atteint en partant
    // des sources (les Ètats initiaux)
    // Puis on inverse le graphe et on refait un parcours en profondeur en partant
    // des Ètats terminaux.
    // Les Ètat atteignables par les Ètats initiaux et qui ne sont pas atteint
    // par les Ètats terminaux du graphe inverse sont des Ètats ‡ supprimer.
    // --- Fin de l'explication de l'algorithme ---

    // d'abord le parcours inverse
    // 1 - On crÈe le graphe inverse
    Graph Inv; // The inverse graph
    Graph::iterator q; // some state
    LettreSFunc::iterator a; // some letter
    SFunc::iterator s, r;
    StateSet Qir; // set of states from which we can reach a terminal state
    StateSet Q2; // set of reachable states
    StateSet W;  // set for algorithm
    StateSet BadStates; // set of bad states
    State etat_cour; // current state

    // On inverse le graphe Inv(r,a,q) = phi (q,a,r)
    for (q = phi.begin (); q != phi.end (); ++q)
    {
        for (a = q->second.begin (); a != q->second.end (); ++a)
        {
            for (s = a->second.begin (); s != a->second.end (); ++s)
            {
                Inv[s->first][a->first][q->first] = s->second;
            }
        }
    }

    // Pour chaque Ètat terminal, on regarde l'ensemble des Ètats que
    // l'on atteint
    for (r = tau.begin (); r != tau.end (); ++r)
    {
        W.insert (r->first);
        Qir.insert (r->first);
    }

    while (!W.empty ())
    {
        etat_cour = *(W.begin ());
        W.erase (W.begin ());
        q = Inv.find (etat_cour);
        if (q != Inv.end ())
        {
            for (a = q->second.begin ();  a != q->second.end (); ++a)
            {
                for (s = a->second.begin (); s != a->second.end (); ++s)
                {
                    if (Qir.find (s->first) == Qir.end ())
                    {
                        Qir.insert (s->first);
                        W.insert (s->first);
                    }
                }
            }
        }
    }
    // Qir contient tous les Ètats atteignables dans le graphe inverse

    // Pour chaque Ètat initial on regarde si on tombe bien dans
    // l'ensemble des Ètats qui atteignent un Ètat terminal.
    W.clear() ; // inutile mais au cas o˘
    for (r = iota.begin (); r != iota.end (); ++r)
    {
        W.insert (r->first);
        Q2.insert (r->first);
    }

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
                    if (Q2.find (s->first) == Q2.end ())
                    {
                        Q2.insert (s->first);
                        W.insert (s->first);
                    }

                    if (Qir.find (s->first) == Qir.end ())
                    {
                        BadStates.insert (s->first);
                    }

                }
            }
        }
    }
    // Q2 contient l'ensemble des Ètats atteignables
    // BadStates contient l'ensemble des Ètats accessibles ‡ partir
    // desquels on ne peut atteindre un Ètat non terminal

    // Suppressions de tous les Ètats de BadStates
    StateSet::iterator bs;
    for (bs = BadStates.begin (); bs != BadStates.end (); ++bs)
    {
        erase (*bs);
    }

    // suppression des Ètats non accessibles
    BadStates=Q;
    for (bs = BadStates.begin() ; bs != BadStates.end() ; ++bs)
    {
        if (Q2.find(*bs) == Q2.end())
        {
            erase(*bs);
        }
    }
	return VAL(0);
}


// Charge le MA de Filename
RESULT MA::load (const char *Filename)
{
    try
    {
        // on modifie le ffa donc on perd l'attribution du ffa
        if (Filename == NULL)
            throw - 1;
        ifstream fp (Filename);
        if (!fp)
            throw - 2;
        fp.precision (30);
        string buf;
        char tmp[64];
        char c;
        bool finfichier;
        State q, r, s;
        Lettre a;
        Dictionnaire::iterator l;
        float proba;

        vide ();	// on vide le MA
        fp >> buf;
        if (buf == "%Alergia")
        {		// ================= cas format alergia =============
            while (fp.good ())
            {	// jusqu'‡ la fin du fichier
                fp >> q;	// on lit l'Ètat de dÈpart
                addState (q);
                fp >> buf;	// on lit la lettre de la transition
                a = addLettre (buf);
                fp >> r;	// on lit l'Ètat d'arrivÈ
                addState (r);
                fp >> proba;	// on lit la probabilitÈ de la transition
                addTransition (q, a, r, proba);
            }

            // L'Ètat 0 est le seul Ètat initial
            iota[0] = 1;

            // On met ‡ jour tau
            StateSet::const_iterator iq;
            Graph::const_iterator ie;
            LettreSFunc::const_iterator ia;
            SFunc::const_iterator ir;

            for (iq = Q.begin (); iq != Q.end (); ++iq)
            {
                proba = 0;	// we use proba to make the sum
                ie = phi.find (*iq);
                if (ie != phi.end ())
                {
                    for (ia = ie->second.begin ();
                            ia != ie->second.end (); ++ia)
                    {
                        for (ir = ia->second.begin ();
                                ir != ia->second.end ();
                                ++ir)
                        {
                            proba += ir->second;
                        }
                    }
                }
                tau[*iq] = 1 - proba;
            }
        }
		else if (buf == "%MDI") 
		{ // ========================== cas format MDI ==================
			float f1, f2;
			int numstates;
			char cbuf[20];
			while (fp.good()) {
				if (*(buf.begin()) == '%') {	// =========================== commentaire =====================
                    finfichier = ((fp >> c) == NULL);
                    while ((!finfichier) && (c != '\n'))
                    {
                        finfichier = ((fp.get (c)) == NULL);
                    }
                    finfichier = ((fp >> buf) == NULL);
                }
				else if (buf == "Name") {
					getline(fp,name);
					fp >> buf;
				}
				else if (buf == "NumStates") {
					fp >> numstates;
					for (int i=0 ; i< numstates ; ++i) {
						addState(i);
					}
					fp >> buf;
				}
				else if (buf == "NumEdges") {
					fp >> buf;
					fp >> buf;
				}
				else if (buf == "State") {					
					fp >> q; // on lit l'etat courant					
					fp.getline(cbuf,20,'=');
					buf = cbuf;
					if (buf == " initial") {
						fp >> iota[q];
						fp.getline(cbuf,20,'=');
						buf = cbuf;
					}
					if (buf == " final") {
						fp >> f1;
						fp.get(c);
						fp >> f2;
						tau[q]=f1/f2;
					}
					else {
						cerr << "PARSE ERROR final : q=" << q << endl;
					}
					fp >> buf; // on recupere le %Ci=0.xxxxxxx
					fp >> buf; // on recupere le numero de l'etat
					// on commence a lire les transition sortantes
					while ((buf != "State") && (fp.good())) {
						fp >> s; // on lit l'etat d'arrive
						fp >> buf; // on lit la lettre
						// on efface les guillemets
						buf.erase(buf.begin());
						buf.erase(--buf.end());
						fp.getline(cbuf,20,'=');
						fp >> f1;
						fp.get(c);
						fp >> f2;
						if (buf.empty()) {
							tau[q]+=f1/f2;
						}
						else {
							a=addLettre(buf);
							addTransition(q,a,s,f1/f2);
							//phi[q][a][s]=f1/f2;
						}
						fp >> buf; // on lit le %Ci=0.xxxxxxx
						fp >> buf;
					}
				}
			}
		}
		else if (buf == "%OLDMA")
        {		// ================= cas format ffa interne =============
            while (fp.good ())
            {	// Tand que l'on est pas ‡ la fin du fichier
                if (*(buf.begin()) == '%')
                {	// =========================== commentaire =====================
                    finfichier = ((fp >> c) == NULL);
                    while ((!finfichier) && (c != '\n'))
                    {
                        finfichier =
                            ((fp.get (c)) ==
                             NULL);
                    }
                    finfichier = ((fp >> buf) == NULL);
                }
                // ============================== State ========================
                else if (buf == "State")
                {
                    // --------------- Lecture de la ligne initiale ----------------
                    fp >> q;	// on lit l'Ètat de dÈpart
                    Q.insert (q);	// on ajoute l'etat q
                    sprintf (tmp, "%d", q);
                    finfichier = ((fp >> buf) == NULL);	// on lit la suite
                    while ((buf != "State")
                            && (buf != tmp)
                            && (!finfichier))
                    {
                        if (buf == "c")
                        {
                            fp >> buf;	// on lit le =
                            fp >> buf;	// on lit le chiffre mais on en fait rien
                        }
                        else if (buf == "i")
                        {	// on lit la probabilitÈ d'Ítre initial
                            fp >> buf;	// on lit le =
                            fp >> proba;	// on lit la probabilitÈ
                            if (proba != 0)
                                iota[q] =
                                    proba;
                        }
                        else if (buf == "f")
                        {	// on lit la probabilitÈ d'Ítre terminal
                            fp >> buf;	// on lit le =
                            fp >> proba;	// on lit la probabilitÈ
                            if (proba != 1)
                                tau[q] = 1 -
                                         proba;
                        }
                        finfichier = ((fp >> buf) == NULL);	// on lit la suite
                    }
                    // ---------------- Lecture des transitions sortantes ----------------
                    r = q;
                    while ((buf != "State")
                            && (!finfichier))
                    {
                        if (PFA_VERBOSE)
                        {
                            if (buf != tmp)
                            {
                                cerr << "load WARNING !" << endl;
                                cerr << "La premiËre lettre d'une colone ne correspond pas ‡ l'Ètat" << endl;
                            }
                        }
                        fp >> s;	// on lit l'Ètat d'arrivÈ
                        fp >> buf;	// on lit la lettre

						a = addLettre(buf);
                        //~ // on cherche la lettre correspondante
                        //~ for (l = alph.begin ();
                                //~ (l != alph.end ())
                                //~ && (l->second != buf);
                                //~ l++)
                            //~ ;
                        //~ if (l == alph.end ())
                        //~ {	// si on ne retrouve pas la lettre
                            //~ a = addLettre ();	// on ajoute une lettre
                            //~ alph[a] = buf;	// ‡ laquelle on fait correspondre la bonne chaine de caractËre
                        //~ }
                        //~ else	// sinon
                            //~ a = l->first;	// on attribu la lettre en cours.

                        fp >> buf;	// on lit 'p'
                        fp >> buf;	// on lit '='
                        fp >> proba;	// on lit la proba de l'arÍte

                        phi[q][a][s] = proba;	// on insËre l'arÍte
                        finfichier = ((fp >> buf) == NULL);	// on lit la lettre suivante.
                    }
                }
                else
                {	// Autres cas, on ne fait rien d'autre qu'avancer
                    fp >> buf;
                }
            }
		        {		// ================= cas format ffa interne =============
            while (fp.good ())
            {	// Tand que l'on est pas ‡ la fin du fichier
                if (*(buf.begin()) == '%')
                {	// =========================== commentaire =====================
                    finfichier = ((fp >> c) == NULL);
                    while ((!finfichier) && (c != '\n'))
                    {
                        finfichier =
                            ((fp.get (c)) ==
                             NULL);
                    }
                    finfichier = ((fp >> buf) == NULL);
                }
                // ============================== State ========================
                else if (buf == "State")
                {
                    // --------------- Lecture de la ligne initiale ----------------
                    fp >> q;	// on lit l'Ètat de dÈpart
                    Q.insert (q);	// on ajoute l'etat q
                    sprintf (tmp, "%d", q);
                    finfichier = ((fp >> buf) == NULL);	// on lit la suite
                    while ((buf != "State")
                            && (buf != tmp)
                            && (!finfichier))
                    {
                        if (buf == "c")
                        {
                            fp >> buf;	// on lit le =
                            fp >> buf;	// on lit le chiffre mais on en fait rien
                        }
                        else if (buf == "i")
                        {	// on lit la probabilitÈ d'Ítre initial
                            fp >> buf;	// on lit le =
                            fp >> proba;	// on lit la probabilitÈ
                            if (proba != 0)
                                iota[q] =
                                    proba;
                        }
                        else if (buf == "f")
                        {	// on lit la probabilitÈ d'Ítre terminal
                            fp >> buf;	// on lit le =
                            fp >> proba;	// on lit la probabilitÈ
                            if (proba != 1)
                                tau[q] = 1 -
                                         proba;
                        }
                        finfichier = ((fp >> buf) == NULL);	// on lit la suite
                    }
                    // ---------------- Lecture des transitions sortantes ----------------
                    r = q;
                    while ((buf != "State")
                            && (!finfichier))
                    {
                        if (PFA_VERBOSE)
                        {
                            if (buf != tmp)
                            {
                                cerr << "load WARNING !" << endl;
                                cerr << "La premiËre lettre d'une colone ne correspond pas ‡ l'Ètat" << endl;
                            }
                        }
                        fp >> s;	// on lit l'Ètat d'arrivÈ
                        fp >> buf;	// on lit la lettre

						a = addLettre(buf);
                        //~ // on cherche la lettre correspondante
                        //~ for (l = alph.begin ();
                                //~ (l != alph.end ())
                                //~ && (l->second != buf);
                                //~ l++)
                            //~ ;
                        //~ if (l == alph.end ())
                        //~ {	// si on ne retrouve pas la lettre
                            //~ a = addLettre ();	// on ajoute une lettre
                            //~ alph[a] = buf;	// ‡ laquelle on fait correspondre la bonne chaine de caractËre
                        //~ }
                        //~ else	// sinon
                            //~ a = l->first;	// on attribu la lettre en cours.

                        fp >> buf;	// on lit 'p'
                        fp >> buf;	// on lit '='
                        fp >> proba;	// on lit la proba de l'arÍte

                        phi[q][a][s] = proba;	// on insËre l'arÍte
                        finfichier = ((fp >> buf) == NULL);	// on lit la lettre suivante.
                    }
                }
                else
                {	// Autres cas, on ne fait rien d'autre qu'avancer
                    fp >> buf;
                }
            }

        }
        return 0;
    }



		else {
			stringstream str;
			Lettre a;
			TrueLettre b;
            while (fp.good ())
            {	// Tand que l'on est pas ‡ la fin du fichier			
                if (*(buf.begin()) == '%')
                {	// =========================== commentaire =====================
					fp.get(c);
                    while ((fp.good()) && (c != '\n')) {
                        fp.get (c);
                    }
                    fp >> buf;
                }
				else if (buf == "name") {
					fp.get(c);
					name+=c;
                    while ((fp.good()) && (c != '\n')) {
                        fp.get (c);
						name+=c;
                    }
                    fp >> buf;
				}
                else if (buf == "alph")
                {	// dictionnaire
					fp >> buf;
					while ((buf != ".") && (fp.good())) {						
				        if (*(buf.begin()) == '%')
						{	// --- commentaire ---
							fp.get(c);
							while ((fp.good()) && (c != '\n')) {
								fp.get (c);
							}
							fp >> buf;
						}
						str << buf;
						str >> a;
						fp >> b;
						Sigma.insert(a);
						alph[a]=b;
						fp >> buf;
					}
					fp >> buf;
                }
                else if (buf == "iota")
                {	// iota
					fp >> buf;
					while ((buf != ".") && (fp.good())) {
					
						if (*(buf.begin()) == '%')
						{	// --- commentaire ---
							fp.get(c);
							while ((fp.good()) && (c != '\n')) {
								fp.get (c);
							}
							fp >> buf;
						}
					
						q=atoi(buf.c_str());
						fp >> proba;
						Q.insert(q);
						iota[q]=proba;
						fp >> buf;
					}
					fp >> buf;
                }
			    else if (buf == "tau")
                {	// iota
					fp >> buf;
					while ((buf != ".") && (fp.good())) {

						if (*(buf.begin()) == '%')
						{	// --- commentaire ---
							fp.get(c);
							while ((fp.good()) && (c != '\n')) {
								fp.get (c);
							}
							fp >> buf;
						}

						q=atoi(buf.c_str());
						fp >> proba;
						Q.insert(q);
						tau[q]=proba;
						fp >> buf;
					}
					fp >> buf;
                }
				else if (buf == "phi") {
					fp >> buf;
					while ((buf != ".") && (fp.good())) {
					
						if (*(buf.begin()) == '%')
						{	// --- commentaire ---
							fp.get(c);
							while ((fp.good()) && (c != '\n')) {
								fp.get (c);
							}
							fp >> buf;
						}
					
						q=atoi(buf.c_str());
						fp >> a;						
						fp >> s;						
						fp >> proba;						
					
						if (alph.find(a) == alph.end()) {
							cerr << "!!! Lettre : " << a << endl;
							throw -3;
						}
					
						Q.insert(q);
						Q.insert(s);
						phi[q][a][s]=proba;
					
						fp >> buf;
					}
					fp >> buf;
				}
                else
                {	// Autres cas, on ne fait rien d'autre qu'avancer
                    fp >> buf;
                }
            }
		}
		return VAL(0);
	}
	catch (int erreur)
    {
        if (PFA_VERBOSE)
        {
            cerr << "MA::load" << endl;
            switch (erreur)
            {
            case -1:
                cerr << "ERROR in MA::Save !!!!" << endl;
                cerr << " \tFilename == NULL !!!" << endl;
                break;
            case -2:
                cerr << "Erreur lors de l'ouverture du fichier : " << Filename << endl;
                break;
			case -3:
				cerr << "Erreur alph n'est pas bien defini, des lettres manques !" << endl;
            default:
                cerr << "Erreur inconnue ?" << endl;
            }
        }
        return erreur;
    }
    catch (...)
    {
        if (PFA_VERBOSE)
        {
            cout << "MA::load(" << Filename << ")" << endl;
            cout << "Erreur inexplicable ?" << endl;
        }
        return ERR (0);
    }
}


/// ----------------------- Modifieurs -------------------------------

// Fonction qui ajoute une transition
RESULT MA::addTransition (const Transition & t, const float val)
{
    if (PFA_SAFE)
    {
        assert (Q.find (t.qdep) != Q.end ());
        assert (Sigma.find (t.a) != Sigma.end ());
        assert (Q.find (t.qarr) != Q.end ());
    }

	if (val != 0) 
	{
		phi[t.qdep][t.a][t.qarr] = val;
	}
	else {
		Graph::const_iterator q;
		q = phi.find(t.qdep);
		if ( q != phi.end()) {
			LettreSFunc::const_iterator a = q->second.find(t.a);
			if (a != q->second.end()) {
				SFunc::const_iterator s = a->second.find(t.qarr);
				if (s != a->second.end()) {
					phi[t.qdep][t.a][t.qarr] = val;
				}
			}		
		}
	}

	if (PFA_SAFE)
    {
        assert (coherent ());
    }
    return VAL (0);
}


// Ajout une lettre par defaut ‡ l'alphabet
Lettre MA::addLettre (void)
{
    Lettre
    a;
    if (!Sigma.empty ())
    {
        a = *(--(Sigma.end ())) + 1;
    }
    else
    {
        a = 'a';
    }

    Sigma.insert (a);
    alph[a] = a;
    return a;
}

// Ajout de la lettre buf
Lettre MA::addLettre (const TrueLettre & buf)
{
    // search whether buf is a new letter
    Lettre a;
    Dictionnaire::iterator l;
    for (l = alph.begin (); (l != alph.end ()) && (l->second != buf); ++l)
        ;
    if (l == alph.end ())
    {			// if buf is new then
        a = addLettre ();	// we add the letter
        alph[a] = buf;	// and we associate it with buf
    }
    else {
        a = l->first;
	}
    return a;
}


// Ajoute un Ètat
State MA::addNewState (const float init, const float term)
{
    State q;
    if (!Q.empty ())
    {
        q = *(--Q.end()) + 1;
        while (Q.find(q) != Q.end())
            ++q;
    }
    else
    {
        q = 1;
    }
    addState(q,init,term);
    return q;
}


// Ajoute l'Ètat q s'il n'existe pas
State MA::addState (const State q, const float init, const float term)
{
    if (Q.find (q) == Q.end ())
    {
        Q.insert (q);
    }
    if (init != 0)
        iota[q]=init;
    if (term != 0)
        tau[q]=term;
    return q;
    if (PFA_SAFE)
        assert(coherent());
}


// Ajout d'une transition
RESULT MA::addTransition (const State qdep, const Lettre a, const State qarr, const float val)
{
    if (PFA_SAFE)
    {
        assert (Q.find (qdep) != Q.end ());
        assert (Sigma.find (a) != Sigma.end ());
        assert (Q.find (qarr) != Q.end ());
    }
    if (val != 0)
    {
        phi[qdep][a][qarr] = val;
    }

    if (PFA_SAFE)
    {
        assert (coherent ());
    }
    return VAL (0);
}

// Cree un MA aleatoire
RESULT MA::becomeRandom (const int nb_etats, const int nb_lettres,
                          const int num_graphe, const float densite,
                          const float prob_init,    const float prob_term,
                          const float min_trans,  const float max_trans)
{
    float init, term;
    int i;
    StateSet::const_iterator q, s;
    Alphabet::const_iterator a;

    // on vide le MA
    vide ();
    // on initialise la variable aleatoire
    if (num_graphe != 0)
        srand (num_graphe);

    // on ajoute les états
    for (i = 0; i < nb_etats; i++)
    {
        if (random () < prob_init)
            init = random (min_trans, max_trans);
        else
            init = 0;
		
        if (random () < prob_term)
            term = random (min_trans, max_trans);
        else
            term = 0;
		
        addNewState (init, term);
    }

    // on ajoute les lettres
    for (i = 0; i < nb_lettres; i++)
    {
        addLettre ();
    }

    // on cree les transitions avec des distributions uniformes
    for (q = Q.begin (); q != Q.end (); q++)
    {			// pour tous les Ètats
        for (a = Sigma.begin (); a != Sigma.end (); a++)
        {
            for (s = Q.begin (); s != Q.end (); s++)
            {
                if (random () < densite)
                {	// si on tire un nombre inferieur ‡ la densitÈ, on cree une arete
                    phi[*q][*a][*s] =
                        random (min_trans, max_trans);
                }
            }
        }
    }
    return VAL (0);
}

// Cree un MA aleatoire
RESULT MA::becomeRandomMax (const int nb_etats, const int nb_lettres,
							const int num_graphe, const int nb_succ,
							const int nb_init,    const int nb_term,
							const float min_trans,  const float max_trans)
{
    int i;
    StateSet::const_iterator q, s;
    Alphabet::const_iterator a;
	
    // on vide le MA
    vide ();
    // on initialise la variable aleatoire
    if (num_graphe != 0)
        srand (num_graphe);
	
    // on ajoute les Ètats
    for (i = 0; i < nb_etats; ++i)
    {
        addState(i,0,0);
    }
	
    // on ajoute les lettres
    for (i = 0; i < nb_lettres; i++)
    {
        addLettre ();
    }
	
    for (i=0 ; i< nb_init ; ++i)
    {
        iota[rand()%nb_etats]=random(min_trans,max_trans);
    }
	
    for (i=0 ; i<nb_term ; ++i)
    {
        tau[rand() % nb_etats]= random(min_trans,max_trans);
    }
	
    // on cree les transitions avec des distributions uniformes
    for (q = Q.begin (); q != Q.end (); ++q)
    {			// pour tous les Ètats
        for (a = Sigma.begin (); a != Sigma.end (); ++a)
        {
            for (i=0 ; i<nb_succ ; ++i)
            {
                phi[*q][*a][rand()%nb_etats] = random (min_trans, max_trans);
            }
        }
    }
    return VAL (0);
}


RESULT MA::renormalize_state(const State q, const float normalization, const float min_tau, const float max_tau, const float min, const float max) {
	const int ERR_INPUT=-1023;
	if (min_tau > max_tau) {
		cerr << "MA::renormalize_state : ";
		cerr << "min_tau > max_tau" << endl;
		throw ERR_INPUT;
	}
	if (min > max) {
		cerr << "MA::renormalize_state : ";
		cerr << "min > max" << endl;
		throw ERR_INPUT;
	}
	
	Graph::iterator r;
	LettreSFunc::iterator a;
	SFunc::iterator s;
	float sum, x;
	bool taumodif=true;
	sum=tau[q];
	r=phi.find(q); 
	if (r != phi.end()) {
		for (a=r->second.begin() ; a != r->second.end() ; a++) {
			for (s=a->second.begin() ; s!=a->second.end() ; s++) {
				sum+=s->second;
			}
		}
	}
	bool good=false;
	StateSet BadStates;
	while ( !good ) {
		good=true;
		x=normalization/sum;
		if (tau[q]*x > max_tau) {
			tau[q]=max_tau;
			if (taumodif) {
				taumodif=false;
				good=false;
			}
		} else if (tau[q]*x < min_tau) {
			tau[q]=min_tau;
			if (taumodif) {
				taumodif=false;
				good=false;
			}
		} else {
			tau[q]=x*tau[q];
		}
		sum=tau[q];
		if (r != phi.end()) {
			for (a=r->second.begin() ; a != r->second.end() ; a++) {
				for (s=a->second.begin() ; s!=a->second.end() ; s++) {
					if (x*s->second > max) {
						s->second = max;
						if (BadStates.find(s->first) == BadStates.end()) {
							BadStates.insert(s->first);
							good=false;
						}
					} 
					else if (x*s->second < min) {
						s->second=min;
						if (BadStates.find(s->first) == BadStates.end()) {
							BadStates.insert(s->first);
							good=false;
						}
					}
					else {
						s->second=x*s->second;
					}
				}
			}
		}
	}
	return VAL(0);
}

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
RESULT MA::becomeRandomControl (const int nb_etats, 
								const int nb_lettres, 
								const int num_graphe, 
								const int max_nb_succ, 
								const int max_nb_init, 
								const int max_nb_term, 
								const float min_trans,  
								const float max_trans, 
								const float min_iota, 
								const float max_iota, 
								const float min_tau, 
								const float max_tau, 
								const float normalization, 
								const float prob_init, 
								const float prob_trans, 
								const float prob_term)
{
	
    int i;
    StateSet::const_iterator q, s;
	SFunc::iterator r;
    Alphabet::const_iterator a;
	
	const int ERR_INPUT=-1023;
	if (min_iota > max_iota) {
		cerr << "min_iota > max_iota !!!" << endl;
		throw ERR_INPUT;
	}
	if (min_trans > max_trans) {
		cerr << "min_trans > max_trans !!!" << endl;
		throw ERR_INPUT;
	}
	if (min_tau > max_tau) {
		cerr << "min_tau > max_tau !!!" << endl;
		throw ERR_INPUT;
	}
	
    // on vide le MA
    vide ();
    // on initialise la variable aleatoire
    if (num_graphe != 0)
        srand (num_graphe);

    // on ajoute les Ètats
    for (i = 0; i < nb_etats; i++)
    {
        addState(i,0,0);
    }

    // on ajoute les lettres
    for (i = 0; i < nb_lettres; i++) {
        addLettre ();
    }
	
	/// --- iota ---
	i=0;
	for (i=0 ; i<max_nb_init ; i++) {
		if (float(rand())/RAND_MAX < prob_init) {
			iota[rand()%nb_etats]=random(min_iota, max_iota);
		}
	}
	// en mode normalization
	if (normalization != 0) {
		/// renormalization
		// calculus of the sum of iota
		float sum;
		for (sum=0, r=iota.begin() ; r != iota.end() ; r++) {
			sum+=r->second;
		}
		float x;
		bool good=false;
		StateSet BadStates;
		while ( !good ) {
			x=normalization/sum;
			good=true;
			for (r=iota.begin() ; r != iota.end() ; r++) {
				if (x * r->second > max_iota) { 
					sum+=r->second=max_iota;
					if (BadStates.find(r->first) == BadStates.end()) {
						BadStates.insert(r->first);
						good=false;
					}
				} 
				else if (x * r->second < min_iota) {
					sum+=r->second=min_iota;
					if (BadStates.find(r->first) == BadStates.end()) {
						BadStates.insert(r->first);
						good=false;
					}
				}
				else {
					sum+=r->second = x*r->second;
				}
			}
		}
	}
	
	/// --- Tau + Transitions ---
	i=0;
	for (i=0 ; i<max_nb_term ; i++) {
		if (float(rand())/RAND_MAX < prob_term) {
			tau[rand()%nb_etats]=random(min_tau, max_tau);
		}
	}

    // on cree les transitions avec des distributions uniformes
    for (q = Q.begin (); q != Q.end (); ++q)
    {			// pour tous les Ètats
        for (a = Sigma.begin (); a != Sigma.end (); ++a)
        {
            for (i=0 ; i<max_nb_succ ; ++i)
            {
				if (float(rand())/RAND_MAX < prob_trans) {
					phi[*q][*a][rand()%nb_etats] = random (min_trans, max_trans);
				}
            }
        }
		if (normalization != 0) {
			renormalize_state(*q,normalization,min_tau, max_tau, min_trans,max_trans);
		}
    }
    return VAL (0);
}


RESULT MA::becomeRandomPrefix (const int nb_etats,
                                const int nb_lettres,
                                const int num_graphe,
                                const float densite,
                                const float prob_init,
                                const float prob_term,
                                const float min_trans,
                                const float max_trans)
{
    float
    term;
    int
    i;
    StateSet
    R,
    S;		// l'ensemble des Ètats restants
    StateSet::const_iterator q, s;
    Alphabet::const_iterator a;

    // on vide le MA
    vide ();
    // on initialise la variable aleatoire
    if (num_graphe != 0)
        srand (num_graphe);

    // on ajoute les Ètats
    for (i = 0; i < nb_etats; i++)
    {
        if (random () < prob_term)
            term = random (min_trans, max_trans);
        else
            term = 0;
        addState (0, term);
    }
    // on ajoute l'unique Ètat terminal
    iota[*Q.begin ()] = random (min_trans, max_trans);

    R = Q;
    S.clear ();		// on initialise R et S

    // on ajoute les lettres
    for (i = 0; i < nb_lettres; i++)
    {
        addLettre ();
    }

    // on supprime l'Ètat correspondant ‡ epsilon de R et on l'ajoute ‡ S
    R.erase (R.begin ());
    S.insert (*Q.begin ());
    // on cree les transitions avec des distributions uniformes
    for (q = Q.begin (); q != Q.end (); q++)
    {			// pour tous les Ètats
        for (a = Sigma.begin (); a != Sigma.end (); a++)
        {
            if ((random () < densite) || (R.empty ()))
            {	// on rabat les arÍtes sur les Ètats dÈja affectÈs
                for (s = S.begin (); s != S.end (); s++)
                {
                    if (random () < densite)
                    {	// si on tire un nombre inferieur ‡ la densitÈ, on cree une arete
                        phi[*q][*a][*s] =
                            random (min_trans,
                                    max_trans);
                    }
                }
            }
            else
            {
                s = R.begin ();
                R.erase (s);
                S.insert (*s);
                phi[*q][*a][*s] =
                    random (min_trans, max_trans);
            }
        }
    }

    return VAL (0);
}

/// --------------------- Methodes ------------------------------
void MA::normaliseVecteur(SFunc &V) const {
    SFunc::iterator q;
    double sum;
    sum = 0;
    for (q=V.begin() ; q != V.end() ; q++) {
        sum += q->second>0?q->second:-q->second;
    }
    sum = 1/sum;
    for (q=V.begin() ; q != V.end() ; q++) {
        q->second *= sum;
        if (q->second == 0) {
            V.erase(q);
        }
    }
    // // Affichage de debogage
    // cout << "\tMA::normaliseVecteur\n";
    // for (q=V.begin() ; q != V.end() ; q++) {
    //     cout << "\tV[" << q->first << "]=" << q->second << endl;
    // }
    // cout << "\t---\n";
}


// genere un mot en fonction du modèle. En particulier, il faut récupérer
// la valeur des P_q(\Se) pour tous les états q.
// Description de l'algorithme :
//   Soit r un langage rationnel, on construit P_r de la façon suivante ;
//  Neg(u)=\sum_{x \in \Sigma}r(ux\Se)<0?r(ux\Se):0 + r(u)<0?r(u):0
//  P_r(ux\Se) = r(ux\Se)/Neg(u) si r(a\Se)>0 et 0 sinon
//  P_r(\eps) = r(\eps)/Neg(u)

RESULT MA::genere_mot(Word &w, const bool safe){
    // Dans un premier temps, on rÈcupËre les valeurs des 
    // P_q(\Se)
	SFunc PSe; // contiendra la valeur des P_q(\Se)

	if (val_PSe(PSe) != VAL(0)) {
		if ( safe ) {
			return ERR(0);
		}
	}
	
    // on a w donnÈ, on calcule P(w), P(wx\Se) pour tout x et on conserve
    // les w^-1P(\eps) et w^-1P(x\Se) positifs
	w="";    // le mot courant initialisÈ ‡ epsilon
    Lettre c;     // la lettre courante
    bool termine=false;
    double Negw, Posw;
    Alphabet::const_iterator x;
    SFunc V, Vtmp;
    SFunc::const_iterator q;
    LettreSFunc::const_iterator a;
    SFunc::const_iterator s;
	Graph::const_iterator r;
    map<Lettre,double> px; // va contenir px[0] = P(w), px[1] = P(wx_1\Se), ... , px[n] = P(wx_n\Se)
    map<Lettre,double>::const_iterator l;
    bool lettreTrouve;
	double val;
	double sum;

    V = iota;

    while (!termine) {
        // normaliser le vecteur V, fait en sorte de mettre
        // la somme de ses valeurs ‡ 1 et supprimes les entrÈes
        // nulles (qui sont inutiles)
        normaliseVecteur(V);
        Negw=0; Posw=0;
        // calcul de w^-1P(\eps)
        px[0] = 0;
        for (q = V.begin() ; q != V.end() ; q++) {
            px[0]+= q->second * val_tau(q->first);
        }
        // on met ‡ jour Neg(w) et Pos(w)
        if (px[0] < 0) 
            Negw += px[0];
        else
            Posw += px[0];

        // // AFFICHAGE DEBUG
        // cout << "px[.]:" << endl;
        // for (l=px.begin() ; l != px.end() ; l++) {
        //     cout << "px[";
        //     cout << l->first==0?'0':l->first;
        //     cout << "]=";
        //     cout << l->second << endl;
        // }
        // cout << "Negw=" << Negw << ", Posw=" << Posw << endl;
        // cout << " -- px[.]" << endl;
        // // END AFFICHAGE DEBUG

        // initialisation de px
        for (x=Sigma.begin() ; x != Sigma.end() ; x++)
            px[*x]=0;
        // calcul de w^-1P(x\Se)
        for (q = V.begin() ; q != V.end() ; q++) {
            r=phi.find(q->first);
            if (r != phi.end()) {
                for (a=r->second.begin() ; a != r->second.end() ; a++) {
                    for (s=a->second.begin() ; s != a->second.end() ; s++) {
                        // // AFFICHAGE DEBUG
                        // cout << "\t\t**********" << endl;
                        // cout << "q->first=" << q->first << endl;
                        // cout << "q->second=" << q->second << endl;
                        // cout << "r->first=" << r->first << endl;
                        // cout << "a->first=" << a->first << endl;
                        // cout << "s->first=" << s->first << endl;
                        // cout << "s->second=" << s->second << endl;
                        // cout << "PSe[s->first]=" << PSe[s->first] << endl;
                        // cout << "\t\t**********" << endl;
                        // // FIN AFFICHAGE DEBUG
                        px[a->first] += q->second * s->second * PSe[s->first];
                    }
                }
            }
        }
        for (x=Sigma.begin() ; x != Sigma.end() ; x++) {
            // mise ‡ jour de Neg(w) et Pos(w)
            if (px[*x] < 0) 
                Negw += px[*x];
            else
                Posw += px[*x];
        }

        // // AFFICHAGE DEBUG
        // cout << "px[.]: " << endl;
        // for (l=px.begin() ; l != px.end() ; l++) {
        //     cout << "px[";
        //     cout << l->first==0?'0':l->first ;
        //     cout << "]=";
        //     cout << l->second ;
        //     cout << endl;
        // }
        // cout << "Negw=" << Negw << ", Posw=" << Posw << endl;
        // cout << " -- px[.]" << endl;
        // // END AFFICHAGE DEBUG

        // on a fini de calculer les valeurs possibles. On choisit seulement
        // parmis les valeurs positives
        if (Posw == 0) {
            cerr << "ERREUR !!! MA::genere_mot" << endl;
			cerr << "prefixe = ";
			cerr << affiche(w) << endl;
            cerr << "Impossible d'utiliser cette methode avec ce prefixe dans ce MA" << endl;
            cerr << "Pos(\\epsilon)=0 !!!" << endl;
            throw ERR(21);
        }
        val = random(0,Posw);
        sum = 0;
        lettreTrouve=false;
        for (l=px.begin() ; l != px.end() && !lettreTrouve ; l++) {
            // // AFFICHAGE DE DEBUG
            // cout << "\tval=" << val;
            // cout << ", sum=" << sum;
            // cout << ", Posw=" << Posw;
            // cout << ", Negw=" << Negw;
            // cout << "\tpx[";
            // cout << l->first==0?'0':l->first;
            // cout << "]="; 
            // cout << l->second << endl;
            // // AFFICHAGE DE DEBUG
            if (l->second > 0) {
                sum+= l->second;
                if (sum >= val) {
                    if (l->first != 0) {
                        c = l->first;
                        w += c;
                        // // AFFICHAGE DEBUG
                        // cout << "c=";
                        // cout << c;
                        // cout << ", w=" ;
                        // cout << affiche(w);
                        // cout << endl;
                        // // AFFICHAGE DEBUG
                    }
                    else {
                        termine=true;
                        // // AFFICHAGE DEBUG
                        // cout << "termine" << endl ;
                        // // AFFICHAGE DEBUG
                    }
                    lettreTrouve=true;
                }
            }
        }

        if (!termine) {
            Vtmp.clear();
            // Mise ‡ jour du vecteur d'initialisation
            for (q=V.begin() ; q != V.end() ; q++) {
                r=phi.find(q->first);
                if (r != phi.end()) {
                    a = r->second.find(c);
                    if (a != r->second.end()) {
                        for (s=a->second.begin() ; s != a->second.end() ; s++) {
                            Vtmp[s->first] += q->second * s->second;
                        }
                    }
                }
            }
            V = Vtmp;
        }
    }
    return VAL(0);
}

// genere un Èchantillon de mots
RESULT MA::genere_echantillon ( int taille, 
                                Sample & S,
                                const int num_echantillon)
{
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

// genere un Èchantillon de mots : l'echantillon reste identique si lot est identique
RESULT MA::genere_echantillon (const int taille,
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


// Calcule les valeurs de P_q(\S_e) pour
// tous les Ètats q
RESULT MA::val_PSe(SFunc &PSe) const {
	Simplex simplex;
	StateSet::const_iterator q,r;
	Alphabet::const_iterator a;
	int j;
	float tmp;

	simplex.set_number_of_unknown_to(Q.size());
	
	for (q=Q.begin() ; q != Q.end() ; q++) {
		simplex.setval(-val_tau(*q));
		for (r=Q.begin(), j=1 ; r != Q.end() ; r++, j++) {
			if (*r == *q) {
				tmp=-1;
			}
			else { 
				tmp=0;
			}
			for (a=Sigma.begin() ; a != Sigma.end() ; a++) {
				tmp += val_phi(*q,*a,*r); 
			}
			simplex.setparam(j,tmp);
		}
		simplex.add_equality_constraint();
	}
	
	simplex.set_mode(realNotBounded);     // on utilise des réels pouvant êtres négatifs
	
	float epsilon;
	if (! simplex.has_solution(PSe,epsilon) ) {
		PSe.clear();
		simplex.affiche();
		cerr << "ERREUR INEXPLICABLE !!!!\n";
		cerr << "Le système d'équation doit être mal posé !\n";
		throw ERR(0);
	} else {
		return VAL(0);
	}
}

/// --------------------- Testeurs ------------------------------

// Teste la cohÈrence de la structure de donnÈe.
bool MA::coherent (void) const
{
    bool
    res = true;
    // on cherche une transition contenant un Ètat n'appartenant pas
    // ‡ l'ensemble des Ètats.
    Graph::const_iterator q;
    LettreSFunc::const_iterator a;
    SFunc::const_iterator s;
    for (q = phi.begin (); q != phi.end (); ++q)
    {
        for (a = q->second.begin (); a != q->second.end (); ++a)
        {
            for (s = a->second.begin (); s != a->second.end (); ++s)
            {
                if ( (Q.find(q->first) == Q.end ())
                        || (Q.find (s->first) == Q.end ()))
                {
                    res = false;
                    if (PFA_VERBOSE)
                    {
                        cerr << "Erreur de coherence (phi) : (";
                        cerr << q->first << "," << a->first << ",";
                        cerr << s->first << ") = " << s->second << " !!! " << endl;
                    }
                }
            }
        }
    }
    return res;
}

// Teste si le MA est coherent et est un SPFA. valeur de iota, tau et phi dans [0,1]
bool MA::isSPFA (void) const
{

    if (!coherent ())
    {			// si le ffa n'est pas coherent
        return false;
    }
    // si le ffa est coherent
    // Test de iota
    SFunc::const_iterator s;
    float sum = 0;
    // cerr.precision(30);
    for (s = iota.begin (); s != iota.end (); ++s)
    {
        if ((s->second < 0) || (s->second > 1))
        {
            //~ cerr << "iota " << s->first << "<0 ou >1 !!! : " << s->second << endl;
            return false;
        }
        sum += s->second;
    }
    if (!simeqone (sum))
    {
        return false;
    }
    //cerr << "iota OK" << endl;

    // Test de tau
    for (s = tau.begin (); s != tau.end (); ++s)
    {
        if ((s->second < 0) || (s->second > 1))
        {
            //~ cerr << "tau " << s->first << "<0 ou >1 : " << s->second << endl;
            return false;
        }
    }

    // Test of  phi and
    // of the property : \forall q \in Q and \forall a \in \Sigma
    // tau(q) + \sum_{s\in Q} phi(q,a,s)=1
    Graph::const_iterator q;
    LettreSFunc::const_iterator a;
    for (q = phi.begin (); q != phi.end (); ++q)
    {
        sum = val_tau (q->first);
        for (a = q->second.begin (); a != q->second.end (); ++a)
        {
            for (s = a->second.begin (); s != a->second.end (); ++s)
            {
                if ((s->second < 0) || (s->second > 1))
                {
                    //~ cerr << "phi(" << q->first << ",";
                    //~ cerr << a->first << ",";
                    //~ cerr << s->first << ")=";
                    //~ cerr << s->second << " FALSE !!!" << endl;
                    return false;
                }
                //cerr << "phi(" << q->first << "," << a->first << "," << s->first << ")=" << s->second << endl;
                //cerr << "sum = " << sum << endl;
                sum += s->second;
            }
        }
        if (!simeqone (sum))
        {
            //cerr << "phi(" << q->first << ",*,*)=" << sum << " FALSE !!!" << endl;
            return false;
        }
    }
    // on a passÈ toutes les Èpreuves
    return true;
}

// Fonction qui rend le MA vide (pas d'Ètats pas de lettre)
void
MA::vide ()
{
    name.clear ();
    Q.clear ();
    Sigma.clear ();
    iota.clear ();
    tau.clear ();
    phi.clear ();
    alph.clear ();
}

// renvoi vrai si le MA est un PFA
bool MA::isPFA () const
{
    if (!isSPFA ())
        return false;

    // --- Explication de l'algorithme ------
    // On fait un parcour en profondeur en notant tous les Ètats atteint en partant
    // des sources (les Ètats initiaux)
    // Puis on inverse le graphe et on refait un parcours en profondeur en partant
    // des Ètats terminaux.
    // Si l'ensemble des Ètats atteint par les Ètats terminaux dans le graphe
    // inverse contient l'ensemble des Ètats atteint en partant des Ètats initiaux
    // alors le SPFA  considÈrÈ est un PFA.
    // --- Fin de l'explication de l'algorithme ---

    // d'abord le parcours inverse
    // 1 - On crÈe le graphe inverse
    Graph
    Inv;
    Graph::const_iterator q;
    LettreSFunc::const_iterator a;
    SFunc::const_iterator s, r;
    StateSet
    Qir,
    Q2,
    W;
    State
    etat_cour;

    // On inverse le graphe Inv(r,a,q) = phi (q,a,r)

    for (q = phi.begin (); q != phi.end (); q++)
    {
        for (a = (q->second).begin (); a != (q->second).end (); a++)
        {
            for (s = (a->second).begin ();
                    s != (a->second).end (); s++)
            {
                Inv[s->first][a->first][q->first] = s->second;
            }
        }
    }

    // Pour chaque Ètat terminal, on regarde l'ensemble des Ètats que
    // l'on atteint
    for (r = tau.begin (); r != tau.end (); r++)
    {
        W.insert (r->first);
        Qir.insert (r->first);
    }

    while (!W.empty ())
    {
        etat_cour = *(W.begin ());
        W.erase (W.begin ());
        q = Inv.find (etat_cour);
        if (q != Inv.end ())
        {
            for (a = (q->second).begin ();
                    a != (q->second).end (); a++)
            {
                for (s = (a->second).begin ();
                        s != (a->second).end (); s++)
                {
                    if (Qir.find (s->first) == Qir.end ())
                    {
                        Qir.insert (s->first);
                        W.insert (s->first);
                    }
                }
            }
        }
    }

    // Pour chaque Ètat initial on regarde si on tombe bien dans
    // l'ensemble des Ètats qui atteignent un Ètat terminal.
    for (r = iota.begin (); r != iota.end (); ++r)
    {
        W.insert (r->first);
        Q2.insert (r->first);
    }

    while (!W.empty ())
    {
        etat_cour = *W.begin ();
        W.erase (W.begin ());
        q = phi.find (etat_cour);
        if (q != phi.end ())
        {
            for (a = (q->second).begin ();
                    a != (q->second).end (); ++a)
            {
                for (s = (a->second).begin ();
                        s != (a->second).end (); ++s)
                {
                    if (Q2.find (a->first) == Q2.end ())
                    {
                        Q2.insert (a->first);
                        W.insert (a->first);
                    }
                    if (Qir.find (s->first) == Qir.end ())
                    {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}


string MA::affiche (const Word & w) const
{
    Word::const_iterator a;
    string res;
	Dictionnaire::const_iterator x;
    for (a = w.begin() ; a != w.end () ; ++a)
    {
		x=alph.find(*a);
		if (x == alph.end()) {
			res = "Dictionnaire ERROR ";
			res += *a;
			return res;
		}
        res += x->second;
        res += ' ';
    }
    return res;
}

void MA::next (Word & w) const
{
    Word::reverse_iterator ll;	// last letter
    Word::iterator z;
    Lettre a;
    Lettre first = *Sigma.begin ();
    int i;
    int level = w.size ();
    bool end = false;
    while ((w.size () > 0) && (!end))
    {
        ll = w.rbegin ();
        a = *ll;
        ++a;
        if (Sigma.find (a) != Sigma.end ())
        {
            ++(*ll);
            end = true;
        }
        else
        {
            z=w.end();
            --z;
            w.erase(z); // w.pop_back ();
        }
    }
    if (!end)
        ++level;
    for (i = w.size (); i < level; ++i)
    {
        w += first; // w.push_back (first);
    }

}


// erase the state badstate and all transitions from and to it
RESULT MA::erase (State badstate)
{
    Graph::iterator q, qtmp;
    LettreSFunc::iterator a, atmp;
    SFunc::iterator s, stmp;
    StateSet::iterator r;

    // Suppression dans l'ensemble d'Ètats
    r=Q.find (badstate);
    if ( r == Q.end () )
    {
        cerr << "MA::erase (State badstate) : demande de suppression d'un etat inexistant" << endl;
        exit (-1);
    }
    Q.erase (r);
    // Suppression dans iota
    s = iota.find(badstate);
    if (s != iota.end ())
        iota.erase (s);
    // Suppression dans tau
    s = tau.find(badstate);
    if (s != tau.end ())
        tau.erase (s);
    // Suppression dans phi
    if (phi.find (badstate) != phi.end ())
        phi.erase (phi.find (badstate));

    for (q = phi.begin (); q != phi.end ();)
    {
        for (a = q->second.begin (); a != q->second.end ();)
        {
            for (s = a->second.begin (); s != a->second.end ();)
            {
                if (s->first == badstate)
                {
                    stmp = s;
                    ++stmp;
                    a->second.erase (s);
                    s = stmp;
                }
                else
                    ++s;
            }
            if (a->second.empty ())
            {
                atmp = a;
                ++atmp;
                q->second.erase (a);
                a = atmp;
            }
            else
                ++a;
        }
        if (q->second.empty ())
        {
            qtmp = q;
            ++qtmp;
            phi.erase (q);
            q = qtmp;
        }
        else
            ++q;
    }

    return VAL (0);
}

StateSet MA::delta (const State q, const Word &w) const
{
    StateSet R;
    R.insert(q);
    return delta(R,w);
}

StateSet MA::delta (const StateSet &R, const Word &w) const
{
    StateSet V; // le vecteur correspondant aux Ètats atteints
    StateSet X; // un vecteur temporaire
    StateSet::iterator xi;
    Graph::const_iterator q;
    LettreSFunc::const_iterator b;
    SFunc::const_iterator r;
    Word::const_iterator a;

    //~ cout << endl;
    V=R;
	for (a=w.begin() ;a != w.end() ; a++)    
    {   
        X.clear();
        for (xi=V.begin() ; xi != V.end() ; ++xi)
        {
            q=phi.find(*xi);
            if (q != phi.end())
            {
                b=q->second.find(*a);
                if (b != q->second.end())
                {
                    for (r = b->second.begin() ; r  != b->second.end() ; ++r)
                    {
                        X.insert(r->first);
                    }
                }
            }
        }
        V=X;
    }
    return V;
}


void MA::allTransitions(TransitionFunction &T) const {
	T.clear();
	Graph::const_iterator q;
	LettreSFunc::const_iterator a;
	SFunc::const_iterator r;
	Transition t;
	for (q=phi.begin() ; q != phi.end() ; ++q) {
		t.qdep= q->first;
		for (a=q->second.begin() ; a!=q->second.end() ; ++a) {
			t.a=a->first;
			for (r=a->second.begin() ; r != a->second.end() ; ++r) {
				t.qarr=r->first;
				T[t]=0;
			}
		}
	}
}

void MA::allStates(SFunc &S) const {
	StateSet::const_iterator q;
	S.clear();
	for (q=Q.begin() ; q!=Q.end() ;++q) {
		S[*q]=0;
	}
}

// Become the MA with maximum number of transitions
RESULT MA::becomeComplete (const int nb_etats, const int nb_lettres) {
	clear();	
	float valinit = 1/ float(nb_etats);
	float val = 1/float(nb_etats + 1);
	int i;

	for (i=0 ; i<nb_lettres ; i++) {
		addLettre();
	}

	for (i=0 ; i<nb_etats ; i++) {
		addNewState(valinit, val);
	}

	StateSet::const_iterator q1,q2;
	Alphabet::const_iterator a;

	for (q1=Q.begin() ; q1 != Q.end() ; q1++) {
		for (a=Sigma.begin() ; a != Sigma.end() ; a++) {
			for (q2=Q.begin() ; q2 != Q.end() ; q2++) {
				addTransition(*q1,*a,*q2,val);
			}
		}
	}
	return VAL(0);
}

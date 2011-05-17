/***************************************************************************
                          main.C  -  description
                             -------------------
    begin                : mer oct  9 19:14:09 CEST 2002
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

using namespace std;
#include "interface.H"
#include "test.H"
#include "main.H"

// #define MAXFLOAT 3.40282347e+38F

/// ########################################################################## ///
// --- Les variables globales d'entrées qui parametre l'expérimentation

// Les fichier

string temporary_directory= "/tmp/dees/"; // le répertoire temporaire

string ficffaA   = temporary_directory + "A";	// Fichier du ffa A
string ficffaB   = temporary_directory + "B";	// Fichier du ffa B
string ficdotA   = temporary_directory + "A.dot";	// Le fichier dot associé à A
string ficdotB   = temporary_directory + "B.dot";	// Le fichier dot associé à B
string ficech    = temporary_directory + "A.ech";	// Le fichier échantillon généré par A
string ficechpre = temporary_directory + "A.ecp";	// Le fichier échantillon préfixe généré par A
string ficpsA    = temporary_directory + "A.ps";	// La fichier de visualisation de A
string ficpsB    = temporary_directory + "B.ps";	// Le fichier de visualisation de B

// Les Èchantillons
Sample S;			// main sample
Sample Stest;		// test sample
int num_ech = 0;		// Le numero de l'Èchantillon
int num_ech_test = 1;		// Le numbero de l'Èchantillon test
int taille_ech = 1000;		// La taille de l'Èchantillon
int taille_ech_test = 1000;	// La taille de l'Èchantillon test

// Le PFA alÈatoire
int numPFA = 0;			// Le numero du PFA gÈnÈrÈ
int nb_etats = 5;		// Le nombre d'Ètats du PFA gÈnÈrÈ
int nb_lettres = 2;		// Le nombre de lettres du PFA gÈnÈrÈ

// Les prÈcisions
double precision = 0.5;		// la prÈcision (0 -> prefix tree, 1 -> single state)
double precpow = -1;		// La prÈcision en mode exponentiel
double epsprime = 0;		// La valeur de suppression des transitions (0 delete nothing, 1 delete all)
unsigned int seuil = 0;		// le seuil en dessous duquel un prefixe n'est plus pris en compte dans l'apprentissage
double seuilbm = .0001 ;    // Le seuil pour le test de convergence de BM
unsigned int maxmots = INT_MAX; // le nombre maximal de mots à ajouter à l'ensemble de test à chaque étapes

// L'affichage
string nom_A = "A";		// le nom du PFA A
string nom_B = "B";		// le nom du PFA B
bool affiche_A = true;		// affiche l'automate source ou non
bool affiche_B = true;		// affiche l'automate destination
bool verbose = false;		// mode verbose pour le test
bool affiche_etats = false;  // affiche les Ètats dans le .dot
int  viewprec=3;			// precision des parametres a l'affichage
bool viewstates=false;		// on met le numero des etats
bool viewmultipletransitions=false; // vrai pour ne plus regrouper les transitions
bool interactive=false;     // interactive mode

// L'apprentissage
T_ModeReturn moderet = end;	// Le mode de retour des arÍtes (end ou begin)
T_ModeEpsilon modeeps = variable;	// Le mode d'apprentissage avec epsilon fixe ou variable (epsfixed ou variable)
float test_proportion = 0.2;	// la proportion de l'echantillon attribuee au test
int nb_tours = 100;	// Le nombre d'itÈrations pour Baum Welch
int nb_essais = 10; // le nombre d'essais pour le best learning
int max_states = INT_MAX; // le nombre maximal d'etats pour DEES
bool stepssave = false; // vrai si on veut enregistrer toutes les étapes de l'apprentissage.

// L'exportation de fichier
string format = "ffa";

// les distances
bool smooth = false;

string cmd;
PPRFA A;
PPRFA B;

stringstream convert;		// variable qui permet de convertir les stream en entiers ou flottants

/// ########################################################################## ///


// affiche un mot
string
affiche (const Word & w, const Dictionnaire & alph)
{
    Word::const_iterator a;
    string res;
    for (a = w.begin (); a != w.end (); ++a)
    {
        res += alph.find (*a)->second;
        res += ' ';
    }
    return res;
}


void show(MA A) {	
	string dotfile;
	string pdffile;
	stringstream trad;
	trad << temporary_directory << "pfa"<< int(time(NULL));
	trad >> pdffile;
	dotfile=pdffile;
	pdffile += ".pdf";
	dotfile += ".dot";
	A.save_dot(dotfile, viewprec, viewstates, viewmultipletransitions);
	string cmd;
    cmd += "dot -Tepdf ";
    cmd += dotfile;
	cmd += ">";
	cmd += pdffile;
	system(cmd.c_str());
    show_ps(pdffile);
}



// La procÈdure qui initialise le random
int
randomize (void)
{
    // mise a jour du random
    unsigned int i = 0;
    char c;
    ifstream dev_rnd ("/dev/urandom");
    if (dev_rnd != NULL)
    {
        for (int k = 0; k < 4; k++)
        {
            dev_rnd.get (c);
            i = i << 8 | c;
        }
        srand (i);
        return i;
    }
    else
    {
        srand (num_ech);
        return num_ech;
    }
}

int
affiche_utilisation (void)
{
    string fichelp = "./help/";
    fichelp += "help.txt";
    ifstream fh (fichelp.c_str ());
    if (fh == NULL)
    {
        cout << "DEES : ";
        cout << "WARNING: " << fichelp << " not found." << endl;
    }
    else
    {
        cout << fh.rdbuf ();
        fh.close ();
        cout << endl;
    }
    return 0;
}

int
affiche_options (void)
{
    string fichelp = "/Users/esposito/dees/help/";
    fichelp += "/help_options.txt";
    ifstream fh (fichelp.c_str ());
    if (fh == NULL)
    {
        cout << "DEES : ";
        cout << "WARNING: " << fichelp << " not found." << endl;
    }
    else
    {
        cout << fh.rdbuf ();
        fh.close ();
        cout << endl;
    }
    return 0;
}


RESULT
export_to_pdf (const MA & A,
              const string & fic_pdf,
              const string & fic_tmp,
              const int precision,
              const bool noms_etats,
              const bool multiples_trans)
{
    A.save_dot (fic_tmp, precision, noms_etats, multiples_trans);	// On l'enregistre dans un fichier dot
    // on Ècrit la commande
    string cmd;
    cmd = "dot -Tepdf ";
    cmd += fic_tmp;
    cmd += " > ";
	cmd += fic_pdf;
	
    // on execute la transformation
    system (cmd.c_str ());
    return 0;
}

void
show_ps (string & fic_pdf)
{
    string cmd;
    cmd += "open ";
    cmd += fic_pdf;
    cmd += "&";
    system (cmd.c_str ());
}


// Affiche toutes les distances entre A et B
void
affiche_distance (PFA & A, const PFA & B, const int taille_ech_test,
                  const string & flag)
{
    double pA, pB;
    A.genere_echantillon (taille_ech_test, S, num_ech_test);
	cerr << "Test_Sample_size\t";
	cerr << "precision\t";
	cerr << "perplex A\t";
	cerr << "perplex B\t";
	cerr << "KL div\t";
	cerr << "d_2\t";
	cerr << "L_1\t";
	cerr << "d_2(nlog)\t";
	cerr << "L_1(nlog)\t";
    cout << taille_ech_test << "\t";
    cout << precision << "\t";
    cout << (pA = A.perplexite (S)) << "\t";
    cout << (pB = B.perplexite (S)) << "\t";
    cout << pB - pA << "\t";
    cout << A.d_2 (B, S) << "\t";
    cout << A.L_1 (B, S) << "\t";
    cout << A.d_2nlog (B, S) << "\t";
    cout << A.L_1nlog (B, S) << "\t";
    cout << flag << endl;
}

void
affiche_distance (PFA & A, const PFA & B, const int taille_ech,
                  const char *flag)
{
    string s_flag = flag;
    affiche_distance (A, B, taille_ech, s_flag);
}

void
affiche_classe (PFA & A)
{
    if (!A.coherent ())
    {
        cout << "uncoherent MA" << endl;
    }
    else if (!A.isSPFA ())
    {
        cout << "MA" << endl;
    }
    else if (!A.isPFA ())
    {
        cout << "SPFA" << endl;
    }
    else if (!A.isPRFA ())
    {
        cout << "PFA" << endl;
    }
    else if (!A.isPDFA ())
    {
        cout << "PRFA" << endl;
    }
    else
    {
        cout << "PDFA" << endl;
    }
}


// Crossed validation learning
// Apprentissage avec validation croisÈe

void Alergia(PFA &A, float prec, string ficsample, bool verbose) {
	stringstream AL;
	string ficsampletmp=temporary_directory + "tmpsample.alergia";
	Sample S;
	S.load(ficsample);
	S.save(ficsampletmp, alergia);
	AL << "alergia";
	AL << " -a " << prec;
	AL << " -o " << ficffaA;
	AL << " -f " << ficsampletmp;
	if (verbose) {
		cout << ":: " << AL.str() << endl;
	}
	system(AL.str().c_str());

	stringstream MOD;
	string fictmp=temporary_directory + "alergia.tmp";
	MOD << "echo %Alergia > " << fictmp;
	MOD << " && cat " << ficffaA << " >> " << fictmp;
	if (verbose) {
		cout << "transformation de l'automate avec :\n";
		cout << MOD.str() << endl;
	}
	system(MOD.str().c_str()); // on ajoute %Alergia au dÈbut du fichier

	if (verbose) {
		cout << "on charge le fichier : " << fictmp << endl;
	}

	A.load(fictmp);
}

void MDI(PFA &A, float prec, string ficsample, bool verbose) {
	string ficsampletmp=temporary_directory + "tmpsample.mdi";
	Sample S;
	S.load(ficsample);
	S.save(ficsampletmp, mdi);
	stringstream MDI;
	MDI << "mdi";
	MDI << " -b";
	if (verbose) {
		MDI << " -v";
	}
	MDI << " -c " << prec / double(taille_ech);
	MDI << " -o " << ficffaA;
	MDI << " -i " << ficsampletmp;
	if (verbose) {
		cout << ":: " << MDI.str() << endl;
	}
	system(MDI.str().c_str());	

	stringstream MOD;
	string fictmp=temporary_directory + "mdi.tmp";
	MOD << "echo %MDI > " << fictmp;
	MOD << " && cat " << ficffaA << " >> " << fictmp;

	if (verbose) {
		cout << "transformation de l'automate avec :\n";
		cout << MOD.str() << endl;
	}
	system(MOD.str().c_str()); // on ajoute %MDI au dÈbut du fichier

	if (verbose) {
		cout << "on charge le fichier : " << fictmp << endl;
	}

	A.load(fictmp);
}

// --- La procedure pour debbugger ---
float fonction (float val, float moy, float ecart)
{
    float constante = 1/(ecart * sqrt(2 * 3.1415));
    return constante * exp(- ((val - moy)*(val - moy)) / ecart);
}

void BaumWelch(PFA A, Sample S, int nb_turns, string output) {
	TransitionFunction T;
	SFunc Iota, Tau;
	double val, oldval;
	A.allTransitions(T);
	A.allStates(Iota);
	A.allStates(Tau);
	oldval=A.Likelihood(S);
	for (int i=0 ; i < nb_turns ; i++) {
		A.BaumWelch(S,T,Iota,Tau,3,true);
		val=A.Likelihood(S);
		//if (abs(val-oldval)<.1) {
		//	break;
		//}
		oldval=val;
		if (verbose) {
			cout << "Likelihood: " << A.Likelihood(S) << "\n";
		}
	}
	cout << endl;
	
	A.erase_bad_states();
	A.save(output);
}


// génère un PA complet

RESULT generate_complete_random_pfa(int argc, char **argv, int i) {
	
	if (argc < i + 2) {
		cerr << "dees --grp nbstates nbletters file" << endl;
		cerr << "generate a complete random automaton with 'nbstates' states\n";
		cerr << "'nbletters' letters and save it in file\n";
		cerr << "usefull option:\n";
		cerr << "--numPFA : the number of the ma by default random number is choosen.\n";
		return -1;
	}
	
	int nbstates=atoi(argv[i++]);
	int nbletters=atoi(argv[i++]);
	string output=argv[i++];
	
	int num_graphe=0;      //cout << "num graphe: ";    cin >> num_graphe;
	float min_trans=0.01;     //cout << "min_trans: ";     cin >> min_trans;
	float max_trans=1;     //cout << "max_trans: ";     cin >> max_trans;
	float prob_init=1;     //cout << "prob_init: ";     cin >> prob_init;
	float prob_trans=1;    //cout << "prob_trans: ";    cin >> prob_trans;
	float prob_term=1;     //cout << "prob_term: ";     cin >> prob_term;

	if (verbose) {
		cout << "num graphe: "    << num_graphe    << endl;
		cout << "prob_trans: "   << prob_trans    << endl;
		cout << "prob_init: "    << prob_init     << endl;
		cout << "prob_term: "    << prob_term     << endl;
		cout << "min_trans: "     << min_trans     << endl;
		cout << "max_trans: "     << max_trans     << endl;
	}				
	
	A.becomeRandom (nbstates, nbletters, numPFA, prob_trans, prob_init, prob_term, min_trans, max_trans);
				
	return A.save (output);
}

RESULT generate_random_pfa(int argc, char **argv, int i) {
	if (argc < i+3) {
		cerr << "dees --random nbstates nbletters file\n";
		cerr << "generate a random automaton of class 'class' with 'nbstates' states\n";
		cerr << "'nbletters' letters and save it in file\n";
		cerr << "usefull option:\n";
		cerr << "--numPFA : the number of the ma by default random number is choosen.\n";
		cerr << "--interactive : ask you to enter values\n";
		return -1;
	}
	int nbstates=atoi(argv[i++]);
	int nbletters=atoi(argv[i++]);
	string output=argv[i++];
	
	int num_graphe=0;      
	int max_nb_succ=nbstates;
	int max_nb_init=nbstates;
	int max_nb_term=nbstates;
	float min_trans=0;     
	float max_trans=1;     
	float min_iota=0;      
	float max_iota=1;      
	float min_tau=0;       
	float max_tau=1;       
	float normalization=1; 
	float prob_init=.3;    
	float prob_trans=.3;   
	float prob_term=.3;    
	
	if (interactive) {
		cout << "num graphe (0=random): ";    cin >> num_graphe;
		cout << "max_nb_succ (width of the graph): ";   cin >> max_nb_succ;
		cout << "max_nb_init: ";   cin >> max_nb_init;
		cout << "max_nb_term: ";   cin >> max_nb_term;
		cout << "min_trans (minimal value of a transition): ";     cin >> min_trans;
		cout << "max_trans (maximal value of a transition): ";     cin >> max_trans;
		cout << "min_iota  (minimal initial value): ";      cin >> min_iota;
		cout << "max_iota  (maximal initial value): ";      cin >> max_iota;
		cout << "min_tau   (minimal terminaison value): ";       cin >> min_tau;
		cout << "max_tau   (maximal terminaison value): ";       cin >> max_tau;
		cout << "normalization (value of normalisation mostly 1): "; cin >> normalization;
		cout << "prob_init  (probability for a state to be initial: ";     cin >> prob_init;
		cout << "prob_trans (probability for a terminaison to araise): ";    cin >> prob_trans;
		cout << "prob_term  (probability for a state to be terminal): ";     cin >> prob_term;
	}
	
	A.becomeRandomControl (nbstates, nbletters, num_graphe, max_nb_succ, max_nb_init, max_nb_term, min_trans, max_trans, min_iota, max_iota, min_tau, max_tau, normalization, prob_init, prob_trans, prob_term);
	
	return A.save(output);
}

RESULT show_sample_forward(int argc, char**argv, int i) {
	if (argc < i+2) {
		cerr << "usage : --forward PFA sample\n";
		cerr << "return a list of couple word and its probability for the PFA\n";
		throw 2;
	}
	A.load(argv[i++]);
	S.load(argv[i++]);
	if (A.dictionnaire() != S.dictionnaire()) {
		cerr << "Alphabets of the PFA and of the sample differ !" << endl;
		Dictionnaire D;
		Dictionnaire::const_iterator x;
		cerr << "dictionnaire of the PFA\n";
		D=A.dictionnaire();
		for (x=D.begin() ; x != D.end() ; x++) {
			cout << x->first << ":" << x->second << "\n";
		}
		
		cerr << "dictionnaire of the Sample\n";
		D=S.dictionnaire();
		for (x=D.begin() ; x != D.end() ; x++) {
			cout << x->first << ":" << x->second << "\n";
		}
		throw 2;
	}
	Sample::const_iterator w;
	double somme=0 ;
	double logsomme = 0;
	double val ;
	for (w=S.begin() ; w != S.end() ; w++) {
		cout << A.affiche(w->first) <<  " : ";
		cout.flush();
		val = A.p_bar(w->first);
		somme += val;
		cout << val ;
		val = A.plog_bar(w->first);
		logsomme += val;
		cout << ", log forward : " << val;
        cout << ", exp(log forward) :" << exp(val);
		cout << ", pBarDirect : " << A.p_bar_directe(w->first);
        cout << endl;
	}
	cout << "sum = " << somme << ", somme de exp(log(forward)) : " << logsomme << endl;
	return VAL(0);
}


RESULT show_sample_proba(int argc, char **argv, int i){
	if (argc < i+2) {
		cerr << "usage : --proba PFA sample\n";
		cerr << "return a list of couple word and its probability for the PFA\n";
		throw 2;
	}
	A.load(argv[i++]);
	S.load(argv[i++]);
	if (A.dictionnaire() != S.dictionnaire()) {
		cerr << "Alphabets of the PFA and of the sample differ !" << endl;
		Dictionnaire D;
		Dictionnaire::const_iterator x;
		cerr << "dictionnaire of the PFA\n";
		D=A.dictionnaire();
		for (x=D.begin() ; x != D.end() ; x++) {
			cout << x->first << ":" << x->second << "\n";
		}
		
		cerr << "dictionnaire of the Sample\n";
		D=S.dictionnaire();
		for (x=D.begin() ; x != D.end() ; x++) {
			cout << x->first << ":" << x->second << "\n";
		}
		throw 2;
	}
	Sample::const_iterator w;
	double somme=0 ;
	double logsomme = 0;
	double val ;
	for (w=S.begin() ; w != S.end() ; w++) {
		cout << A.affiche(w->first) <<  " : ";
		cout.flush();
		val = A.p_directe(w->first);
		somme += val;
		cout << val ;
		val = A.plog(w->first);
		logsomme += val;
		cout << ", log p : " << val;
        cout << ", exp(log p) :" << exp(val);
        cout << endl;
	}
	cout << "sum = " << somme << ", somme de exp(log(p)) : " << logsomme << endl;
	return VAL(0);
}

void affiche_PSe(void) {
	SFunc PSe;
	SFunc::const_iterator q;
	int res;
		
	res = A.val_PSe(PSe);
	cout << "res = " << res << endl;
	
	for (q=PSe.begin() ; q != PSe.end() ; q++) {
		cout << "P_" << q->first << "(\Se) = " << q-> second << endl;
	}
}


/// ############################################################################# ///
/// ###                           INIT COMMAND LINE + MAIN                   #### ///
/// ############################################################################# ///

// initialisation des arguments (ligne de commande)
int
initialise_arguments (int argc, char *argv[])
{
    string arg;
    int i = 1;
    while ((i < argc) && (argv[i][0] == '-'))
    { // pour chaque argument commenÁant par -
        arg = argv[i];
        if ((arg == "-h") || (arg == "--help"))
        {		// affiche l'aide
            throw 1;
        }
		// options muettes
        else if (
				 (arg == "--export") || // exporte vers des formats d'automates differents
				 (arg == "--convert") || // convertit des echantillons
				 (arg == "--dist") || // affiche la distance entre deux ffa
				 (arg == "--class") || // affiche la classe d'un pfa
				 (arg == "--sample") || // genere un Èchantillon
				 (arg == "--grp") ||  // genere un automate probabiliste
				 (arg == "--gcrp") || // genere un automate probabiliste complet
				 (arg == "--mdi") || // algorithme mdi
				 (arg == "--alergia") || // algorithme alergia
				 (arg == "--bm") || // algorithme Baum Welch 
				 (arg == "--affiche") || // affiche un automate sur la sortie standard
				 (arg == "-H") || // affiche l'aide Ètendue (options)
				 (arg == "-Y" || arg == "-I") || // envoie la fonction Y (mode interactif)
				 (arg == "--test") || // fait un test ‡ partir de l'automate
				 (arg == "--random") || // generate a random automaton
				 (arg == "--proba") || // show probabilities of words of a sample
				 (arg == "--forward") || // show forward values of words of a sample
				 (arg == "--deletetransitions") || // delete transitions of some value
				 ((arg == "-P") || (arg == "--showps")) || // affichage postscript
				 (arg == "--PSe")       // Affichage des valeurs de P(\Se) pour chaque état
				 
				 )
		{
		}
		else if (
				 (arg == "--dees")    || // algorithme dees (cible PRFA)
				 (arg == "--deesha")  || // algorithme deesha (cible MA)
				 (arg == "--deesdet") )  // algorithme basé sur dees mais se limitant aux solution déterministes
		{
			precision=1;
		} else if (arg == "--deesbm") // algorithme dees Baum Welch
		{
			precision=2; 
		}
        else if ((arg == "--DEBUG")
                 || (arg == "-d") || (arg == "--debug"))
        {
			// PFA_DEBUG=true;
        }
		else if (arg == "--interactive") 
		{
			interactive=true; 
			++i;
		}
        else if (arg == "--format")
        {
            format = argv[++i];
        }
        else if ((arg == "--VERBOSE") || (arg == "-v") || (arg == "--verbose"))
        {
            verbose = true;
            //      PFA_VERBOSE=true;
        }
        else if ((arg == "--MUTE") || (arg == "-m") || (arg == "--mute"))
        {
            //      PFA_VERBOSE=false;
            verbose = false;
        }
        else if ((arg == "--SAFE") || (arg == "--safe"))
        {
            //      PFA_SAFE=true;
        }
        else if ((arg == "--UNSAFE") || (arg == "--unsafe"))
        {
            //      PFA_SAFE=false;
        }
        else if ((arg == "--taille_ech") || (arg == "-te"))
        {
            convert << argv[++i];
            convert >> taille_ech;            
        }
		else if (arg == "--maxmots") {
			convert << argv[++i];
			convert >> maxmots;
		}
        else if (arg == "--num_ech")
        {
            num_ech = atoi (argv[++i]);
        }
        else if (arg == "--num_ech_test")
        {
            num_ech_test = atoi (argv[++i]);
        }
        else if ((arg == "--ficffaA") || (arg == "-i"))
        {
            ficffaA = argv[++i];
        }
        else if (arg == "--ficffaB")
        {
            ficffaB = argv[++i];
        }
        else if (arg == "-o")
        {
            ficffaA = ficech = ficffaB = ficdotB = ficpsB = argv[++i];
        }
        else if (arg == "--moderet")
        {
            arg = argv[++i];
            if (arg == "begin")
                moderet = begin;
            else
                moderet = end;
        }
        else if (arg == "--modeeps")
        {
            arg = argv[++i];
            if (arg == "fixed")
                modeeps = epsfixed;
            else if (arg == "variable")
                modeeps = variable;
			else
				cerr << "option for --modeeps can be 'fixed' or 'variable', here it is '" << arg << endl;
        }
        else if (arg == "--ficdotA")
        {
            ficdotA = argv[++i];
        }
        else if (arg == "--ficdotB")
        {
            ficdotB = argv[++i];
        }
        else if (arg == "--ficech")
        {
            ficech = argv[++i];
        }
        else if (arg == "--ficechpre")
        {
            ficechpre = argv[++i];
        }
        else if ((arg == "--precision") || (arg == "-p"))
        {
            precision = atof (argv[++i]);
        }
        else if (arg == "--epsprime")
        {
            epsprime = atof (argv[++i]);
        }
		else if (arg == "--seuil")
		{
			seuil = atoi (argv[++i]);
		}
        else if (arg == "--numPFA")
        {
            numPFA = atoi (argv[++i]);
        }
		else if (arg == "--seuilbm") 
		{
			seuilbm = atof (argv[++i]);
		}
		else if (arg == "--nb_tours") 
		{
			nb_tours = atoi (argv[++i]);
		}
		else if (arg == "--ficps")
        {
            ficpsA = argv[++i];
        }
        else if (arg == "--ficpsB")
        {
            ficpsB = argv[++i];
        }
		else if (arg == "--stepssave")
		{
			stepssave=true ;
		}
        else if (arg == "--nb_etats")
        {
            nb_etats = atoi (argv[++i]);
        }
		else if (arg == "--nb_essais") {
			nb_essais = atoi (argv[++i]);
		}
		else if (arg == "--max_states") {
			max_states = atoi (argv[++i]);
		}
        else if (arg == "--nb_lettres")
        {
            nb_lettres = atoi (argv[++i]);
        }
        else if (arg == "--precpow")
        {
            precpow = atof (argv[++i]);
        }
        else if (arg == "--affiche_A")
        {
            affiche_A = true;
        }
        else if (arg == "--affiche_B")
        {
            affiche_B = true;
        }
        else if (arg == "--blind")
        {
            affiche_A = affiche_B = false;
        }
        else if ((arg == "--taille_ech_test") || (arg == "-tt"))
        {
            taille_ech_test = atoi (argv[++i]);
        }
        else if (arg == "--name")
        {
            nom_A = argv[++i];
        }
        else if (arg == "--nom_B")
        {
            nom_B = argv[++i];
        }
        else if (arg == "--test_proportion")
        {
            test_proportion = atof (argv[++i]);
        }
		else if (arg == "--smooth") {
			smooth=true;
		}
		else if (arg == "--rand") {
			srand( atoi(argv[++i]) );			
		}
		else if (arg == "--viewprec") {
			viewprec=atoi(argv[++i]);
		}
		else if (arg == "--viewstates") {
			viewstates=true;
		}
		else if (arg == "--viewmultipletransitions") {
			viewmultipletransitions=true;
		}
        else
        {
            cerr << "unknow option : " << arg << endl;
        }
        i++;
    }
    return i;
}

// La procedure principale
int
main (int argc, char *argv[])
{
    int i = 0;
	string commande;
	commande = "mkdir " + temporary_directory + ">/dev/null 2>&1";
	system(commande.c_str());
    try
    {
        if (argc < 2)
            throw 1;
        randomize ();
        i = initialise_arguments (argc, argv);

        string arg;
        arg = argv[1];
        if (arg == "-h" || arg == "--help")
        {		// ---------- help ---------------
            throw 1;
        }
        else if (arg == "-H")
        {		// ----------------- HELP with options -------
            affiche_utilisation ();
            affiche_options ();
            cout << "Report bugs or suggestions to Yann Esposito <esposito@cmi.univ-mrs.fr>." << endl;
        }
        else if ((arg == "-Y") || (arg == "-I"))
        { // -------------- interactive mode -------------
            Y ();
        }
		else if (arg == "--test")
		{
			test();
		}
		else if (arg == "--PSe")
		{
            ficffaA = argv[i];
			if (!OK(A.load(ficffaA))) {
				throw string ("Erreur d'ouverture du fichier ") + ficffaA;
			}
			affiche_PSe();
		}
        else if (arg == "--export")
        {		// ---------------- exportation----------------
            if (argc < i + 3) {
				cerr << "export command export some automaton file to another format" << endl;
				cerr << "available format are dot, pdf and ffa." << endl;
				cerr << "usage: dees --export [OPTIONS] format input output" << endl;
                throw 2;
			}
			format = argv[i];
            ficffaA = argv[i+1];
            ficffaB = argv[i + 2];
            if (!OK(A.load (ficffaA))) {
				throw string ("Erreur d'ouverture du fichier ") + ficffaA;
			}
            A.name = nom_A;
            if (format == "dot")
            {	// exportation to dot format
                A.save_dot (ficffaB);
            }
            else if (format == "pdf")
            {	// exportation to ps format
                export_to_pdf (A, ficffaB, temporary_directory + "tmp.dot");
            }
            else if (format == "ffa")
            {	// exportation ffa
                A.save (ficffaB);
            }
            else
            {
				cerr << "export command export some automaton file to another format" << endl;
				cerr << "available format are dot, ps and ffa." << endl;
				cerr << "usage: dees --export [OPTIONS] format input output" << endl;
                throw 4;
            }
        }
        else if (arg == "--affiche")
        {		// ------------------------ affichage -------------------
            if (argc < i + 1) {
				cerr << "--affiche command show the automaton on standart output" << endl;
				cerr << "usage: dees --affiche [OPTIONS] automaton" << endl;
                throw 2;
			}
            A.load (argv[i]);
            A.affiche ();
        }
        else if ((arg == "--showps") || (arg == "-P"))
        {		// --- affichage du postscript ------
            if (argc < i + 1) {
				cerr << "showps or -P command show the automaton in gv" << endl;
				cerr << "you must have dot (part of Graphviz project) and a script named open" << endl;
				cerr << "which is a pdf viewer installed" << endl;
				cerr << "usage: dees <--showps or -P> [OPTIONS] automaton" << endl;
				cerr << "usefull option: --name name" << endl;
                throw 2;	// pas assez d'arguments
			}
            A.load (argv[i]);
            A.name = nom_A;
			show(A);
        }
        else if (arg == "--dist")
        {		// ----------------------- calcul des distances ------------
            if (argc < i + 2) {
				cerr << "dist command show some distances between two automata." << endl;
				cerr << "usage: dees --dist [OPTIONS] automaton_A automaton_B" << endl;
				cerr << "usefull option: --taille_ech_test size" << endl;
				cerr << " --smooth or --no_smooth" << endl;
                throw 2;
			}
            A.load (argv[i]);
            B.load (argv[i + 1]);
			if (smooth == true) {
				B.lisse ();
			}
            affiche_distance (A, B, taille_ech_test, "");
        }
		else if (arg == "--deletetransitions") {
			if (argc < i + 3) {
				cerr << "dees --deletetransitions input_pfa output_pfa min max\n";
				cerr << "delete transitions of the input_pfa which values are between min and max\n";
				cerr << "renormalise and save it in output_pfa file.\n";
				return -1;
			}
			PFA A;
			string input=argv[i++];
			string output=argv[i++];
			double min=atof(argv[i++]);
			double max=atof(argv[i++]);
			A.load(input);
			if (max < min) {
				double tmp = max;
				max=min;
				min=tmp;
				cerr << "Warning max < min !" << endl;
			}
			A.erase_transitions(max,min);
			A.rend_PFA();
			A.save(output);
		}
        else if (arg == "--sample")
        {		// --------------------- generation d'un echantillon ---
            if (argc < i + 3) {
				cerr << "sample command generate a sample from an automaton." << endl;
				cerr << "usage: dees --sample [OPTIONS] automaton size sample" << endl;
				cerr << "usefull options:" << endl;
				cerr << "--format <ffa, alergia or mdi>" << endl;
				cerr << "--num_ech num_sample" << endl;
                throw 2;
			}
            ficffaA = argv[i];			
            taille_ech = atoi (argv[i + 1]);
            ficech = argv[i + 2];
            A.load (ficffaA);
            A.genere_echantillon (taille_ech, S, num_ech);
            if (format == "ffa")
                S.save (ficech, ffa);
            else if (format == "alergia")
                S.save (ficech, alergia);
            else if (format == "mdi")
                S.save (ficech, mdi);
            else
                throw 4;
        }
		else if (arg == "--mdi") {
            if (argc > i)
            {
                ficech = argv[i];
            }
            else
            {
                throw 1;
            }
            if (verbose)
                cout << "loading " << ficech << "..." << endl;

			MDI(A,precision,ficech,verbose);
					
            if (verbose)
                cout << "saving A in " << ficffaA << endl;
            A.save (ficffaA);
            // affichage de l'automate
            if (affiche_A)
            {
                show(A);
            }
		}
		else if (arg == "--alergia") {
			if (argc <= i) 
				throw 1;
			
			ficech = argv[i];
            if (verbose)
                cout << "loading " << ficech << "..." << endl;

			Alergia(A,precision,ficech,verbose);
		
            if (verbose)
                cout << "saving A in " << ficffaA << endl;
			
            A.save (ficffaA);
            // affichage de l'automate
            if (affiche_A)
            {
                show(A);
            }
		}
        else if (arg == "--deesdet")
        {
			
            if (argc > i)
            {
                ficech = argv[i];
            }
            else
            {
				cerr << "--dees command learn an automaton using dees." << endl;
				cerr << "usage: dees --dees [OPTIONS] sample" << endl;
				cerr << "usefull options:" << endl;
				cerr << "-p or --precision float\tthe precision parameter" << endl;
				cerr << "--epsprime float\tnumber under which a transition is deleted" << endl;
				cerr << "-v or --verbose" << endl;
				cerr << "--moderet [begin or end]" << endl;
				cerr <<  "--max_states number" << endl;
				cerr << "--seuil seuil" << endl;
				cerr << "--stepssave\tsave steps of the algorithm" << endl;
                throw 2;
            }
            if (verbose)
                cout << "loading " << ficech << "..." << endl;
            if (!OK (S.load (ficech.c_str ())))
            {
                cerr << "impossible d'ouvrir " << ficech << endl;
                throw 5;
            }
			
			epsprime = precision*pow(double(S.size()),-0.25)/5;
			epsprime = min(sqrt(precision),0.1);
			
            A.DEES (determinist, S, precision, epsprime, verbose, moderet, modeeps, max_states, seuil,10,0,true,stepssave);
            if (verbose)
                cout << "saving A in " << ficffaA << endl;
            A.save (ficffaA);
            // affichage de l'automate
            if (affiche_A)
            {
                show(A);
            }
        }		
        else if (arg == "--dees")
        {

            if (argc > i)
            {
                ficech = argv[i];
            }
            else
            {
				cerr << "--dees command learn an automaton using dees." << endl;
				cerr << "usage: dees --dees [OPTIONS] sample" << endl;
				cerr << "usefull options:" << endl;
				cerr << "-p or --precision float\tthe precision parameter" << endl;
				cerr << "--epsprime float\tnumber under which a transition is deleted" << endl;
				cerr << "-v or --verbose" << endl;
				cerr << "--moderet [begin or end]" << endl;
				cerr <<  "--max_states number" << endl;
				cerr << "--seuil seuil" << endl;
				cerr << "--stepssave\tsave steps of the algorithm" << endl;
                throw 2;
            }
            if (verbose)
                cout << "loading " << ficech << "..." << endl;
            if (!OK (S.load (ficech.c_str ())))
            {
                cerr << "impossible d'ouvrir " << ficech << endl;
                throw 5;
            }

			epsprime = precision*pow(double(S.size()),-0.25)/5;
			epsprime = min(sqrt(precision),0.1);

            A.DEES (positive, S, precision, epsprime, verbose, moderet, modeeps, max_states, seuil,10,0,true,stepssave);
            if (verbose)
                cout << "saving A in " << ficffaA << endl;
            A.save (ficffaA);
            // affichage de l'automate
            if (affiche_A)
            {
                show(A);
            }
        }
		else if (arg == "--deesbm") {
			// --- DEES version Baum Welch ---
            if (argc > i)
            {
                ficech = argv[i];
            }
            else
            {
				cerr << "--deesbm command learn an automaton using deesbm." << endl;
				cerr << "usage: dees --deesbm [OPTIONS] sample" << endl;
				cerr << "usefull options:" << endl;
				cerr << "-p or --precision number" << endl;
				cerr << "-v or --verbose" << endl;
				cerr << "--max_states number\tmaximal number of states" << endl;
				cerr << "--seuil number\tminimal number of suffix do use residual" << endl;
				cerr << "--nb_tours number\tmax number of turns for Baum Welch" << endl;
				cerr << "--seuilbm double\tprecision under which Baum Welch is considered to found the Max Likelihood Model" << endl;
				cerr << "-o or --ficffaA file\ttarget file" << endl;
                throw 2;
            }
            if (verbose)
                cout << "loading " << ficech << "..." << endl;
            if (!OK (S.load (ficech.c_str ())))
            {
                cerr << "impossible d'ouvrir " << ficech <<
                endl;
                throw 5;
            }
						
			epsprime = precision*pow(double(S.size()),-0.25)/5;

            A.DEESBM (S, precision, epsprime, verbose, max_states, seuil, seuilbm, nb_tours);
            if (verbose)
                cout << "saving A in " << ficffaA << endl;
			// A.emmonde();
            A.save (ficffaA);
            // affichage de l'automate
            if (affiche_A)
            {
                show(A);
            }		
		}
		else if (arg == "--deesha")
        {
			
            if (argc > i)
            {
                ficech = argv[i];
            }
            else
            {
                throw 1;
            }
            if (verbose)
                cout << "loading " << ficech << "..." << endl;
            if (!OK (S.load (ficech.c_str ())))
            {
                cerr << "impossible d'ouvrir " << ficech <<
                endl;
                throw 5;
            }
			
			// epsprime = precision/(2*sqrt(sqrt(n)))
			// epsprime est donc un o(precision)
			epsprime = precision*pow(double(S.size()),-0.25)/5;
			
            A.DEES (nonconstrained, S, precision, epsprime, verbose, moderet,modeeps, max_states, seuil, maxmots);
            if (verbose)
                cout << "saving A in " << ficffaA << endl;
			// A.emmonde();
            A.save (ficffaA);
            // affichage de l'automate
            if (affiche_A)
            {
                show(A);
            }
        }		
        else if ((arg == "--generate_complete_random_pfa") || (arg == "--gcrp"))
        {		// --- generation d'un pfa aléatoire complet ---
			generate_complete_random_pfa(argc, argv, i);
        }
        else if ((arg == "--generate_random_pfa") || (arg == "--grp"))
        {
			generate_random_pfa(argc, argv, i);
		}
		else if (arg == "--proba") { // show values of words for an MA
			show_sample_proba(argc, argv, i);
		}
		else if (arg == "--forward") { // show prefix values of words for an MA
			show_sample_forward(argc, argv, i);
		}
        else if (arg == "--class")
        {		// --- show the class of the automaton
            if (argc < 3)
                throw 2;
            A.load (argv[2]);
            affiche_classe (A);
        }
		else if (arg == "--bm") { // Baum Welch Algorithm
			if (argc < i+4) {
				cerr << "dees --bm pfa sample number_of_turns output_pfa\n";
				cerr << "work baum welch using the pfa for initialisation (structure + parameters)\n";
				cerr << "remember Baum Welch never change a 0 valued parameter\n";
				return -1;
			}
			string pfa=argv[i++];
			string sample=argv[i++];
			int nb_turns=atoi(argv[i++]);
			string output=argv[i++];
			
			A.load(pfa);
			S.load(sample);
			BaumWelch(A,S,nb_turns,output);
		}
		else if (arg == "--convert") {
			if (argc < i+3) {
				cerr << "dees --convert format sample output\n";
				cerr << "Change the format of the sample to output.\n";
				cerr << "format could be dees for internal, alergia or mdi\n";
				return -1;
			}
			string format=argv[++i];
			string input=argv[++i];
			string output=argv[++i];
			S.load(input);
			if (format == "alergia") {
				S.save(output, alergia);
			}
			else if (format == "mdi") {
				S.save(output, mdi);
			}
			else {
				S.save(output);
			}
		}
        else
        {		// ---------- choix non traitÈ -------------
            throw 1;
        }
        return EXIT_SUCCESS;
    }
    catch (int err)
    {
        switch (err)
        {
        case 1:
            affiche_utilisation ();
            cout << "Report bugs or suggestions to Yann Esposito <esposito@cmi.univ-mrs.fr>." << endl;
            break;
        case 2:
            cerr << "\nNo enougth arguments" << endl;
            break;
        case 3:
            cerr << "Input error, I cannot read the file." <<
            argv[i] << endl;
            break;
        case 4:
            cerr << format << "is an unknown format." << endl;
            break;
        case 5:
            cerr << "I/O Error !!!" << endl;
            break;
        default:
            cerr << "unknown error occured !" << endl;
        }

        if (err != 1)
        {
            cerr << "use -h, --help or -H option to display help" << endl;
        }
    }
	catch (string erreur) {
		cerr << erreur << endl;
	}
    catch (...)
    {
        cerr << "Unknown error !!!" << endl;
        affiche_utilisation ();
        return -1;
    }
}
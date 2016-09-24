#include "interface.H"
extern PPRFA A;
extern Sample S;
extern int nb_tours;
extern float epsprime;
extern bool verbose;
void genere_aleatoire_unique()
{
    int nbetats;
    int nblettres;
    int numgraphe;
    int nb_succ, nb_init, nb_term;
    PFA A;

    cout << "Nombre d'Ètat : ";
    cin >> nbetats;
    cout << "Nombre de lettres : ";
    cin >> nblettres;
    cout << "Numero du graphe : ";
    cin >> numgraphe;
    cout << "nombre de successeurs par couple (etat-lettre) : ";
    cin >> nb_succ;
    cout << "nombre maximal d'etats initiaux : ";
    cin >> nb_init;
    cout << "nombre maximal d'etats terminaux : ";
    cin >> nb_term;
    A.becomeRandomMax (nbetats, nblettres, numgraphe, nb_succ, nb_init, nb_term);
    string fichier_ps = "_tmp.ps";
    string fichier_tmp = "_temporary_file";
    export_to_pdf (A, fichier_ps, fichier_tmp);
    show_ps (fichier_ps);

    if (A.isPFA ())
    {
        if (A.isPRFA (true))
            cout << "is PRFA" << endl;
        else
            cout << "is PFA" << endl;
    }
    else
    {
        cout << "not PFA" << endl;
    }

}

void proportion_PRFA()
{
    //------------ Proportion de PRFA -----------
    int max_nb_states;	// max number of states
    int nb_try;			// number of tests
    map < int, int >Tableau;
    map < int, int >nbPFA;
    int i;		// taille du tableau
    int j;		// nombre d'essais
    int nb_moy_succ; // nombre moyen de successeurs
    unsigned int taille_auto;	// taille de l'automate
    
    cout << "Nb d'Ètat maximal : ";
    cin >> max_nb_states;
    cout << "Nb de PFAs par taille : ";
    cin >> nb_try;
    cout << "Nb moyen de successeurs par Ètat et par lettre : ";
    cin >> nb_moy_succ;
    for (i = 1; i <= max_nb_states; ++i)
    {
        nbPFA[i]=0;
        Tableau[i]=0;
        cout << "\ni = " << i << endl;
        for (j = 0; j < nb_try; ++j)
        {
            A.becomeRandomMax(i,2,0,nb_moy_succ, 1, i);
            taille_auto = (int) A.size ();
            if (taille_auto == (unsigned int) i)
            {
                if (A.isPFA ())
                {
                    ++nbPFA[taille_auto];
                    if (A.isPRFA (false))
                    {
                        ++Tableau [taille_auto];
                        cout << "R" << taille_auto << ",";
                    }
                    else
                    {
                        cout << "P" << taille_auto << ",";
                    }
                }
                else
                {
                    cout << "X";
                }
            }
        }
    }

    // compte rendu
    cout << "Nb d'Ètats\t% de PRFA\t(nb d'essais)" << endl;
    map < int, int >::iterator t;
    map < int, int >::iterator nbPFAi;
    for (t = Tableau.begin (); t != Tableau.end (); ++t)
    {
        if (nbPFA[t->first] != 0)
        {
            cout << t->first << "\t\t" << (t->second * 100) /
            nbPFA[t->first] << "%\t(" << t->
            second << "/" << nbPFA[t->first] << ")" << endl;
        }
    }
    //----------  Fin de la procedure de test de la proportion de PRFA

}

void prefixialisa_echantillon()
{
    string ficin, ficout;
    cout << "Entrez le fichier echantillon : ";
    cin >> ficin;
    cout << "Entrez le nom du fichier de sortie : ";
    cin >> ficout;
    cout << "load en cours..." << endl;
    S.load (ficin);
    cout << "save en cours..." << endl;
    S.prefixialise();
    S.save (ficout);
    cout << "Voil‡ c'est fini" << endl;
}

void Choix_cysteines()
{
    // = 1 ===================== on entre les noms des fichiers
    string nom_base_pos;
    cout << "Veuillez entrer les noms des fichiers de pfa,\n";
    cout << "de l'echantillon positif associÈ,\n";
    cout << "de l'Èchantillon nÈgatif associÈ,\n";
    cout << "et enfin le nom de l'Èchantillon test associÈ.\n\n";

    cout << "Finissez en entrant le nom ." << endl;
    ;

    map < int, PPRFA > Modele;
    map < int, Sample > PosSample;
    map < int, Sample > NegSample;
    map < int, Sample > TestSample;
    string nom;
    Sample S;
    PPRFA A;
    cout << "Nom du fichier pfa : ";
    cin >> nom;
    int i;
    i = -1;
    while (nom != ".")
    {
        ++i;
        A.load (nom);
        Modele[i] = A;
        cout << "Èchantillon positif associe : ";
        cin >> nom;
        S.load (nom);
        PosSample[i] = S;
        cout << "Èchantillon negatif associe : ";
        cin >> nom;
        S.load (nom);
        NegSample[i] = S;
        cout << "Èchantillon test associe : ";
        cin >> nom;
        S.load (nom);
        TestSample[i] = S;
        cout << "\nNom du fichier pfa : ";
        cin >> nom;
    }

    // = 2 =============== on calcule les moyennes et ecart-types
    map < int, float >moyenne_pos;
    map < int, float >ecartype_pos;
    map < int, float >moyenne_neg;
    map < int, float >ecartype_neg;
    Dictionnaire dico;

    map < int, float >probabilite_pos;
    map < int, float >probabilite_neg;
    map < int, float >::const_iterator prob_it;
    map < int, float >decalage_pos;
    map < int, float >decalage_neg;
    Sample::const_iterator w;

    int j;
    float min;
    float val;
    for (i = 0; i < (int) Modele.size (); ++i)
    {		// pour tous les modËles

        Modele[i].lisse();

        // 1a - on calcule la probabilite minimale les mots positifs
        dico = PosSample[i].dictionnaire ();
        j=0;
        min = 1;
        for (w = PosSample[i].begin ();
                w != PosSample[i].end (); ++w)
        {
            probabilite_pos[j] = Modele[i].plog (w->first, &dico);
            if (probabilite_pos[j] < min)
                min = probabilite_pos[j];
            ++j;
        }
        decalage_pos[i] = min;

        // 1b - on calcule la probabilite minimale les mots negatifs
        dico = NegSample[i].dictionnaire ();
        j=0;
        min = 1;
        for (w = NegSample[i].begin (); w != NegSample[i].end (); ++w)
        {
            probabilite_neg[j] = Modele[i].plog (w->first, &dico);
            if (probabilite_neg[j] < min)
                min = probabilite_neg[j];
            ++j;
        }
        decalage_neg[i]= min;



        if (decalage_neg[i]<decalage_pos[i])
            decalage_pos[i]=decalage_neg[i];
        else
            decalage_neg[i]=decalage_pos[i];

        // --- on a calcule le decalage ---

        // 2 - on calcule la moyenne pour les mots positifs
        moyenne_pos[i] = 0;
        for (prob_it  = probabilite_pos.begin() ;
                prob_it != probabilite_pos.end(); ++prob_it)
        {
            moyenne_pos[i] += exp (prob_it->second - decalage_pos[i]);
        }
        moyenne_pos[i] /= PosSample[i].size();

        // 3 - on calcule l'ecartype pour les mots positifs
        ecartype_pos[i] = 0;
        for (prob_it  = probabilite_pos.begin() ;
                prob_it != probabilite_pos.end() ; ++prob_it)
        {
            ecartype_pos[i] +=
                pow (moyenne_pos[i] - exp (prob_it->second - decalage_pos[i]), 2);
        }
        ecartype_pos[i] /= PosSample[i].size();
        ecartype_pos[i] = sqrt(ecartype_pos[i]);

        // Maintenant, moyenne_pos[i]  contient la valeur moyenne de la probabilite
        // des mots positifs et ecartype_pos[i] contient l'ecart-type des mots pos.

        // puis pareil pour les nÈgatifs

        // 2 - on calcule la moyenne
        moyenne_neg[i] = 0;
        for (prob_it  = probabilite_neg.begin() ;
                prob_it != probabilite_neg.end(); ++prob_it)
        {
            moyenne_neg[i] += exp (prob_it->second - decalage_neg[i]);
        }
        moyenne_neg[i] /= NegSample[i].size();

        // 3 - on calcule l'ecartype pour les mots positifs
        ecartype_neg[i] = 0;
        for (prob_it  = probabilite_neg.begin() ;
                prob_it != probabilite_neg.end() ; ++prob_it)
        {
            ecartype_neg[i] +=
                pow (moyenne_neg[i] - exp (prob_it->second - decalage_neg[i]), 2);
        }
        ecartype_neg[i] /= NegSample[i].size();
        ecartype_neg[i] = sqrt( ecartype_neg[i] );

        // on a fini de calculer les moyennes et les ecartype pour le modele i

        cout << "µ pos : " << moyenne_pos[i];
        cout << ", sigma pos : " << ecartype_pos[i];
        cout << ", decalage pos : " << decalage_pos[i] << "\n";
        cout << "µ neg : " << moyenne_neg[i];
        cout << ", sigma neg : " << ecartype_neg[i];
        cout << ", decalage neg : " << decalage_neg[i] << endl;
    }

    // = 3 ==================================== evaluation

    double p1, p2;
    bool res;
    map < int, float >probapos;
    Word u;
    res = 0;
    i = 0;
    for (i = 0; i < (int) Modele.size (); ++i)
    {
        u = TestSample[i].begin ()->first;
        dico = TestSample[i].dictionnaire ();
        val = Modele[i].plog (u, &dico);
        cout << "val = " << exp( val - decalage_pos[i] ) << endl;
        if (moyenne_pos[i] > moyenne_neg[i])
        {
            p1 = fonction (exp (val - decalage_pos[i]), moyenne_pos[i], ecartype_pos[i]);	// ECRIRE fonction
            p2 = fonction (moyenne_neg[i], exp (val - decalage_neg[i]) , ecartype_neg[i]);	// !!!!!!!!!!!!!!!
        }
        else
        {
            p1 = fonction (moyenne_pos[i], exp (val - decalage_pos[i]), ecartype_pos[i]);	// ECRIRE fonction
            p2 = fonction (exp (val - decalage_neg[i]), moyenne_neg[i] , ecartype_neg[i]);	// !!!!!!!!!!!!!!!
        }
        cout << "p_pos = " << p1;
        cout << ", p_neg = " << p2;
        probapos[i] = (p1 + (1 - p2)) / 2;	// ???????
        cout << ", probapos = " << probapos[i] << endl;
    }
    float probaposgenerale=0;
    for (i=0 ; i< (int) probapos.size() ; ++i)
        probaposgenerale += probapos[i];
    probaposgenerale /= i;
    cout << "probabilite d'etre positif = " << probaposgenerale << endl;

}

void Apprentissage_unique()
{
    string ficech;
    cout << "Veuillez entrer le nom du fichier contenant l'echantillon : ";
    cin >> ficech;
    Sample S;
    if (!OK (S.load (ficech.c_str ())))
    {
        cerr << "impossible d'ouvrir " << ficech << endl;
        throw 5;
    }

    PPRFA A;
    string reponse;
    unsigned int seuil;
    string str;
    float tmp;
    float precision;
    T_ModeReturn moderet=::begin;
    T_ModeEpsilon modeeps=epsfixed;

    int max_etats=INT_MAX;
    cout << "Les parametres necessaires a l'apprentissage sont : \n";
    while (reponse != "0")
    {
        cout << "1 - La precision .................. : "
        << precision << "\n";
        cout << "2 - Le mode de retour ............. : "
        << (moderet==::begin?"begin":"end") << "\n";
        cout << "3 - Le mode de gestion de epsilon . : "
        << (modeeps==epsfixed?"epsfixed":
            (modeeps==variable?"variable":"word_variable"))
        << "\n";
        cout << "4 - Seuil de suppression des aretes : "
        << epsprime << "\n";
        cout << "5 - Verbose ....................... : "
        << verbose << "\n";
        cout << "6 - La taille maximale de la cible  : "
        << max_etats << "%\n";
        cout << "7 - Nombre de successeurs par mot   : "
        << seuil << "\n";
        cout << "\n";
        cout << "0 - Debut de l'apprentissage\n";
        cout << "\n";
        cout << "Votre choix : ";
        cin >> reponse ;

        if (reponse == "1")
        { // precision
            cout << "nouvelle precision : ";
            cin >> precision;
        }
        else if (reponse == "2")
        { // mode de retour
            cout << "mode de retour (begin|end|b|e) : ";
            cin >> str;
            if ((str=="end") || (str == "e"))
                moderet=::end;
            else
                moderet=::begin;
        }
        else if (reponse == "3")
        { // mode de gestion de epsilon
            cout << "mode de gestion d'epsilon (fixed|variable|word_variable|f|v|w) : ";
            cin >> str;
            if ((str == "word_variable") || (str == "w"))
            {
                modeeps=word_variable;
            }
            else if ((str == "variable") || (str == "v"))
            {
                modeeps=variable;
            }
            else
            {
                modeeps=epsfixed;
            }
        }
        else if (reponse == "4")
        { //seuil de suppression des aretes
            cout << "seuil de suppression des aretes ([0;1]) : ";
            cin >> tmp;
            if ((tmp >=0) && (tmp <= 100))
                epsprime = tmp;
        }
        else if (reponse == "5")
        { // verbose
            verbose = !verbose;
        }
        else if (reponse == "6")
        { // nombre maximal d'etats
            cout << "Nombre maximal d'etats : ";
            cin >> max_etats;
        }
        else if (reponse == "7")
        { // nombre de successeur par mots
            cout << "Nombre de successeurs necessaires par mot : ";
            cin >> seuil;
        }
        else if (reponse == "0")
        {}
        else
        {
            cout << "Choix inconnu" << endl;
        }
    }

    A.DEES (nonconstrained, S, precision, epsprime, verbose, moderet, modeeps, max_etats, seuil);
    string ficffaA="test/A.ffa";
    if (verbose)
        cout << "saving A in " << ficffaA << endl;
    A.save (ficffaA);
    // affichage de l'automate
    string ficpsA="test/A.ps";
    string ficdotA="test/A.dot";
    export_to_pdf (A, ficpsA, ficdotA);
    show_ps (ficpsA);
}

void makePTA()
{
    string ficech;
    PPRFA A;
    Sample S;
    cout << "Donnez le fichier contenant l'echantillon : ";
    cin >> ficech;
    S.load(ficech);
    A.becomePrefixTree(S);
    string ficpsA = "test/tmp.ps";
    string ficdotA= "test/tmp.dot";
    export_to_pdf(A,ficpsA, ficdotA);
    show_ps(ficpsA);
}

void TST_FD_BD() {
	string ficpfa;
	cout << "Donnez le fichier contenant le PFA : ";
	cin >> ficpfa;
	A.load(ficpfa);

	Word w;
	cout << "Donnez un mot : ";
	cin >> w;

	list<PreciseSFunc> F, B;
	A.logforward(F,w);
	A.logbackward(B,w);

	list<PreciseSFunc>::const_iterator f;
	list<PreciseSFunc>::const_reverse_iterator b;
	PreciseSFunc::const_iterator x, y;
	cout << "P(" << w << ")=" << A.p(w) << endl;
	for (f = F.begin(), b=B.rbegin() ; f != F.end() ; ++f, ++b) {
		for (x=f->begin() ; x != f->end() ; ++x) {
			cout << "F[" << x->first << "]=" << exp(x->second) << ", ";
		}
		cout << "\n";
		for (y=b->begin() ; y != b->end() ; ++y) {
			cout << "B[" << y->first << "]=" << exp(y->second) << ", ";
		}
		cout << "\n\t----" << endl;
	}

}

void TST1() {
	S.load("sample");
	SFunc solution;
	WordSFunc XR;
	XR[""]=1;
	XR["a"]=2;
	XR["b"]=3;
	solution[1]=0.5;
	solution[2]=0.5;
	A.becomePrefixTree(S);
	SFunc Iota,Tau;
	A.allStates(Tau);
	TransitionFunction T;
	A.allTransitions(T);
	A.BaumWelch(S,T,Iota,Tau,10 , false);
	float x = S.AutoLikelihood();
	float y = A.Likelihood(S);
	A.save("test/res.ffa");		
	cout << S.size() << " "<< x << " "<< y << " :: " << (x-y)/S.size() << endl;
}

void TST(){
	string entree;
	string sortie;
	cout << "Echantillon d'entrée : ";
	cin >> entree;
	cout << "Echantillon de sortie : ";
	cin >> sortie;
	Sample S;
	S.load(entree);
	S.save(sortie);
}

void TST4() {
	SFunc solution;
	SFunc Iota,Tau;
	TransitionFunction T;
	Transition t1,t2;
	float x, y, x1, x2, y1, y2;
	float xsum, val;
	int i;

	S.load("sample");	
	solution[1]=0.5;
	solution[2]=0.5;
	S.prefixialise();
	int nb_sucb=S["b"];

	cout << "prefixialisation..." << endl;
	S.deprefixialise();

	cout << "becomePrefixTree..." << endl;
	A.becomePrefixTree(S);
	A.save("test/pta");

	cout << "becomeQPFT..." << endl;
	A.becomeQuasiPrefixTree(S,"b",solution);
	A.save("test/qpta");

	A.allStates(Tau);
	A.allTransitions(T);

	t1.qdep=t2.qdep=1;
	t1.a=t2.a='b';
	t1.qarr=1;
	t2.qarr=2;
	

	A.BaumWelch(S,T,Iota,Tau,30, false);	
	x1=T[t1]; x2=T[t2];
	xsum = x1+x2;
	x1 /= xsum;
	x2 /= xsum;

	A.BaumWelch(S,T,Iota,Tau,1, false);
	y1=T[t1]; y2=T[t2];
	xsum = y1+y2;
	y1 /= xsum;
	y2 /= xsum;
	
	val=abs(x1-y1)+abs(x2 - y2);

	i=20;	
	while ((val > 0.001) && (--i != 0)) {			
		cout << ".";
		x1=y1; x2=y2;
		A.BaumWelch(S,T,Iota,Tau,1, true);
		y1=T[t1]; y2=T[t2];
		xsum = y1+y2;
		y1 /= xsum;
		y2 /= xsum;	
		val=abs(x1-y1)+abs(x2 - y2);
	}
	A.save("test/qptafinal.ffa");

	A.save("test/res.ffa");	
	x = S.AutoLikelihood();
	y = A.Likelihood(S);
	val=abs((y-x)/float(S.size()));
	cout << S.size() << " "<< x << " "<< y << " :: " << val << " " << nb_sucb;
	cout << " " << val - 5/sqrt(float(nb_sucb)) << endl;
}

void TST3() {
	S.load("sample");
	A.load("test/tst3.ffa");
	TransitionFunction T;
	SFunc Iota, Tau;
	A.allTransitions(T);
	A.allStates(Iota);
	A.allStates(Tau);
	for (int i=1 ; i<200 ; ++i) {
		A.BaumWelch(S,T,Iota,Tau,1,true);
		cout << "." << flush;
		cout << "\t" << A.Likelihood(S) << "\n";
	}
	cout << endl;

	A.erase_bad_states();
	A.save("test/res.ffa");
}

void TSTprime()
{

    string ficech="sample";
//    cout << "Donnez le fichier contentant l'echantillon : ";
//    cin >> ficech;
    S.load(ficech);

    float precision, epsprime;

    cout << "Donnez la precision : ";
    cin >> precision;
    cout << "Donnez le seuil de suppression des aretes : ";
    cin >> epsprime;
	cout << "Donnez le nombre de tours de BM : ";
	int nb_tours;
	cin >> nb_tours;

    A.DEESBM(S,precision,epsprime,true,INT_MAX,0,nb_tours);

    string ficpsA = "test/tmp.ps";
    string ficdotA= "test/tmp.dot";
    export_to_pdf(A,ficpsA, ficdotA);
    show_ps(ficpsA);
}

void proba(void) {
	string nomfic;
	Word w;
	cout << "entrez le fichier : ";
	cin >> nomfic;
	A.load(nomfic);
	cout << "entrez le mot : ";
	cin >> w;
	cout << "p(" << w << ")= " << A.p(w) << endl;
}

void
Y (void)
{
    int choix;
    cout << "\t\t1 - BaumWelch sur test/tst3.ffa\n";
    cout << "\t\t2 - Probabilite d'un mot\n";
    cout << "\t\t3 - Prefixialise un echantillon\n";
    cout << "\t\t4 - Choix pour les cysteines\n";
    cout << "\t\t5 - Apprend ‡ partir d'un Èchantillon avec optimisation du paramËtre\n";
    cout << "\t\t6 - Apprend ‡ partir d'un Èchantillon\n";
    cout << "\t\t7 - Cree un PTA ‡ partir d'un Èchantillon\n";
    cout << "\t\t8 - TST\n";
    cout << "\n";
    cout << "choix : ";
    cin >> choix;

    if (choix == 1)	// Une generation aleatoire
    {
        TST3();
    } else if (choix == 2)	// Proportion des PRFAs
    {
        proba();
    } else if (choix == 3)	// Prefixialise un echantillon
    {
        prefixialisa_echantillon();
    } else if (choix == 4)	// Test pondere "intelligent" ?
    {
        Choix_cysteines();
    } else if (choix == 5)	//
    {
        cout << "Ce choix n'est plus valide !!!" << endl;
    }
    else if (choix == 6)
    {
        Apprentissage_unique();
    }
    else if (choix == 7)
    {
        makePTA();
    }
    else if (choix == 8)
    {
        TST();
    }
    else			// mauvais choix
    {
        cout << "entree inconnue" << endl;
    }
}

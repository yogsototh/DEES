#ifndef _MAIN_H
#define _MAIN_H
#include "pprfa.H"
#include <sstream>

// affiche un mot
string affiche (const Word & w, const Dictionnaire & alph);
// La fonction qui va charger les info à partir d'un fichier
RESULT charge (char *nomFic);
// initialisation des arguments (ligne de commande)
int initialise_arguments (int argc, char *argv[]);
// La procédure qui initialise le random
int randomize (void);
int affiche_utilisation (void);
int affiche_options (void);
RESULT export_to_pdf (const MA & A,
                     const string & fic_ps,
                     const string & fic_tmp,
                     const int precision = 3,
                     const bool noms_etats = false,
                     const bool multiples_trans = false);
void show_ps (string & fic_ps);
// Affiche toutes les distances entre A et B
void affiche_distance (PFA & A, const PFA & B, const int taille_ech_test,
                       const string & flag);
void affiche_distance (PFA & A, const PFA & B, const int taille_ech,
                       const char *flag);
void affiche_classe (PFA & A);
// --- La procedure pour debbugger ---

float fonction (float val, float moy, float ecart);
// La procedure principale
int main (int argc, char *argv[]);
#endif

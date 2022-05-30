/* Quentin MOREAU
* TIPE Bioinformatique
* Reconstruction d'arbres phylogenetiques
*
* fichier d'entete construction_arbre_phylogenetique
*
* 2021-2022
*/



#ifndef CONSTRUCTION_ARBRE_PHYLOGENETIQUE_H_INCLUDED
#define CONSTRUCTION_ARBRE_PHYLOGENETIQUE_H_INCLUDED


#include "main.h"


typedef struct Arbre Arbre;
struct Arbre
{
	char* nom;
	Arbre** arbres;
};


void afficherArbre(Arbre* arbre);
void afficherArbreRecursif(Arbre* arbre, int espace);
Arbre* constructionArbreUPGMA(int** distances, char** noms, int taille);
Arbre* constructionArbreUPGMARecursif(int** distances, Arbre** arbres, int taille);



#endif // CONSTRUCTION_ARBRE_PHYLOGENETIQUE_H_INCLUDED

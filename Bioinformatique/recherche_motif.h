/* Quentin MOREAU
* TIPE Bioinformatique
* Reconstruction d'arbres phylogenetiques
*
* fichier entete recherche_motif
*
* 2021-2022
*/



#ifndef RECHERCHE_MOTIF_H_INCLUDED
#define RECHERCHE_MOTIF_H_INCLUDED

#include "main.h"






typedef struct AutomateSuffixes AutomateSuffixes;
struct AutomateSuffixes
{
	int etatFinal;
	int** transitions;
};





File* rechercheNaive(Nucleotide sequence[], Nucleotide motif[], int tailleSequence, int tailleMotif);
File* boyerMoore(Nucleotide sequence[], Nucleotide motif[], int tailleSequence, int tailleMotif);
AutomateSuffixes* constructionAutomateSuffixes(Nucleotide motif[], int taille);
File* rechercheAutomateSuffixes(AutomateSuffixes* automate, Nucleotide sequence[], int tailleSequence);



#endif // RECHERCHE_MOTIF_H_INCLUDED

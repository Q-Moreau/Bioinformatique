/* Quentin MOREAU
* TIPE Bioinformatique
* Reconstruction d'arbres phylogenetiques
*
* fichier d'entete main
*
* 2021-2022
*/



#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

#include "construction_arbre_phylogenetique.h"



typedef enum Nucleotide Nucleotide;
enum Nucleotide
{
	A, C, G, T
};



typedef struct Element Element;
struct Element
{
	int element;
	Element* precedent;
	Element* suivant;
};


typedef struct File File;
struct File
{
	Element* premierElement;
	Element* dernierElement;
};



void enfiler(File* file, int element);
int defiler(File* file);
void supprimer(File* file);
Arbre* arbrePhylogenetique(char* nomsFichiers[], char* noms[], int nombreSequence);



#endif // MAIN_H_INCLUDED

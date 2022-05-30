/* Quentin MOREAU
* TIPE Bioinformatique
* Reconstruction d'arbres phylogenetiques
*
* fichier d'entete alignement
*
* 2021-2022
*/


#ifndef ALIGNEMENT_H_INCLUDED
#define ALIGNEMENT_H_INCLUDED

#define TRANSITION 1
#define TRANSVERSION 3
#define COUT_INSERTION_DELETION 10


int alignementRecursif(Nucleotide sequence1[], Nucleotide sequence2[], int tailleSequence1, int tailleSequence2);
int coutRecursif(Nucleotide sequence1[], Nucleotide sequence2[], int indice1, int indice2);
int needlemanWunsch(Nucleotide sequence1[], Nucleotide sequence2[], int tailleSequence1, int tailleSequence2);
int** distancesLevenshtein(Nucleotide** sequences, int* tailles, int nombreSequence);


#endif // ALIGNEMENT_H_INCLUDED
/* Quentin MOREAU
* TIPE Bioinformatique
* Reconstruction d'arbres phylogenetiques
*
* fichier d'entete prediction_regions_codantes
*
* 2021-2022
*/



#ifndef PREDICTION_REGIONS_CODANTES_H_INCLUDED
#define PREDICTION_REGIONS_CODANTES_H_INCLUDED

#define TAILLE 300

#define STOP0 0
#define STOP1 1
#define STOP2 2
#define START 3


File* rechercheCodonsStart(Nucleotide sequence[], int taille, AutomateSuffixes** automates);
File* rechercheCodonsStop(Nucleotide sequence[], int taille, AutomateSuffixes** automates);
File* predictionORF(Nucleotide sequence[], int taille, AutomateSuffixes** automates);
Nucleotide* extraireSequencesCodantes(Nucleotide sequence[], int taille, AutomateSuffixes** automates);


#endif // PREDICTION_REGIONS_CODANTES_H_INCLUDED
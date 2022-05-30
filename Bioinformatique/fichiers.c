/* Quentin MOREAU
* TIPE Bioinformatique
* Reconstruction d'arbres phylogenetiques
*
* fichier source fichier
*
* 2021-2022
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "main.h"
#include "fichiers.h"



Nucleotide* lireSequence(char* nom)
{
	FILE* fichier = NULL;
	fopen_s(&fichier, nom, "r");
	
	if (fichier == NULL)
	{
		printf("impossible d'ouvrir %s", nom);
		exit(EXIT_FAILURE);
	}

	fseek(fichier, 0, SEEK_END);
	int tailleFichier = ftell(fichier);
	rewind(fichier);

	char test = 'a';
	int curseur = 0;
	while (test != '\n')
	{
		curseur++;
		test = fgetc(fichier);
	}

	int lettre = 0;

	int taille = (tailleFichier-1 - curseur) - (tailleFichier-1 - curseur) / 61; // 60 carac par ligne + 1 \n
	Nucleotide* sequence = malloc(sizeof(Nucleotide) * (taille + 1)); // on enregistre aussi la taille
	

	if (sequence == NULL)
	{
		printf("echec (lireSequence 1)");
		exit(EXIT_FAILURE);
	}


	sequence[0] = taille;
	int i = 1;

	while (lettre != EOF)
	{
		lettre = getc(fichier);
		if (lettre != '\n' && lettre != EOF)
		{
			switch (lettre)
			{
			case 'A':
				sequence[i] = A;
				break;
			case 'C':
				sequence[i] = C;
				break;
			case 'G':
				sequence[i] = G;
				break;
			case 'T':
				sequence[i] = T;
				break;
			default:
				printf("%c",lettre);
				printf("\n\nERROR\n\n");
			}
			i += 1;
		}
	}


	fclose(fichier);


	return sequence;
}



void afficherSequence(Nucleotide* sequence, int taille)
{
	for (int i = 0; i < taille; i++)
	{
		switch (sequence[i])
		{
		case A:
			printf("A");
			break;
		case C:
			printf("C");
			break;
		case G:
			printf("G");
			break;
		case T:
			printf("T");
			break;
		default:
			printf(" erreur " );
			break;
		}
	}
}
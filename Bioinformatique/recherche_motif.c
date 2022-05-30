/* Quentin MOREAU
* TIPE Bioinformatique
* Reconstruction d'arbres phylogenetiques
*
* fichier source recherche_motif
*
* 2021-2022
*/



#include <stdlib.h>
#include <stdio.h>

#include "recherche_motif.h"









/////////////////////////////////////////// RECHERCHE NAIVE ///////////////////////////////////////////





File* rechercheNaive(Nucleotide sequence[], Nucleotide motif[], int tailleSequence, int tailleMotif)
{
	File* indices = malloc(sizeof(File));

	if (indices == NULL)
	{
		exit(EXIT_FAILURE);
	}

	indices->dernierElement = NULL;
	indices->premierElement = NULL;

	for (int i = 0; i < tailleSequence - tailleMotif; i++)
	{
		int indiceMotif = 0;
		int test = 1;

		while (test)
		{
			if (indiceMotif<tailleMotif)
			{
				test = sequence[i + indiceMotif] == motif[indiceMotif];
				indiceMotif++;
			}
			else
			{
				test = 0;
			}
		}

		if (indiceMotif >= tailleMotif)
		{
			enfiler(indices, i);
		}
	}


	return indices;
}









/////////////////////////////////////////// BOYER MOORE ///////////////////////////////////////////






File* boyerMoore(Nucleotide sequence[], Nucleotide motif[], int tailleSequence, int tailleMotif)
{
	File* indices = malloc(sizeof(File));

	if (indices == NULL)
	{
		exit(EXIT_FAILURE);
	}

	indices->dernierElement = NULL;
	indices->premierElement = NULL;

	int table[4] = { -1, -1, -1, -1 };

	for (int i = 0; i < tailleMotif; i++)
	{
		table[motif[i]] = i;
	}


	for (int i = 0; i < tailleSequence + 1 - tailleMotif; i++)
	{
		int indiceMotif = tailleMotif - 1;
		int test = 1;

		while (test)
		{
			if (indiceMotif >= 0)
			{
				test = sequence[i + indiceMotif] == motif[indiceMotif];
				indiceMotif--;
			}
			else
			{
				test = 0;
			}
		}

		if (indiceMotif < 0) // le motif est apparu
		{
			enfiler(indices, i);
			if (i + tailleMotif < tailleSequence)
			{
				i += tailleMotif - table[sequence[i + tailleMotif]] - 1;
			}
		}

		else // le motif n'est pas apparu
		{
			i += max(0, indiceMotif-table[sequence[i + indiceMotif]] - 1);
		}


	}


	return indices;
}












/////////////////////////////////////////// AUTOMATE DES SUFFIXES ///////////////////////////////////////////








AutomateSuffixes* constructionAutomateSuffixes(Nucleotide motif[], int taille)
{

	AutomateSuffixes* automate = malloc(sizeof(AutomateSuffixes));
	int** transitions = malloc(sizeof(int*) * (taille+1));

	if (automate == NULL || transitions == NULL)
	{
		exit(EXIT_FAILURE);
	}

	automate->etatFinal = taille;
	automate->transitions = transitions;

	for (int i = 0; i <= taille; i++)
	{
		transitions[i] = malloc(sizeof(int) * 4);
		if (transitions[i] == NULL)
		{
			exit(EXIT_FAILURE);
		}
	}



	for (int i = 0; i <= taille; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			transitions[i][j] = 0;
		}
	}



	for (int i = 1; i <= taille; i++)
	{

		int ancienneTransition = transitions[i - 1][motif[i - 1]]; // enregistrement de l'ancienne transition


		transitions[i - 1][motif[i - 1]] = i;

		for (int j = 0; j < 4; j++)
		{
			transitions[i][j] = transitions[ancienneTransition][j];
		}

	}



	return automate;
}








File* rechercheAutomateSuffixes(AutomateSuffixes* automate, Nucleotide sequence[], int tailleSequence)
{

	File* indices = malloc(sizeof(File));

	if (indices == NULL)
	{
		exit(EXIT_FAILURE);
	}

	indices->dernierElement = NULL;
	indices->premierElement = NULL;





	int etat = 0;

	for (int i = 0; i < tailleSequence; i++)
	{
		etat = automate->transitions[etat][sequence[i]];

		if (etat == automate->etatFinal)
		{
			enfiler(indices, i + 1 - automate->etatFinal); // automate->etatFinal == taille du motif
		}
	}

	return indices;
}











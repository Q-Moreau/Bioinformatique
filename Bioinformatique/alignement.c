/* Quentin MOREAU
* TIPE Bioinformatique
* Reconstruction d'arbres phylogenetiques
*
* fichier source alignement
*
* 2021-2022
*/


#include <stdlib.h>
#include <stdio.h>

#include "main.h"
#include "alignement.h"



int alignementRecursif(Nucleotide sequence1[], Nucleotide sequence2[], int tailleSequence1, int tailleSequence2)
{
	return coutRecursif(sequence1, sequence2, tailleSequence1 - 1, tailleSequence2 - 1);
}


int coutRecursif(Nucleotide sequence1[], Nucleotide sequence2[], int indice1, int indice2)
{

	int const COUTS_SUBSTITUTIONS[4][4] = { {0,TRANSVERSION,TRANSITION,TRANSVERSION},
											{TRANSVERSION,0,TRANSVERSION,TRANSITION},
											{TRANSITION,TRANSVERSION,0,TRANSVERSION},
											{TRANSVERSION,TRANSITION,TRANSVERSION,0} };

	if (indice1 == 0 && indice2 == 0) return 0;
	if (indice1 == 0) return indice2 * COUT_INSERTION_DELETION;
	if (indice2 == 0) return indice1 * COUT_INSERTION_DELETION;

	return min(coutRecursif(sequence1, sequence2, indice1 - 1, indice2 - 1) + COUTS_SUBSTITUTIONS[sequence1[indice1]][sequence2[indice2]],
		min(coutRecursif(sequence1, sequence2, indice1, indice2 - 1) + COUT_INSERTION_DELETION,
			coutRecursif(sequence1, sequence2, indice1 - 1, indice2) + COUT_INSERTION_DELETION));
}






int** distancesLevenshtein(Nucleotide** sequences, int* tailles, int nombreSequence)
{
	int** distances = malloc(sizeof(int*) * nombreSequence);

	if (distances == NULL)
	{
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < nombreSequence; i++)
	{
		distances[i] = malloc(sizeof(int) * nombreSequence);
		if (distances[i] == NULL)
		{
			exit(EXIT_FAILURE);
		}

		for (int j = 0; j < nombreSequence; j++)
		{
			distances[i][j] = 0;
		}
	}

	for (int i = 1; i < nombreSequence; i++)
	{
		for (int j = 0; j < i; j++)
		{
			//distances[i][j] =  alignementRecursif(sequences[i], sequences[j], tailles[i], tailles[j]);
			distances[i][j] = needlemanWunsch(sequences[i], sequences[j], tailles[i], tailles[j]);
			distances[j][i] = distances[i][j];
		}
	}

	return distances;;
}









int needlemanWunsch(Nucleotide sequence1[], Nucleotide sequence2[], int tailleSequence1, int tailleSequence2)
{

	if (tailleSequence2 < tailleSequence1)
	{
		return needlemanWunsch(sequence2, sequence1, tailleSequence2, tailleSequence1);
	}


	int* colonnes[2] = { NULL,NULL };

	colonnes[0] = malloc(sizeof(int) * (tailleSequence1 + 1));
	colonnes[1] = malloc(sizeof(int) * (tailleSequence1 + 1));

	int const COUTS_SUBSTITUTIONS[4][4] = { {0,TRANSVERSION,TRANSITION,TRANSVERSION},
											{TRANSVERSION,0,TRANSVERSION,TRANSITION},
											{TRANSITION,TRANSVERSION,0,TRANSVERSION},
											{TRANSVERSION,TRANSITION,TRANSVERSION,0} };


	if (colonnes[0] == NULL || colonnes[1] == NULL)
	{
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i <= tailleSequence1; i++)
	{
		colonnes[0][i] = i * COUT_INSERTION_DELETION;
	}
	for (int i = 0; i <= tailleSequence1; i++)
	{
		colonnes[1][i] = i * COUT_INSERTION_DELETION;
	}



	int k = 0;

	while (k < tailleSequence2)
	{

		colonnes[k % 2][0] = (k+1) * COUT_INSERTION_DELETION;

		for (int i = 0; i < tailleSequence1; i++)
		{
			colonnes[k % 2][i + 1] = min(colonnes[(k + 1) % 2][i] + COUTS_SUBSTITUTIONS[sequence1[i]][sequence2[k]], /* \ */
				min(colonnes[k % 2][i] + COUT_INSERTION_DELETION, /* | */
					colonnes[(k + 1) % 2][i + 1] + COUT_INSERTION_DELETION)); /* - */
		}

		k++;
	}




	int cout = colonnes[(k + 1) % 2][tailleSequence1];



	free(colonnes[0]);
	free(colonnes[1]);


	return cout;

}
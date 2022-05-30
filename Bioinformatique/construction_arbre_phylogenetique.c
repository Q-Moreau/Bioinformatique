/* Quentin MOREAU
* TIPE Bioinformatique
* Reconstruction d'arbres phylogenetiques
*
* fichier source construction_arbre_phylogenetique
*
* 2021-2022
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "main.h"
#include "construction_arbre_phylogenetique.h"






/////////////////////////////////////////// AFFICHAGE ARBRE ///////////////////////////////////////////


void afficherArbre(Arbre* arbre)
{
	afficherArbreRecursif(arbre, 0);
}

void afficherArbreRecursif(Arbre* arbre, int espace)
{
	if (arbre != NULL)
	{
		for (int i = 0; i < espace; i++) printf(" ");
		printf(" -> ");
		printf(arbre->nom);

		if (arbre->arbres != NULL)
		{
			printf("\n");
			afficherArbreRecursif(arbre->arbres[0], espace + 4 + strlen(arbre->nom));
			printf("\n");
			afficherArbreRecursif(arbre->arbres[1], espace + 4 + strlen(arbre->nom));
			printf("\n");
		}

	}
}










/////////////////////////////////////////// UPGMA ///////////////////////////////////////////




Arbre* constructionArbreUPGMA(int** distances, char** noms, int taille)
{
	Arbre** arbres = malloc(sizeof(Arbre*) * taille);

	if (arbres == NULL)
	{
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < taille; i++)
	{
		arbres[i] = malloc(sizeof(Arbre));
		arbres[i]->arbres = NULL;
		arbres[i]->nom = noms[i];
	}
	
	return constructionArbreUPGMARecursif(distances, arbres, taille);
}


Arbre* constructionArbreUPGMARecursif(int** distances, Arbre** arbres, int taille)
{
	if (taille <= 1)
	{
		return arbres[0];
	}

	int indice1 = 1, indice2 = 0;
	int minimum = distances[1][0];

	// recherche de minimum
	for (int i = 1; i < taille; i++)
	{
		for (int j = 0; j < i; j++)
		{
			if (distances[i][j] < minimum) // on a indice2<indice1
			{
				indice1 = i;
				indice2 = j;
				minimum = distances[i][j];
			}
		}
	}

	// distances recalculees

	for (int i = 0; i < taille; i++)
	{
		distances[i][indice1] = (distances[i][indice2] + distances[i][indice1]) / 2;
		distances[indice1][i] = distances[i][indice1];
		distances[i][indice2] = distances[i][indice1];
		distances[indice2][i] = distances[i][indice1];
	}

	// insertion de la ligne taille-1 dans la ligne indice1
	for (int i = 0; i < taille; i++)
	{
		distances[indice1][i] = distances[taille - 1][i];
		distances[i][indice1] = distances[indice1][i];
	}
	


	if (arbres[indice1] == NULL || arbres[indice2] == NULL)
	{
		exit(EXIT_FAILURE);
	}

	Arbre** sousArbres = malloc(2 * sizeof(Arbre*));
	Arbre* nouvelArbre = malloc(sizeof(Arbre));

	if (sousArbres == NULL || nouvelArbre == NULL)
	{
		exit(EXIT_FAILURE);
	}

	sousArbres[0] = arbres[indice1];
	sousArbres[1] = arbres[indice2];

	nouvelArbre->nom = "ancetre commun";
	nouvelArbre->arbres = sousArbres;

	arbres[indice2] = nouvelArbre;
	arbres[indice1] = arbres[taille - 1];


	return constructionArbreUPGMARecursif(distances, arbres, taille - 1);



}





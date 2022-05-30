/* Quentin MOREAU
* TIPE Bioinformatique
* Reconstruction d'arbres phylogenetiques
*
* fichier source prediction_regions_codantes
*
* 2021-2022
*/



#include <stdio.h>
#include <stdlib.h>

#include "main.h"
#include "recherche_motif.h"
#include "prediction_regions_codantes.h"




File* rechercheCodonsStart(Nucleotide sequence[], int taille, AutomateSuffixes** automates)
{

	if (automates != NULL)
	{
		return rechercheAutomateSuffixes(automates[3], sequence, taille);
	}
	else
	{
		Nucleotide codonStart[] = { A,T,G };
		return boyerMoore(sequence, codonStart, taille, 3);
	}
}



File* rechercheCodonsStop(Nucleotide sequence[], int taille, AutomateSuffixes** automates)
{

	File* codonsStop0 = NULL;
	File* codonsStop1 = NULL;
	File* codonsStop2 = NULL;

	if (automates != NULL)
	{
		codonsStop0 = rechercheAutomateSuffixes(automates[STOP0], sequence, taille);
		codonsStop1 = rechercheAutomateSuffixes(automates[STOP1], sequence, taille);
		codonsStop2 = rechercheAutomateSuffixes(automates[STOP2], sequence, taille);
	}
	else
	{
		Nucleotide codonStop0[] = { T, A, A };
		Nucleotide codonStop1[] = { T, A, G };
		Nucleotide codonStop2[] = { T, G, A };
		codonsStop0 = boyerMoore(sequence, codonStop0, taille, 3);
		codonsStop1 = boyerMoore(sequence, codonStop1, taille, 3);
		codonsStop2 = boyerMoore(sequence, codonStop2, taille, 3);
	}


	File* codonsStop = malloc(sizeof(File));

	if (codonsStop == NULL)
	{
		exit(EXIT_FAILURE);
	}

	codonsStop->dernierElement = NULL;
	codonsStop->premierElement = NULL;


	while (codonsStop0->premierElement != NULL)
	{
		if (codonsStop1->premierElement != NULL)
		{
			if (codonsStop2->premierElement != NULL)
			{
				if (codonsStop0->premierElement->element < codonsStop1->premierElement->element)
				{
					if (codonsStop0->premierElement->element < codonsStop2->premierElement->element)
					{
						enfiler(codonsStop, defiler(codonsStop0));
					}
					else
					{
						enfiler(codonsStop, defiler(codonsStop2));
					}
				}
				else
				{
					if (codonsStop1->premierElement->element < codonsStop2->premierElement->element)
					{
						enfiler(codonsStop, defiler(codonsStop1));
					}
					else
					{
						enfiler(codonsStop, defiler(codonsStop2));
					}
				}
			} // codonStopActuel2 = NULL

			else
			{
				if (codonsStop0->premierElement->element < codonsStop1->premierElement->element)
				{
					enfiler(codonsStop, defiler(codonsStop0));
				}
				else
				{
					enfiler(codonsStop, defiler(codonsStop1));
				}
			}
		} // codonStopActuel1 = NULL

		else
		{
			if (codonsStop2->premierElement != NULL)
			{
				if (codonsStop0->premierElement->element < codonsStop2->premierElement->element)
				{
					enfiler(codonsStop, defiler(codonsStop0));
				}
				else
				{
					enfiler(codonsStop, defiler(codonsStop2));
				}
			} // codonStopActuel2 = NULL
			else
			{
				enfiler(codonsStop, defiler(codonsStop0));
			}
		}
	} // codonStopActuel0 = NULL

	while (codonsStop1->premierElement != NULL)
	{
		if (codonsStop2->premierElement != NULL)
		{
			if (codonsStop1->premierElement->element < codonsStop2->premierElement->element)
			{
				enfiler(codonsStop, defiler(codonsStop1));
			}
			else
			{
				enfiler(codonsStop, defiler(codonsStop2));
			}
		}

		else
		{
			enfiler(codonsStop, defiler(codonsStop1));
		}
	}

	while (codonsStop2->premierElement != NULL)
	{
		enfiler(codonsStop, defiler(codonsStop2));
	}

	return codonsStop;
}






File* predictionORF(Nucleotide sequence[], int taille, AutomateSuffixes** automates)
{
	File* codonsStart = rechercheCodonsStart(sequence, taille, automates);
	File* codonsStop = rechercheCodonsStop(sequence, taille, automates);


	File* orfs = malloc(sizeof(File));
	orfs->dernierElement = NULL;
	orfs->premierElement = NULL;


	int tailleSequenceCodante = 0;



	while (codonsStart->premierElement != NULL && codonsStop->premierElement != NULL)
	{
		// recherche de deux codons stops en phase successifs distants d'au moins TAILLE
		int indiceCodonStop0 = defiler(codonsStop); // premier codon stop
		int indiceCodonStop1 = defiler(codonsStop); // deuxieme codon stop
		int invalide = (indiceCodonStop1 - indiceCodonStop0) % 3 != 0 || (indiceCodonStop1 - indiceCodonStop0) <= TAILLE; // = 1 !
		int fileVide = 0;


		while (invalide)
		{
			if (codonsStop->premierElement == NULL)
			{
				invalide = 0;
				fileVide = 1;
			}
			else if ((indiceCodonStop1 - indiceCodonStop0) % 3 == 0 && (indiceCodonStop1 - indiceCodonStop0) < TAILLE)
			{
				int indiceCodonStopSuivant = defiler(codonsStop);
				indiceCodonStop0 = indiceCodonStop1;
				indiceCodonStop1 = indiceCodonStopSuivant; // il ne doit pas y avoir de codon stop entre les deux !
			}
			else if ((indiceCodonStop1 - indiceCodonStop0) % 3 != 0)
			{
				indiceCodonStop1 = defiler(codonsStop);
			}
			else
			{
				invalide = 0;
			}
		}


		if (!fileVide) // on a trouve deux codons stops ideaux
		{

			if (codonsStart->premierElement != NULL)
			{
				int indiceCodonStart = defiler(codonsStart);

				invalide = indiceCodonStart > indiceCodonStop0 && indiceCodonStart < indiceCodonStop1 && (indiceCodonStart - indiceCodonStop0) % 3 == 0 && indiceCodonStop1 - indiceCodonStart >= TAILLE;
				fileVide = 0;

				while (invalide)
				{
					if (codonsStart->premierElement == NULL)
					{
						invalide = 0;
						fileVide = 1;
					}
					else if (indiceCodonStart > indiceCodonStop1)
					{
						invalide = 0; // le codon start doit etre entre les deux codons stop donc on sort
					}
					else if ((indiceCodonStart - indiceCodonStop0) % 3 != 0 || (indiceCodonStop1 - indiceCodonStart) < TAILLE || indiceCodonStart < indiceCodonStop0)
					{
						indiceCodonStart = defiler(codonsStart);
					}
					else
					{
						invalide = 0;
					}
				}


				// ajout
				if (!fileVide)
				{

					enfiler(orfs, indiceCodonStart);
					enfiler(orfs, indiceCodonStop1);

					tailleSequenceCodante += indiceCodonStop1 - 3 - indiceCodonStart;



					// decalage de codonStartActuel pour que l'indice du codon stop soit superieur a celui du codon start

					int codonStartInvalide = codonsStart->premierElement != NULL;

					while (codonStartInvalide)
					{
						if (codonsStart->premierElement->element > indiceCodonStop1) // le codon start suivant fonctionne
						{
							codonStartInvalide = 0;
						}
						else if (codonsStart->premierElement->suivant != NULL)
						{
							if (codonsStart->premierElement->suivant->element > indiceCodonStop1)
							{
								defiler(codonsStart);
								codonStartInvalide = 0;
							}
							else
							{
								defiler(codonsStart);
							}
						}
						else
						{
							defiler(codonsStart);
							codonStartInvalide = 0;
						}
					}
				}
			}
		}

	}

	// ajout de la taille en tete de file
	
	Element* premierElement = orfs->premierElement;
	Element* nouvelElement = malloc(sizeof(Element));
	nouvelElement->element = tailleSequenceCodante;
	nouvelElement->precedent = NULL;
	nouvelElement->suivant = NULL;

	if (premierElement != NULL)
	{
		premierElement->precedent = nouvelElement;
		nouvelElement->suivant = premierElement;
	}
	else
	{
		orfs->dernierElement = nouvelElement;
	}

	orfs->premierElement = nouvelElement;

	return orfs;
}



Nucleotide* extraireSequencesCodantes(Nucleotide sequence[], int taille, AutomateSuffixes** automates)
{
	File* orfsFile = predictionORF(sequence, taille, automates);

	int tailleSequenceCodante = defiler(orfsFile);

	Nucleotide* orfs = malloc(sizeof(int) * (tailleSequenceCodante + 1));
	int indice = 1;

	if (orfs == NULL)
	{
		exit(EXIT_FAILURE);
	}

	orfs[0] = tailleSequenceCodante;



	while (orfsFile->premierElement != NULL)
	{
		int indiceStartCodon = defiler(orfsFile);
		int indiceStopCodon = defiler(orfsFile);

		for (int i = indiceStartCodon; i < indiceStopCodon - 3; i++)
		{
			orfs[indice] = sequence[i];
			indice++;
		}
	}




	return orfs;
}
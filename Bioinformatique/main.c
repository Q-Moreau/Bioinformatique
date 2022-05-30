/* Quentin MOREAU
* TIPE Bioinformatique
* Reconstruction d'arbres phylogenetiques
*
* fichier source main
*
* 2021-2022
*/



#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "main.h"
#include "fichiers.h"
#include "recherche_motif.h"
#include "prediction_regions_codantes.h"
#include "alignement.h"
#include "construction_arbre_phylogenetique.h"




int main()
{


	char* virusFichiers[30] = { "sequences/Coronaviridae/Middle East respiratory syndrome-related coronavirus.fasta",
		"sequences/Coronaviridae/Human coronavirus 229E.fasta",
		"sequences/Coronaviridae/Bat coronavirus HKU10.fasta",
		"sequences/Coronaviridae/Bat coronavirus CDPHE15.fasta",
		"sequences/Coronaviridae/Human coronavirus NL63.fasta",
		"sequences/Coronaviridae/Human coronavirus HKU1.fasta",
		"sequences/Coronaviridae/Severe acute respiratory syndrome-related coronavirus.fasta",
		"sequences/Roniviridae/Gill-associated virus.fasta",
		"sequences/Roniviridae/Yellow head virus.fasta",
		"sequences/Sphaerolipoviridae/Haloarcula hispanica icosahedral virus 2.fasta",
		"sequences/Sphaerolipoviridae/Haloarcula hispanica virus PH1.fasta",
		"sequences/Sphaerolipoviridae/Haloarcula virus HCIV1.fasta",
		"sequences/Adnaviria/Acidianus filamentous virus 3.fasta",
		"sequences/Adnaviria/Acidianus filamentous virus 6.fasta",
		"sequences/Adnaviria/Acidianus filamentous virus 8.fasta",
		"sequences/Adnaviria/Sulfolobales Beppu filamentous virus 2.fasta",
		"sequences/Adnaviria/Sulfolobus filamentous virus 1.fasta",
		"sequences/Alphasatellitidae/Ageratum yellow vein Singapore alphasatellite.fasta",
		"sequences/Alphasatellitidae/Coconut foliar decay alphasatellite.fasta",
		"sequences/Anelloviridae/Simian torque teno virus 30.fasta",
		"sequences/Anelloviridae/Simian torque teno virus 31.fasta",
		"sequences/Anelloviridae/Simian torque teno virus 32.fasta",
		"sequences/Anelloviridae/Simian torque teno virus 34.fasta",
		"sequences/Anelloviridae/Torque teno virus 10.fasta",
		"sequences/Anelloviridae/Torque teno virus 11.fasta",
		"sequences/Anelloviridae/Torque teno virus 12.fasta",
		"sequences/Anelloviridae/Torque teno virus 13.fasta",
		"sequences/Anelloviridae/Torque teno virus 14.fasta",
		"sequences/Anelloviridae/Torque teno virus 15.fasta",
		"sequences/Anelloviridae/Torque teno virus 16.fasta"
	};

	char* virusNoms[30] = {"Middle East respiratory syndrome-related coronavirus",
		"Human coronavirus 229E",
		"Bat coronavirus HKU10",
		"Bat coronavirus CDPHE15",
		"Human coronavirus NL63",
		"Human coronavirus HKU1",
		"Severe acute respiratory syndrome-related coronavirus",
		"Gill-associated virus",
		"Yellow head virus",
		"Haloarcula hispanica icosahedral virus 2",
		"Haloarcula hispanica virus PH1",
		"Haloarcula virus HCIV1",
		"Acidianus filamentous virus 3",
		"Acidianus filamentous virus 6",
		"Acidianus filamentous virus 8",
		"Sulfolobales Beppu filamentous virus 2",
		"Sulfolobus filamentous virus 1",
		"Ageratum yellow vein Singapore alphasatellite",
		"Coconut foliar decay alphasatellite",
		"Simian torque teno virus 30",
		"Simian torque teno virus 31",
		"Simian torque teno virus 32",
		"Simian torque teno virus 34",
		"Torque teno virus 10",
		"Torque teno virus 11",
		"Torque teno virus 12",
		"Torque teno virus 13",
		"Torque teno virus 14",
		"Torque teno virus 15",
		"Torque teno virus 16"
	};

	arbrePhylogenetique(virusFichiers, virusNoms, 30);







	return 0;
}








Arbre* arbrePhylogenetique(char* nomsFichiers[], char* noms[], int nombreSequence)
{
	Nucleotide** sequencesLues = malloc(sizeof(Nucleotide*) * nombreSequence);
	Nucleotide** sequences = malloc(sizeof(Nucleotide*) * nombreSequence);
	int* tailles = malloc(sizeof(int) * nombreSequence);

	if (sequencesLues == NULL || sequences == NULL || tailles == NULL)
	{
		printf("echec (arbrePhylogenetique)");
		exit(EXIT_FAILURE);
	}

	Nucleotide codonStop0[] = { T, A, A };
	Nucleotide codonStop1[] = { T, A, G };
	Nucleotide codonStop2[] = { T, G, A };
	Nucleotide codonStart[] = { A, T, G };


	AutomateSuffixes* automates[4] = { constructionAutomateSuffixes(codonStop0, 3),
									constructionAutomateSuffixes(codonStop1, 3),
									constructionAutomateSuffixes(codonStop2, 3),
									constructionAutomateSuffixes(codonStart, 3)
	};


	for (int i = 0; i < nombreSequence; i++)
	{
		sequencesLues[i] = lireSequence(nomsFichiers[i]);
		
		Nucleotide* sequenceCodanteLue = extraireSequencesCodantes(sequencesLues[i] + 1, sequencesLues[i][0], automates);
		sequences[i] = sequenceCodanteLue + 1;
		tailles[i] = sequenceCodanteLue[0];

		if (tailles[i] == 0)
		{
			printf("pas de region codante trouvee");
			exit(EXIT_FAILURE);
		}
	}

	int** distances = distancesLevenshtein(sequences, tailles, nombreSequence);
	
	Arbre* arbre = constructionArbreUPGMA(distances, noms, nombreSequence);

	afficherArbre(arbre);

	return NULL;

}















/////////////////////////////////////////// FILES ///////////////////////////////////////////





void enfiler(File* file, int element)
{

	Element* dernierElement = file->dernierElement;
	Element* nouvelElement = malloc(sizeof(Element));

	if (nouvelElement == NULL)
	{
		printf("echec (enfiler)");
		exit(EXIT_FAILURE);
	}

	nouvelElement->element = element;
	nouvelElement->precedent = NULL;
	nouvelElement->suivant = NULL;

	if (dernierElement != NULL)
	{
		dernierElement->suivant = nouvelElement;
		nouvelElement->precedent = dernierElement;
	}
	else
	{
		file->premierElement = nouvelElement;
	}

	file->dernierElement = nouvelElement;
}


int defiler(File* file)
{
	Element* premierElement = file->premierElement;
	if (premierElement == NULL)
	{
		return -1;
	}
	int element = premierElement->element;
	file->premierElement = premierElement->suivant;
	if (premierElement->suivant != NULL)
	{
		premierElement->suivant->precedent = NULL;
	}
	free(premierElement);
	return element;
}



void supprimer(File* file)
{
	while (file->premierElement != NULL)
	{
		defiler(file);
	}
}
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <math.h>

using namespace std;

#define MaxLineLength 10000

char egFileIdx[MaxLineLength], egFileSeq[MaxLineLength], egFileOut[MaxLineLength];
int egKMer, egWindowSize;
char egMode[MaxLineLength];

int KMPSyntax(void)
{
	cout << "K-mer Profiler             version 1.2 (c) Fengfeng Zhou (FengfengZhou@gmail.com)" << endl;
	cout << "                2009-05-10" << endl;
	cout << "" << endl;
	cout << "Syntax:" << endl;
	cout << "     ./KMerProfiler [-index Index.file] -i Input.file.fasta -o Output.KMerVectors.txt -k KMer" << endl;
	cout << "Note:" << endl;
	cout << "	[-index Index.file]           (Optional) The ID of each K-mer will be saved into this file." << endl;
	cout << "	-i Input.file.fasta           The input sequence file in FASTA format." << endl;
	cout << "	-o Output.KMerVectors.txt     Output.KMerVectors.txt." << endl;
	cout << "	-k KMer                       The length of the KMer." << endl;
	cout << "	-w WindowSize                 The length of the window (the last window of length < WindowSize will be omitted)." << endl;
	cout << "Example command line:" << endl;
	cout << "	./KMerProfiler -index Example.index.txt -i Example.fasta -o Example.Average.kmp -k 5 -w 1000" << endl;
	cout << "Some requirements on the input FASTA sequence:" << endl;
	cout << "   >SeqLength|SeqName" << endl;
	cout << "   Sequence" << endl;
}

#define egNUNum 4
char egNU[egNUNum] = {'A', 'C', 'T', 'G'};
int *egNU2Int;  // map all 256 characters to integers
int egNU2IntNumber = 256;

// Keep a KMer and its reverse complement together
int *egFreq; float *egFreqFlag; int egFreqSize = 0; int egFreqNum = 0;
int Weight = 1;

int rvID, rvRCID;
char *rvKMer, *rvKMerRC;


void KMerRC2ID(void)
{
	int tID = 0;
	int weight = 1;
	if( strlen(rvKMerRC)!=egKMer )
	{
		tID = -1;
	}
	else
	{
		for(int i=0;i<egKMer;i++)
		{
			tID = tID + egNU2Int[(int)rvKMerRC[egKMer-i-1]]*weight;
			weight = weight*egNUNum;
		}
	}
	rvRCID = tID;
}

void KMer2ID(void)
{
	int tID = 0;
	int weight = 1;
	if( strlen(rvKMer)!=egKMer )
	{
		tID = -1;
	}
	else
	{
		for(int i=0;i<egKMer;i++)
		{
			tID = tID + egNU2Int[(int)rvKMer[egKMer-i-1]]*weight;
			weight = weight*egNUNum;
		}
	}
	rvID = tID;
}

void ID2KMer(void)
{
	if( rvID>=egFreqSize )
	{
		printf("-Error! KMer id out of range!\n");
		exit(-3);
	}
	
	int tID = rvID;
	int tPos = 0;
	memset(rvKMer, 0, (egKMer+1)*sizeof(char));
	for(int i=0;i<egKMer;i++)
	{
		tPos = tID % egNUNum; tID = (int)((tID-tPos)/egNUNum);
		rvKMer[egKMer-i-1] = egNU[tPos];
	}
}

void KMerRC(void)
{
	memset(rvKMerRC, 0, (egKMer+1)*sizeof(char));
	for(int i=0;i<egKMer;i++)
	{
		if( (rvKMer[i]=='A')||(rvKMer[i]=='a') ) { rvKMerRC[egKMer-i-1] = 'T'; }
		else if( (rvKMer[i]=='T')||(rvKMer[i]=='t') ) { rvKMerRC[egKMer-i-1] = 'A'; }
		else if( (rvKMer[i]=='C')||(rvKMer[i]=='c') ) { rvKMerRC[egKMer-i-1] = 'G'; }
		else if( (rvKMer[i]=='G')||(rvKMer[i]=='g') ) { rvKMerRC[egKMer-i-1] = 'C'; }
		else { rvKMerRC[egKMer-i-1] = rvKMer[i]; }
	}
}


int InitKMP(char *fidx)
{
	if( (egNU2Int=(int*)malloc(egNU2IntNumber*sizeof(int)))==NULL )
	{
		printf("-Error when allocating the memory!\n");
		exit(-2);
	}
	memset(egNU2Int, '\0', egNU2IntNumber*sizeof(int));

	egNU2Int[(int)('A')] = 0;
	egNU2Int[(int)('a')] = 0;
	egNU2Int[(int)('C')] = 1;
	egNU2Int[(int)('c')] = 1;
	egNU2Int[(int)('T')] = 2;
	egNU2Int[(int)('t')] = 2;
	egNU2Int[(int)('G')] = 3;
	egNU2Int[(int)('g')] = 3;

	egFreqSize = (int) pow(4, egKMer);
	Weight = (int)pow(4, egKMer);
	if( egKMer%2==1 )
	{
		egFreqNum = (int)(pow(4, egKMer)/2);
	}
	else
	{
		egFreqNum = (int)( ( pow(4, egKMer)+pow(4, (int)(egKMer/2)) )/2 );
	}

	if( (egFreq=(int*)malloc(egFreqSize*sizeof(int)))==NULL )
	{
		printf("-Error when allocating the memory for the frequency vector!\n");
		exit(-2);
	}
	memset(egFreq, 0, egFreqSize*sizeof(int));

	if( (egFreqFlag=(float*)malloc(egFreqSize*sizeof(float)))==NULL )
	{
		printf("-Error when allocating the memory for the frequency vector!\n");
		exit(-2);
	}
	memset(egFreqFlag, 0, egFreqSize*sizeof(float));

	rvID = 0;
	if( (rvKMer=(char*)malloc((egKMer+1)*sizeof(char)))==NULL )
	{
		printf("-Error when allocating the memory!\n");
		exit(-2);
	}
	memset(rvKMer, 0, (egKMer+1)*sizeof(char));

	if( (rvKMerRC=(char*)malloc((egKMer+1)*sizeof(char)))==NULL )
	{
		printf("-Error when allocating the memory!\n");
		exit(-2);
	}
	memset(rvKMerRC, 0, (egKMer+1)*sizeof(char));

	FILE *fhidx;
	if( strncmp(fidx, "", MaxLineLength)!=0 )
	{
		if( (fhidx=fopen(fidx, "w"))== 0 )
		{
			printf("-Error when saving the index file!\n");
		}
		else
		{
			for(int i=0;i<egFreqSize;i++)
			{
				rvID = i;
				ID2KMer(); KMerRC(); KMerRC2ID();
				if( rvRCID==rvID )
				{
					fprintf(fhidx, "%d	%s\n", rvID, rvKMer);
					egFreqFlag[rvID] = 1.0;
				}
				else if( rvRCID>rvID )
				{
					fprintf(fhidx, "%d	%s	%s\n", rvID, rvKMer, rvKMerRC);
					egFreqFlag[rvID] = 0.5;
					egFreqFlag[rvRCID] = 0;
				}
			}
			fclose(fhidx);
		}
	}

}

void ProcessSequence(char *qSeq)
{
	int klen = strlen(qSeq);
	int kid;
	int ki, kj;
	for(ki=0;ki<klen-egKMer;ki++)
	{
		for(kj=0;kj<egKMer;kj++)
		{
			rvKMer[kj] = qSeq[ki+kj];
		}
		rvKMer[kj] = '\0';
		KMerRC();
		KMer2ID(); KMerRC2ID();
		if( rvID>rvRCID ) { kid = rvRCID; }
		else { kid = rvID; }
		egFreq[kid]++;
	}
}

int KMP(char *fi, char *fo)
{
	FILE *fhi, *fho;

    fhi=fopen(fi, "r");
	if( fhi==NULL )
	{
		printf("-Error when reading the sequence file!\n");
		exit(-4);
	}
	fho=fopen(fo, "w");
	if( fho==NULL )
	{
		printf("-Error when reading the sequence file!\n");
		exit(-4);
	}

	char pSeqLine[MaxLineLength];
	char pSeqName[MaxLineLength];
	int pTotalLength = 0;
	int pSeqLength = 0;
	char *pSequence;

	while( !feof(fhi) )
	{
		fgets(pSeqLine, MaxLineLength, fhi);
		pSeqLength = strlen(pSeqLine);
		if( (pSeqLine[pSeqLength-1]=='\r') || (pSeqLine[pSeqLength-1]=='\n') ) { pSeqLine[pSeqLength-1] = '\0'; }
		if( (pSeqLine[pSeqLength-2]=='\r') || (pSeqLine[pSeqLength-2]=='\n') ) { pSeqLine[pSeqLength-2] = '\0'; }
		pSeqLength = strlen(pSeqLine);

		if( pSeqLine[0]=='>' )
		{
			sscanf(pSeqLine, ">%d|%s", &pTotalLength, pSeqName);
			if( (pSequence=(char*)malloc((pTotalLength+5)*sizeof(char)))==NULL )
			{
				printf("-Error when allocating the memory!\n");
				exit(-2);
			}
			memset(pSequence, 0, (pTotalLength+5)*sizeof(char));
			fgets(pSequence, pTotalLength+5, fhi);
			pSeqLength = strlen(pSequence);
			if( (pSequence[pSeqLength-1]=='\r') || (pSequence[pSeqLength-1]=='\n') ) { pSequence[pSeqLength-1] = '\0'; }
			if( (pSequence[pSeqLength-2]=='\r') || (pSequence[pSeqLength-2]=='\n') ) { pSequence[pSeqLength-2] = '\0'; }
			ProcessSequence(pSequence);
			pTotalLength = strlen(pSequence);
			if( (pTotalLength>0)&&(pSequence[0]!='\0')&&(pSeqName[0]!='\0') )
			{
                fprintf(fho, "%s	%d", pSeqName, Weight);
                for(int j=0;j<egFreqSize;j++)
                {
					if( egFreqFlag[j]!=0 )
					{
                    	fprintf(fho, "	%.5f", (float)((egFreq[j]*1.0*egFreqFlag[j]*Weight)/pTotalLength));
					}
                }
                fprintf(fho, "\n");
			}
			memset(pSeqName, 0, MaxLineLength*sizeof(char));
			memset(egFreq, 0, egFreqSize*sizeof(int));
			free(pSequence);
		}
	}
	fclose(fhi);
	fclose(fho);
}

int InitARGV(int argc, char **argv)
{
	for(int i=1;i<argc-1;i++)
	{
		if( strncmp(argv[i], "-index", MaxLineLength)==0 )
		{
			strncat(egFileIdx, argv[i+1], MaxLineLength);
		}
		else if( strncmp(argv[i], "-i", MaxLineLength)==0 )
		{
			strncat(egFileSeq, argv[i+1], MaxLineLength);
		}
		else if( strncmp(argv[i], "-o", MaxLineLength)==0 )
		{
			strncat(egFileOut, argv[i+1], MaxLineLength);
		}
		else if( strncmp(argv[i], "-k", MaxLineLength)==0 )
		{
			egKMer = atoi(argv[i+1]);
		}
		else if( strncmp(argv[i], "-w", MaxLineLength)==0 )
		{
			egWindowSize = atoi(argv[i+1]);
		}
	}

	if(
		(egFileIdx[0]=='\0')
	 || (egFileSeq[0]=='\0')
	 || (egFileOut[0]=='\0')
	 || ((egKMer<1)||(egKMer>6))
	)
	{
		printf("-Error in parameter range!\n");
		KMPSyntax();
		exit(-1);
	}

	InitKMP(egFileIdx);

}

int main(int argc, char **argv)
{
	InitARGV(argc, argv);
	KMP(egFileSeq, egFileOut);

}



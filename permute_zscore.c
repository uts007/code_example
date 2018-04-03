#include <windows.h>
#include <wininet.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <TCHAR.H>
#include <time.h>
#include <math.h>

#define	Z_MAX          6.0            /* maximum meaningful z value */

#define LOOP	100000
#define FIRST_COL	3
//#define TOTAL_SAMPLE	100

double				/*VAR returns cumulative probability from -oo to z */
poz (double z)		/*VAR normal z value */
{
	double	y, x, w;
	if (z == 0.0)
		x = 0.0;
	else
	{
		y = 0.5 * fabs (z);
		if (y >= (Z_MAX * 0.5))
			x = 1.0;
		else if (y < 1.0)
		{
			w = y*y;
			x = ((((((((0.000124818987 * w
				-0.001075204047) * w +0.005198775019) * w
				-0.019198292004) * w +0.059054035642) * w
				-0.151968751364) * w +0.319152932694) * w
				-0.531923007300) * w +0.797884560593) * y * 2.0;
		}
		else
		{
			y -= 2.0;
			x = (((((((((((((-0.000045255659 * y
				+0.000152529290) * y -0.000019538132) * y
				-0.000676904986) * y +0.001390604284) * y
				-0.000794620820) * y -0.002034254874) * y
				+0.006549791214) * y -0.010557625006) * y
				+0.011630447319) * y -0.009279453341) * y
				+0.005353579108) * y -0.002141268741) * y
				+0.000535310849) * y +0.999936657524;
		}
	}
	return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
}


main( int argc, char *argv[]){
	int s, t, u, i, j, k, n, m, nn, total, count, found, sig_count, GENE_NO, selected[3000];
	double sig, sig2, pvalue[60000], mean, stdev;
	char line[65535], line2[65535], field[800], gene[3000][20];
	char geneid[80], geneid_old[80], value[80], value_old[80];
	FILE *fin, *fin2, *fout, *ftemp;
	int	TOTAL_SAMPLE;
	char *inputfile, *outputfile;

//	fin = fopen("C:/tools/work/data/pathway/pathways_geneid_download.txt", "r");
//	fin = fopen("pathways_geneid_studio_5.txt", "r");
//	fin = fopen("biocarta_Pathways_geneid.txt", "r");
//	fin = fopen("pathway_mouse.txt", "r");
//	fin = fopen("pathway_human.txt", "r");
	fin = fopen("13datasets.txt", "r");
//	fin = fopen("signature/Murat_signature_gene.txt", "r");
	if (!fin)
	{
		printf("1\n");
		return 0;
	}
//	fin2 = fopen("C:/tools/work/data/MPSS_sigature/brain_CompareSet_cl_1_5.txt", "r");
//	fin2 = fopen("C:/tools/work/stemcell/colon/elin_20080821_3.txt", "r");
//	fin2 = fopen("C:/tools/work/stemcell/ES_MPSS/Supplementary_data_3.txt", "r");
//	fin2 = fopen("C:/tools/work/stemcell/brain/REF/Murat_GSE7696_gene.txt", "r");

	inputfile = (char *) malloc(MAX_PATH);
	outputfile = (char *) malloc(MAX_PATH);

	if (argv[1]!=NULL)
		strcpy(inputfile,  argv[1]);
	else
		strcpy(inputfile,  "expressions.txt");

	fin2 = fopen(inputfile, "r");
	if (!fin2)
	{
		printf("type the name of your id list file \r\ndefault is expressions.txt: ");
		scanf("%s", inputfile);

		fin2 = fopen(inputfile, "r");
		if (!fin2)
		{
			printf("\r\ninput file cannot be opened!\r\n\r\n");
			return 0;
		}
	}

	strcpy(outputfile, inputfile);
	i = 0;
	while (outputfile[i]!=0 && outputfile[i]!='.')
		i++;
	outputfile[i] = 0;
	strcat(outputfile, "_sig_13datasets.txt");
	fout = fopen(outputfile, "w");
	if (!fout)
	{
		printf("\r\noutput file cannot be created!\r\n\r\n");
		return 0;
	}
	fputs("dataset\tsource\tgene #\thit #\tsignificance\t", fout);
	fgets(line2, 65530, fin2);
	k = 0;
	for (j=0; j<FIRST_COL-1; j++)
	{
		while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
			k++;
		if (line2[k]=='\t')
			k++;
	}
	i = 0;
	TOTAL_SAMPLE = 0;
	while (line2[k]!='\n' && line2[k]!=0)
	{
		while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
			line2[i++] = line2[k++];
		TOTAL_SAMPLE++;
		if (line2[k]=='\t')
		{
			k++;
			line2[i++] = '\t';
		}
	}
	line2[i++] = '\n';
	line2[i] = 0;
	fputs(line2, fout);

	count = 0;
	while (fgets(line, 65530, fin))
	{
		if (line[0]=='#' || line[0]=='/' || line[0]=='!')
			continue;
		for (s=0; s<TOTAL_SAMPLE; s++)
		{
			ftemp = fopen("temp.txt", "w");
			rewind(fin2);
			fgets(line2, 65530, fin2);
			fgets(line2, 65530, fin2);
			k = 0;
			while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
				geneid_old[k++] = line2[k];
			geneid_old[k] = 0;
			for (j=0; j<FIRST_COL+s-1; j++)
			{
				while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
					k++;
				if (line2[k]=='\t')
					k++;
			}
			j = 0;
			while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
				value_old[j++] = line2[k++];
			value_old[j] = 0;

			while (fgets(line2, 65530, fin2))
			{
				k = 0;
				while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
					geneid[k++] = line2[k];
				geneid[k] = 0;
				for (j=0; j<FIRST_COL+s-1; j++)
				{
					while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
						k++;
					if (line2[k]=='\t')
						k++;
				}
				j = 0;
				while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
					value[j++] = line2[k++];
				value[j] = 0;
				if (strcmp(geneid_old, geneid))
				{
//mask for non-MPSS data
//					if (atof(value_old)>0)
					{
						fputs(geneid_old, ftemp);
						fputs("\t", ftemp);
						fputs(value_old, ftemp);
						fputs("\n", ftemp);
					}
					strcpy(geneid_old, geneid);
					strcpy(value_old, value);
				}
				else if (atof(value)>atof(value_old))
				{
					strcpy(geneid_old, geneid);
					strcpy(value_old, value);
				}
			}
//mask for non-MPSS data
//			if (atof(value_old)>0)
			{
				fputs(geneid_old, ftemp);
				fputs("\t", ftemp);
				fputs(value_old, ftemp);
				fputs("\n", ftemp);
			}
			fclose(ftemp);

			ftemp = fopen("temp.txt", "r");
			GENE_NO = 0;
			mean = 0.0;
			while (fgets(line2, 65530, ftemp))
			{
				GENE_NO++;
				k = 0;
				for (j=0; j<1; j++)
				{
					while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
						k++;
					if (line2[k]=='\t')
						k++;
				}
				j = 0;
				while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
					field[j++] = line2[k++];
				field[j] = 0;
				mean += atof(field);
			}
			mean /= GENE_NO;

			stdev = 0.0;
			rewind(ftemp);
			while (fgets(line2, 65530, ftemp))
			{
				k = 0;
				for (j=0; j<1; j++)
				{
					while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
						k++;
					if (line2[k]=='\t')
						k++;
				}
				j = 0;
				while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
					field[j++] = line2[k++];
				field[j] = 0;
				stdev += (atof(field)-mean)*(atof(field)-mean);
			}
			stdev = sqrt(stdev/(double)GENE_NO);

			m = 0;
			rewind(ftemp);
			while (fgets(line2, 65530, ftemp))
			{
				k = 0;
				for (j=0; j<1; j++)
				{
					while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
						k++;
					if (line2[k]=='\t')
						k++;
				}
				j = 0;
				while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
					field[j++] = line2[k++];
				field[j] = 0;
				if (stdev!=0)
					pvalue[m] = (atof(field) - mean) / stdev;
				else
					pvalue[m] = 0;
				pvalue[m] = poz(pvalue[m]);
				if (pvalue[m]>0.9999)
					pvalue[m] = 0.9999;
				else if (pvalue[m]<0.00001)
					pvalue[m] = 0.00001;
				if (pvalue[m]>0.5)
					pvalue[m] = log10(1.0-pvalue[m])/log10(2.0);
				else if (pvalue[m]<0.5)
					pvalue[m] = -log10(pvalue[m])/log10(2.0);
				else
					pvalue[m] = 0;
				m++;
			}

/*			rewind(ftemp);
			m = 0;
			while (fgets(line2, 65530, ftemp))
			{
				k = 0;
				for (j=0; j<1; j++)
				{
					while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
						k++;
					if (line2[k]=='\t')
						k++;
				}
				j = 0;
				while (line2[k]!='\t' && line2[k]!='\n' && line2[k]!=0)
					field[j++] = line2[k++];
				field[j] = 0;
				pvalue[m] = atof(field) - mean;
				m++;
			}
*/

			sig = 0.0;		//significance in a pathway
			total = 0;		//total gene number in the pathway
			n = 0;			//hits gene number in the pathway

			i = 0;
			k = 0;
			for (j=0; j<2; j++)
			{
				while (line[i]!='\t' && line[i]!='\n' && line[i]!=0)
					field[k++] = line[i++];
				field[k++] = '\t';
				if (line[i]=='\t')
					i++;
			}
			if (k<3)
				continue;
			if (k>0 && field[k-1]=='\t')
				k--;
			field[k] = 0;
			if (s==0)
				fputs(field, fout);

			while (line[i]!='\n' && line[i]!=0)
			{
				j = 0;
				while (line[i]!='\t' && line[i]!='\n' && line[i]!=0)
					field[j++] = line[i++];
				field[j] = 0;
				if (line[i]=='\t')
					i++;
				if (j<1)
					continue;
				found = 0;
				for (k=0; k<total; k++)
				{
					if (strcmp(field, gene[k])==0)
					{
						found = 1;
						break;
					}
				}
				if (found)
				{
					printf("repeat in %d pathway: %s\n", count, field);
					continue;
				}
				strcpy(gene[total], field);
				total++;

				rewind(ftemp);
				nn = 0;
				while (fgets(line2, 65530, ftemp))
				{
					k = 0;
					while (field[k]==line2[k] && field[k]!=0)
						k++;
					if (k>0 && field[k]==0 && line2[k]=='\t')
					{
						sig += pvalue[nn];
						n++;
						break;
					}
					else
						nn++;
				}
			}
			fclose(ftemp);
			if (s==0)
			{
				sprintf(field, "\t%d\t%d\t%.20f", total, n, sig);
				fputs(field, fout);
			}

			sig_count = 0;
			for (k=0; k<LOOP; k++)
			{
				sig2 = 0.0;
				for (i=0; i<n; i++)
				{
					selected[i] = GENE_NO;
					//generate a random number between 0 and (GENE_NO-1)
					j = (int)((((double)rand()) / (double)RAND_MAX) * (double)(GENE_NO-1-i) + 0.5);
					for (t=0; t<i; t++)
					{
						if (selected[t]<j)
							continue;
						else if (selected[t]==j)
							j++;
						else if (selected[t]>j)
						{
							for (u=i; u>t; u--)
								selected[u] = selected[u-1];
							selected[t] = j;
							break;
						}
					}
					if (t==i)
						selected[i] = j;
					sig2 += pvalue[j];
				}
				if (sig2<=sig)
					sig_count++;
			}
			sig = (double)sig_count / (double)LOOP;
			sprintf(field, "\t%.20f", sig);
			fputs(field, fout);
			printf("sample = %d\t%.6f\n", s+1, sig);
		}
		fputs("\n", fout);
		printf("dataset # = %d\n", ++count);
	}

	fcloseall();
	return 0;
}

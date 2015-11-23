/* 
 * File:   main.c
 * Author: veraalva
 *
 * Created on November 13, 2015, 11:45 AM
 */

#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "btree.h"
#include "main.h"

char *program_name;

void print_usage(FILE *stream, int exit_code) {
    fprintf(stream, "\n********************************************************************************\n");
    fprintf(stream, "\nUsage: %s \n", program_name);
    fprintf(stream, "\n\n%s options:\n\n", program_name);
    fprintf(stream, "-v,   --verbose                     Print info\n");
    fprintf(stream, "-h,   --help                        Display this usage information.\n");
    fprintf(stream, "-d,   --dir                         Directory with the ENT and OUT files\n");
    fprintf(stream, "-o,   --output                      Output file\n");
    fprintf(stream, "-g,   --gene                        Gene average output file\n");
    fprintf(stream, "-s,   --sum                         Sum file for checking\n");
    fprintf(stream, "-r,   --rep                         File for printing TPM values for gene for replicates\n");
    fprintf(stream, "********************************************************************************\n");
    fprintf(stream, "\n            Roberto Vera Alvarez (e-mail: r78v10a07@gmail.com)\n\n");
    fprintf(stream, "********************************************************************************\n");
    exit(exit_code);
}

int main(int argc, char** argv) {
    int next_option, verbose;
    const char* const short_options = "vhd:o:s:g:r:";
    FILE *s;
    char *line = NULL;
    size_t len = 0;
    ssize_t read = 0;
    char **fields = NULL;
    size_t fieldsSize = 0;
    char **samples = NULL;
    int samplesLen = 0;
    BtreeNode_t *genes = NULL;
    void **genesArray = NULL;
    int genesLen = 0;
    struct dirent **namelist;
    int i, n, l, k;
    BtreeRecord_t *rec = NULL;
    MGene_l gene;
    double cllSum, nbcSum, tmpValue, tmp1Value;
    int cllLen, nbcLen;
    int cllLenTMP, nbcLenTMP;
    FILE *outFile = NULL;
    FILE *geneFile = NULL;
    FILE *sumFile = NULL;
    FILE *repFile = NULL;
    char *dir = NULL;
    int s1Len;
    int s2Len;

    program_name = argv[0];

    const struct option long_options[] = {
        { "help", 0, NULL, 'h'},
        { "verbose", 0, NULL, 'v'},
        { "output", 1, NULL, 'o'},
        { "dir", 1, NULL, 'd'},
        { "gene", 1, NULL, 'g'},
        { "sum", 1, NULL, 's'},
        { "rep", 1, NULL, 'r'},
        { NULL, 0, NULL, 0} /* Required at end of array.  */
    };

    verbose = 0;
    do {
        next_option = getopt_long(argc, argv, short_options, long_options, NULL);

        switch (next_option) {
            case 'h':
                print_usage(stdout, 0);

            case 'v':
                verbose = 1;
                break;

            case 'd':
                dir = strdup(optarg);
                break;

            case 'o':
                outFile = checkPointerError(fopen(optarg, "w"), "Can't open OUTPUT file", __FILE__, __LINE__, -1);
                break;

            case 'g':
                geneFile = checkPointerError(fopen(optarg, "w"), "Can't open GENE file", __FILE__, __LINE__, -1);
                break;

            case 's':
                sumFile = checkPointerError(fopen(optarg, "w"), "Can't open SUM file", __FILE__, __LINE__, -1);
                break;

            case 'r':
                repFile = checkPointerError(fopen(optarg, "w"), "Can't open replicate file", __FILE__, __LINE__, -1);
                break;
        }
    } while (next_option != -1);
    if (!dir || !outFile || !sumFile) {
        print_usage(stderr, -1);
    }

    fprintf(geneFile, "GeneId\tTranscriptId\tNBC_Gene_TPM\tNBC_Exon_TPM\tNBC_Intron_TPM\tNBC_Intron_Exon_Ratio\tNBC_LOG2_Intron_Exon_Ratio\tCLL_Gene_TPM\tCLL_Exon_TPM\tCLL_Intron_TPM\tCLL_Intron_Exon_Ratio\tCLL_LOG2_Intron_Exon_Ratio\tHausdorff_distance\tCLL_NBC_Intron_Ratio\tLOG2_CLL_NBC_Intron_Ratio\tCLL_NBC_Exon_Ratio\tLOG2_CLL_NBC_Exon_Ratio\n");
    fprintf(outFile, "GeneId\tTranscriptId\tIntron_Exon_Number\tType");
    fprintf(repFile, "GeneId\tTranscriptId");

    cllLen = nbcLen = samplesLen = 0;
    n = scandir(dir, &namelist, 0, alphasort);
    if (n < 0)
        perror("scandir");
    else {
        for (i = 0; i < n; i++) {
            if (strbcmp(namelist[i]->d_name, ".ent") == 0) {
                s = fopen(namelist[i]->d_name, "r");
                fseeko(s, 0, SEEK_END);
                if (ftello(s) != 0) {
                    samples = reallocate(samples, sizeof (samples) * (samplesLen + 1), __FILE__, __LINE__);
                    samples[samplesLen] = strdup(namelist[i]->d_name);
                    if (strncmp(samples[samplesLen], "CLL", 3) == 0) cllLen++;
                    if (strncmp(samples[samplesLen], "NBC", 3) == 0) nbcLen++;
                    samplesLen++;
                }
                fclose(s);
            }
            free(namelist[i]);
        }
    }
    free(namelist);

    for (i = 0; i < samplesLen; i++) {
        s = checkPointerError(fopen(samples[i], "r"), "Can't open sample ENT file", __FILE__, __LINE__, -1);
        *(strstr(samples[i], ".ent")) = '\0';
        fprintf(outFile, "\t%s\t%s_Reads", samples[i], samples[i]);
        fprintf(repFile, "\t%s_Exon\t%s_Exon_Reads\t%s_Intron\t%s_Intron_Reads", samples[i], samples[i], samples[i], samples[i]);
        if (verbose) printf("Parsing ENT file: %s\n", samples[i]);

        while ((read = getline(&line, &len, s)) != -1) {
            if (strncmp(line, "Gene_Id", 7) != 0) {
                fieldsSize = splitString(&fields, line, "\t");
                if (fieldsSize == 8) {
                    gene = NewMGene(fields);
                    rec = BTreeFind(genes, gene, GeneKeyCMP);
                    if (rec == NULL) {
                        gene->entities = reallocate(gene->entities, sizeof (MEntity_t) * (gene->entitiesLen + 1), __FILE__, __LINE__);
                        gene->entities[gene->entitiesLen].length = atoi(fields[5]);
                        gene->entities[gene->entitiesLen].type = strdup(fields[3]);
                        gene->NBCIntronExonRatio = allocate(sizeof (double) * nbcLen, __FILE__, __LINE__);
                        for (k = 0; k < nbcLen; k++) {
                            gene->NBCIntronExonRatio[k] = NAN;
                        }
                        gene->CLLIntronExonRatio = allocate(sizeof (double) * cllLen, __FILE__, __LINE__);
                        for (k = 0; k < cllLen; k++) {
                            gene->CLLIntronExonRatio[k] = NAN;
                        }

                        gene->repExonTPM = allocate(sizeof (double) * samplesLen, __FILE__, __LINE__);
                        for (k = 0; k < samplesLen; k++) {
                            gene->repExonTPM[k] = NAN;
                        }
                        gene->repExonCount = allocate(sizeof (int) * samplesLen, __FILE__, __LINE__);
                        for (k = 0; k < samplesLen; k++) {
                            gene->repExonCount[k] = NAN;
                        }

                        gene->repIntronTPM = allocate(sizeof (double) * samplesLen, __FILE__, __LINE__);
                        for (k = 0; k < samplesLen; k++) {
                            gene->repIntronTPM[k] = NAN;
                        }

                        gene->repIntronCount = allocate(sizeof (int) * samplesLen, __FILE__, __LINE__);
                        for (k = 0; k < samplesLen; k++) {
                            gene->repIntronCount[k] = NAN;
                        }

                        gene->entities[gene->entitiesLen].samples = allocate(sizeof (MEntitySample_t) * (samplesLen), __FILE__, __LINE__);
                        for (k = 0; k < samplesLen; k++) {
                            gene->entities[gene->entitiesLen].samples[k].sample = samples[k];
                            gene->entities[gene->entitiesLen].samples[k].count = -1;
                            gene->entities[gene->entitiesLen].samples[k].TPM = INFINITY;
                        }
                        gene->entities[gene->entitiesLen].samples[i].count = atoi(fields[6]);
                        gene->entities[gene->entitiesLen].samples[i].TPM = strtod(fields[7], NULL);
                        gene->entitiesLen++;

                        genes = BtreeInsert(genes, gene, gene, GeneKeyCMP);
                    } else {
                        FreeMGene(gene);
                        gene = rec->value;
                        l = atoi(fields[4]) - 1;
                        if (l >= gene->entitiesLen) {
                            gene->entities = reallocate(gene->entities, sizeof (MEntity_t) * (gene->entitiesLen + 1), __FILE__, __LINE__);
                            gene->entities[l].length = atoi(fields[5]);
                            gene->entities[l].type = strdup(fields[3]);
                            gene->entities[l].samples = allocate(sizeof (MEntitySample_t) * (samplesLen), __FILE__, __LINE__);
                            for (k = 0; k < samplesLen; k++) {
                                gene->entities[gene->entitiesLen].samples[k].sample = samples[k];
                                gene->entities[gene->entitiesLen].samples[k].count = -1;
                                gene->entities[gene->entitiesLen].samples[k].TPM = INFINITY;
                            }
                            gene->entities[l].samples[i].count = atoi(fields[6]);
                            gene->entities[l].samples[i].TPM = strtod(fields[7], NULL);
                            gene->entitiesLen++;
                        } else {
                            gene->entities[l].samples[i].sample = samples[i];
                            gene->entities[l].samples[i].count = atoi(fields[6]);
                            gene->entities[l].samples[i].TPM = strtod(fields[7], NULL);
                        }
                    }
                } else {
                    fprintf(stderr, "\n\n%s\n\n", line);
                    printLog(stderr, "Bad ENT format", __FILE__, __LINE__, -1);
                }
                freeArrayofPointers((void **) fields, fieldsSize);
            }
        }
        fclose(s);
    }

    fprintf(outFile, "\tCLL_Mean\tNBC_Mean\tCLL_NBC_Mean_Diff\n");
    fprintf(repFile, "\n");

    s1Len = s2Len = 0;
    for (i = 0; i < samplesLen; i++) {
        strcat(samples[i], ".out");
        s = checkPointerError(fopen(samples[i], "r"), "Can't open sample OUT file", __FILE__, __LINE__, -1);
        *(strstr(samples[i], ".out")) = '\0';
        if (verbose) printf("Parsing OUT file: %s\n", samples[i]);

        while ((read = getline(&line, &len, s)) != -1) {
            if (strncmp(line, "Gene_Id", 7) != 0) {
                fieldsSize = splitString(&fields, line, "\t");
                if (fieldsSize == 12) {
                    gene = NewMGene(fields);
                    rec = BTreeFind(genes, gene, GeneKeyCMP);
                    if (rec == NULL) {
                        fprintf(stderr, "\n\n%s\n\n", line);
                        printLog(stderr, "Wrong GENE name", __FILE__, __LINE__, -1);
                    }
                    FreeMGene(gene);
                    gene = rec->value;
                    if (gene->length == -1) {
                        gene->length = atoi(fields[3]);
                    } else if (gene->length != atoi(fields[3])) {
                        fprintf(stderr, "\n\n%s\n\n", line);
                        printLog(stderr, "Wrong GENE length", __FILE__, __LINE__, -1);
                    }
                    tmpValue = strtod(fields[5], NULL);
                    if (!isnan(tmpValue) && !isinf(tmpValue)) {
                        if (strncmp(samples[i], "CLL", 3) == 0) {
                            gene->CLLTPM += tmpValue;
                            gene->CLLTMPCount++;
                        } else {
                            gene->NBCTPM += tmpValue;
                            gene->NBCTMPCount++;
                        }
                    }

                    tmpValue = strtod(fields[7], NULL);
                    if (!isnan(tmpValue) && !isinf(tmpValue)) {
                        gene->repExonCount[i] = tmpValue;
                    }

                    tmpValue = strtod(fields[8], NULL);
                    if (!isnan(tmpValue) && !isinf(tmpValue)) {
                        gene->repExonTPM[i] = tmpValue;
                        if (strncmp(samples[i], "CLL", 3) == 0) {
                            gene->CLLTPMExon += tmpValue;
                            gene->CLLexonCount++;
                        } else {
                            gene->NBCTPMExon += tmpValue;
                            gene->NBCexonCount++;
                        }
                    }

                    tmpValue = strtod(fields[10], NULL);
                    if (!isnan(tmpValue) && !isinf(tmpValue)) {
                        gene->repIntronCount[i] = tmpValue;
                    }

                    tmp1Value = strtod(fields[11], NULL);
                    if (!isnan(tmp1Value) && !isinf(tmp1Value)) {
                        gene->repIntronTPM[i] = tmp1Value;
                        if (strncmp(samples[i], "CLL", 3) == 0) {
                            gene->CLLTPMIntron += tmp1Value;
                            gene->CLLintronCount++;

                            gene->CLLIntronExonRatio[s1Len] = tmp1Value;
                        } else {
                            gene->NBCTPMIntron += tmp1Value;
                            gene->NBCintronCount++;

                            gene->NBCIntronExonRatio[s2Len] = tmp1Value;
                        }
                    }
                } else {
                    fprintf(stderr, "\n\n%s\n\n", line);
                    printLog(stderr, "Bad OUT format", __FILE__, __LINE__, -1);
                }
                freeArrayofPointers((void **) fields, fieldsSize);
            }
        }
        fclose(s);
        if (strncmp(samples[i], "CLL", 3) == 0) s1Len++;
        else s2Len++;
    }

    BtreeRecordsToArray(&genesArray, &genesLen, genes);
    for (i = 0; i < genesLen; i++) {
        nbcSum = ((MGene_l) genesArray[i])->NBCTPMExon / ((MGene_l) genesArray[i])->NBCexonCount;
        tmpValue = ((MGene_l) genesArray[i])->NBCTPMIntron / ((MGene_l) genesArray[i])->NBCintronCount;
        fprintf(geneFile, "%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f",
                ((MGene_l) genesArray[i])->geneId,
                ((MGene_l) genesArray[i])->transcriptId,
                ((MGene_l) genesArray[i])->NBCTPM / ((MGene_l) genesArray[i])->NBCTMPCount,
                nbcSum,
                tmpValue,
                tmpValue / nbcSum,
                log2(tmpValue / nbcSum));
        cllSum = ((MGene_l) genesArray[i])->CLLTPMExon / ((MGene_l) genesArray[i])->CLLexonCount;
        tmp1Value = ((MGene_l) genesArray[i])->CLLTPMIntron / ((MGene_l) genesArray[i])->CLLintronCount;
        fprintf(geneFile, "\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                ((MGene_l) genesArray[i])->CLLTPM / ((MGene_l) genesArray[i])->CLLTMPCount,
                cllSum,
                tmpValue,
                tmpValue / cllSum,
                log2(tmpValue / cllSum),
                HausdorffDistance(((MGene_l) genesArray[i])->CLLIntronExonRatio, cllLen, ((MGene_l) genesArray[i])->NBCIntronExonRatio, nbcLen),
                tmp1Value / tmpValue,
                log2(tmp1Value / tmpValue),
                cllSum / nbcSum,
                log2(cllSum / nbcSum));

        for (l = 0; l < ((MGene_l) genesArray[i])->entitiesLen; l++) {
            fprintf(outFile, "%s\t%s\t%d", ((MGene_l) genesArray[i])->geneId, ((MGene_l) genesArray[i])->transcriptId, l + 1);
            cllSum = nbcSum = 0.0;
            cllLenTMP = nbcLenTMP = 0;
            for (k = 0; k < samplesLen; k++) {
                if (strncmp(samples[k], "CLL", 3) == 0 && ((MGene_l) genesArray[i])->entities[l].samples[k].TPM != INFINITY) {
                    cllLenTMP++;
                    cllSum += ((MGene_l) genesArray[i])->entities[l].samples[k].TPM;
                }
                if (strncmp(samples[k], "NBC", 3) == 0 && ((MGene_l) genesArray[i])->entities[l].samples[k].TPM != INFINITY) {
                    nbcLenTMP++;
                    nbcSum += ((MGene_l) genesArray[i])->entities[l].samples[k].TPM;
                }
                if (k == 0) {
                    fprintf(outFile, "\t%s", ((MGene_l) genesArray[i])->entities[l].type);
                }
                fprintf(outFile, "\t%.4f\t%d", ((MGene_l) genesArray[i])->entities[l].samples[k].TPM, ((MGene_l) genesArray[i])->entities[l].samples[k].count);

            }
            cllSum = cllSum / cllLenTMP;
            nbcSum = nbcSum / nbcLenTMP;
            fprintf(outFile, "\t%.4f\t%.4f\t", cllSum, nbcSum);
            if (fabs(cllSum) <= 0.000001 && fabs(nbcSum) <= 0.000001) {
                fprintf(outFile, "nan\n");
            } else {
                fprintf(outFile, "%.4f\n", cllSum - nbcSum);
            }
        }

        fprintf(repFile, "%s\t%s", ((MGene_l) genesArray[i])->geneId, ((MGene_l) genesArray[i])->transcriptId);
        for (k = 0; k < samplesLen; k++) {
            fprintf(repFile, "\t%.4f\t%d\t%.4f\t%d", ((MGene_l) genesArray[i])->repExonTPM[k], ((MGene_l) genesArray[i])->repExonCount[k],
                    ((MGene_l) genesArray[i])->repIntronTPM[k], ((MGene_l) genesArray[i])->repIntronCount[k]);
        }
        fprintf(repFile, "\n");
    }

    for (k = 0; k < samplesLen; k++) {
        cllSum = 0.0;
        for (i = 0; i < genesLen; i++) {
            for (l = 0; l < ((MGene_l) genesArray[i])->entitiesLen; l++) {
                if (((MGene_l) genesArray[i])->entities[l].samples[k].TPM != INFINITY) {
                    cllSum += ((MGene_l) genesArray[i])->entities[l].samples[k].TPM;
                }
            }
        }
        fprintf(sumFile, "%s\t%.4f\n", samples[k], cllSum);
    }

    if (dir) free(dir);
    if (genesArray) free(genesArray);
    BTreeFree(genes, FreeMGene, NULL);
    freeArrayofPointers((void **) samples, samplesLen);
    if (line) free(line);
    fclose(outFile);
    fclose(geneFile);
    fclose(sumFile);
    fclose(repFile);
    return (EXIT_SUCCESS);
}


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
    fprintf(stream, "-d,   --dir                         Directory with the ENT files\n");
    fprintf(stream, "-o,   --output                      Output file\n");
    fprintf(stream, "********************************************************************************\n");
    fprintf(stream, "\n            Roberto Vera Alvarez (e-mail: r78v10a07@gmail.com)\n\n");
    fprintf(stream, "********************************************************************************\n");
    exit(exit_code);
}

int main(int argc, char** argv) {
    int next_option, verbose;
    const char* const short_options = "vhd:o:";
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
    int i, n, l, k, entLen;
    BtreeRecord_t *rec = NULL;
    MGene_l gene;
    double cllSum, nbcSum;
    int cllLen, nbcLen;
    FILE *outFile = NULL;
    char *dir = NULL;
    
    program_name = argv[0];
    
    const struct option long_options[] = {
        { "help", 0, NULL, 'h'},
        { "verbose", 0, NULL, 'v'},
        { "output", 1, NULL, 'o'},
        { "dir", 1, NULL, 'd'},
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
        }
    } while (next_option != -1);
    if (!dir || ! outFile){
        print_usage(stderr, -1);
    }

    n = scandir(dir, &namelist, 0, alphasort);
    if (n < 0)
        perror("scandir");
    else {
        entLen = 0;
        for (i = 0; i < n; i++) {
            if (strbcmp(namelist[i]->d_name, ".ent") == 0) entLen++;
            free(namelist[i]);
        }
    }
    free(namelist);

    n = scandir(dir, &namelist, 0, alphasort);
    if (n < 0)
        perror("scandir");
    else {
        for (i = 0; i < n; i++) {
            if (strbcmp(namelist[i]->d_name, ".ent") == 0) {
                s = fopen(namelist[i]->d_name, "r");
                fseeko(s, 0, SEEK_END);
                if (ftello(s) != 0) {
                    fseeko(s, 0, SEEK_SET);
                    samples = reallocate(samples, sizeof (char **) * (samplesLen + 1), __FILE__, __LINE__);
                    samples[samplesLen] = strdup(namelist[i]->d_name);
                    *(strstr(samples[samplesLen], ".ent")) = '\0';
                    if (verbose) printf("Parsing sample file: %s\n", namelist[i]->d_name);

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
                                    gene->entities[gene->entitiesLen].samples = allocate(sizeof (MEntitySample_t) * (entLen), __FILE__, __LINE__);
                                    for (k = 0; k < entLen; k++) {
                                        gene->entities[gene->entitiesLen].samples[k].sample = samples[k];
                                        gene->entities[gene->entitiesLen].samples[k].count = -1;
                                        gene->entities[gene->entitiesLen].samples[k].TPM = INFINITY;
                                    }
                                    gene->entities[gene->entitiesLen].samples[samplesLen].count = atoi(fields[6]);
                                    gene->entities[gene->entitiesLen].samples[samplesLen].TPM = strtod(fields[7], NULL);
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
                                        gene->entities[l].samples = allocate(sizeof (MEntitySample_t) * (entLen), __FILE__, __LINE__);
                                        for (k = 0; k < entLen; k++) {
                                            gene->entities[gene->entitiesLen].samples[k].sample = samples[k];
                                            gene->entities[gene->entitiesLen].samples[k].count = -1;
                                            gene->entities[gene->entitiesLen].samples[k].TPM = INFINITY;
                                        }
                                        gene->entities[l].samples[samplesLen].count = atoi(fields[6]);
                                        gene->entities[l].samples[samplesLen].TPM = strtod(fields[7], NULL);
                                        gene->entitiesLen++;
                                    } else {
                                        gene->entities[l].samples[samplesLen].sample = samples[samplesLen];
                                        gene->entities[l].samples[samplesLen].count = atoi(fields[6]);
                                        gene->entities[l].samples[samplesLen].TPM = strtod(fields[7], NULL);
                                    }
                                }
                            } else {
                                fprintf(stderr, "\n\n%s\n\n", line);
                                printLog(stderr, "Bad ENT format", __FILE__, __LINE__, -1);
                            }
                            freeArrayofPointers((void **) fields, fieldsSize);
                        }
                    }
                    samplesLen++;
                }

                fclose(s);
            }
            free(namelist[i]);
        }
    }
    free(namelist);

    fprintf(outFile, "GeneId\tTranscriptId\tIntron_Exon_Number");
    for (k = 0; k < samplesLen; k++) {
        fprintf(outFile, "\t%s\t%s_Reads", samples[k], samples[k]);
        if (strncmp(samples[k], "CLL", 3) == 0) cllLen++;
        if (strncmp(samples[k], "NBC", 3) == 0) nbcLen++;
    }
    fprintf(outFile, "\tCLL_Mean\tNBC_Mean\n");

    BtreeRecordsToArray(&genesArray, &genesLen, genes);
    for (i = 0; i < genesLen; i++) {
        for (l = 0; l < ((MGene_l) genesArray[i])->entitiesLen; l++) {
            fprintf(outFile, "%s\t%s\t%d", ((MGene_l) genesArray[i])->geneId, ((MGene_l) genesArray[i])->transcriptId, l + 1);            
            cllSum = nbcSum = 0.0;
            cllLen = nbcLen = 0;
            for (k = 0; k < samplesLen; k++) {
                if (strncmp(samples[k], "CLL", 3) == 0 && ((MGene_l) genesArray[i])->entities[l].samples[k].TPM != INFINITY) {
                    cllLen++;
                    cllSum += ((MGene_l) genesArray[i])->entities[l].samples[k].TPM;
                }
                if (strncmp(samples[k], "NBC", 3) == 0 && ((MGene_l) genesArray[i])->entities[l].samples[k].TPM != INFINITY) {
                    nbcLen++;
                    nbcSum += ((MGene_l) genesArray[i])->entities[l].samples[k].TPM;
                }
                fprintf(outFile, "\t%.4f\t%d", ((MGene_l) genesArray[i])->entities[l].samples[k].TPM, ((MGene_l) genesArray[i])->entities[l].samples[k].count);
            }
            fprintf(outFile, "\t%.4f\t%.4f\n",cllSum/cllLen, nbcSum/nbcLen);
        }
    }

    free(dir);
    free(genesArray);
    BTreeFree(genes, FreeMGene, NULL);
    freeArrayofPointers((void **) samples, samplesLen);
    if (line) free(line);
    fclose(outFile);
    return (EXIT_SUCCESS);
}


#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/unistd.h>
#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "btree.h"
#include "chrentity.h"
#include "gene.h"
#include "chromosome.h"

/********************************************
 * Chromosome Functions
 ********************************************/

int FindSortedIndex(void *self, int key) {
    int i = 0;
    BtreeNode_t * c = ((Chromosome_l) self)->index;
    if (c == NULL) return -1;
    while (!c->is_leaf) {
        i = 0;
        while (i < c->num_keys) {
            if (key >= ((Gene_l) c->keys[i])->pTo) i++; // key >= c->keys[i]
            else break;
        }
        c = (BtreeNode_t *) c->pointers[i];
    }
    return *((int *) ((BtreeRecord_t *) c->pointers[0])->value);
}

Chromosome_l NewChromosome(char *name) {
    Chromosome_l chr = allocate(sizeof (Chromosome_t), __FILE__, __LINE__);
    chr->name = strdup(name);
    chr->genesTree = NULL;
    chr->genesSorted = NULL;
    chr->genesLen = 0;
    chr->index = NULL;
    chr->chromosomeLen = 0;

    chr->FindSortedIndex = &FindSortedIndex;

    return chr;
}

void FreeChromosome(void *self) {
    _CHECK_SELF_P(self);
    if (((Chromosome_l) self)->name) free(((Chromosome_l) self)->name);
    if (((Chromosome_l) self)->genesSorted) free(((Chromosome_l) self)->genesSorted);
    if (((Chromosome_l) self)->genesTree) BTreeFree(((Chromosome_l) self)->genesTree, FreeGenes, NULL);
    if (((Chromosome_l) self)->index) BTreeFree(((Chromosome_l) self)->index, free, NULL);
    free(((Chromosome_l) self));
}

/********************************************
 * Comparison Functions
 ********************************************/

int ChromosomeCMP(void *key1, void *key2) {
    Chromosome_l chrKey1 = (Chromosome_l) key1;
    Chromosome_l chrKey2 = (Chromosome_l) key2;

    return strcmp(chrKey1->name, chrKey2->name);
}

int GeneSortCMP(const void *g1, const void *g2) {
    GeneSorted_t *gene1 = (GeneSorted_t *) g1;
    GeneSorted_t *gene2 = (GeneSorted_t *) g2;

    return gene1->gene->pFrom - gene2->gene->pFrom;
}

int GeneKeyFromCMP(void *key1, void *key2) {
    Gene_l geneKey1 = (Gene_l) key1;
    Gene_l geneKey2 = (Gene_l) key2;

    return geneKey1->pFrom - geneKey2->pFrom;
}

void createGenesSortedIndex(Chromosome_l chr) {
    int i, *index;
    for (i = 0; i < chr->genesLen; i++) {
        index = allocate(sizeof (int), __FILE__, __LINE__);
        *index = i;
        if (BTreeFind(chr->index, chr->genesSorted[i].gene, GeneKeyFromCMP) == NULL) {
            chr->index = BtreeInsert(chr->index, chr->genesSorted[i].gene, index, GeneKeyFromCMP);
        } else {
            free(index);
        }
    }
}

void ChromosomeParse(void *self, FILE *gftFile) {
    char *line = NULL;
    size_t len = 0;
    ssize_t read = 0;
    char **fields = NULL;
    size_t fieldsSize = 0;
    Chromosome_l chr = NULL;
    Gene_l gene = NULL;
    EntityType_l entityType = NULL;
    ChrEntity_l ent, intron;
    BtreeRecord_t *rec = NULL;

    while ((read = getline(&line, &len, gftFile)) != -1) {
        fieldsSize = splitString(&fields, line, "\t");
        if (strchr(fields[0], '_') == NULL && strcmp(fields[2], "exon") == 0) {
            if (!(chr && strncmp(chr->name, fields[0], (strlen(chr->name) > strlen(fields[0]) ? strlen(chr->name) : strlen(fields[0]))) == 0)) {
                if (chr) {
                    BtreeGenesToArray(&chr->genesSorted, &chr->genesLen, chr->genesTree);
                    qsort(chr->genesSorted, chr->genesLen, sizeof (GeneSorted_t), GeneSortCMP);
                    chr->chromosomeLen = chr->genesSorted[chr->genesLen - 1].gene->pTo;
                    createGenesSortedIndex(chr);
                }
                chr = NewChromosome(fields[0]);
                ((Chromosome_f *) self)->chromosomes = BtreeInsert(((Chromosome_f *) self)->chromosomes, chr, chr, ChromosomeCMP);
            }
            entityType = ((Chromosome_f *) self)->getEntityTypeByName(self, fields[2]);
            if (!entityType) {
                ((Chromosome_f *) self)->entityTypes = reallocate(((Chromosome_f *) self)->entityTypes, sizeof (EntityType_l) * (((Chromosome_f *) self)->entityTypesLen + 1), __FILE__, __LINE__);
                ((Chromosome_f *) self)->entityTypes[((Chromosome_f *) self)->entityTypesLen] = NewEntityType(fields[2]);
                entityType = ((Chromosome_f *) self)->entityTypes[((Chromosome_f *) self)->entityTypesLen];
                ((Chromosome_f *) self)->entityTypesLen++;
            }
            fields[8][strlen(fields[8]) - 1] = '\0';
            gene = NewGene(fields[8]);
            rec = BTreeFind(chr->genesTree, gene, GeneKeyCMP);
            if (!rec) {
                chr->genesTree = BtreeInsert(chr->genesTree, gene, gene, GeneKeyCMP);
            } else {
                FreeGenes(gene);
                gene = rec->value;
            }

            ent = NewChrEntity(entityType, fields);

            if (gene->chrEntityLen != 0) {
                if (ent->pFrom - gene->chrEntity[gene->chrEntityLen - 1]->pTo > 0) {
                    gene->chrEntity = reallocate(gene->chrEntity, sizeof (ChrEntity_l) * (gene->chrEntityLen + 1), __FILE__, __LINE__);
                    gene->chrEntity[gene->chrEntityLen] = NewChrEntity(((Chromosome_f *) self)->getEntityTypeByName(self, "intron"), NULL);
                    gene->chrEntity[gene->chrEntityLen]->pFrom = gene->chrEntity[gene->chrEntityLen - 1]->pTo + 1;
                    gene->chrEntity[gene->chrEntityLen]->pTo = ent->pFrom - 1;
                    gene->chrEntity[gene->chrEntityLen]->length = gene->chrEntity[gene->chrEntityLen]->pTo - gene->chrEntity[gene->chrEntityLen]->pFrom + 1;
                    gene->chrEntityLen++;
                }
            }

            gene->chrEntity = reallocate(gene->chrEntity, sizeof (ChrEntity_l) * (gene->chrEntityLen + 1), __FILE__, __LINE__);
            gene->chrEntity[gene->chrEntityLen] = ent;
            if (gene->pFrom == -1 || gene->chrEntity[gene->chrEntityLen]->pFrom < gene->pFrom) {
                gene->pFrom = gene->chrEntity[gene->chrEntityLen]->pFrom;
                gene->length = gene->pTo - gene->pFrom;
            }
            if (gene->pTo == -1 || gene->chrEntity[gene->chrEntityLen]->pTo > gene->pTo) {
                gene->pTo = gene->chrEntity[gene->chrEntityLen]->pTo;
                gene->length = gene->pTo - gene->pFrom;
            }
            gene->chrEntityLen++;
        }
        freeArrayofPointers((void **) fields, fieldsSize);
    }
    if (chr) {
        BtreeGenesToArray(&chr->genesSorted, &chr->genesLen, chr->genesTree);
        qsort(chr->genesSorted, chr->genesLen, sizeof (GeneSorted_t), GeneSortCMP);
        chr->chromosomeLen = chr->genesSorted[chr->genesLen - 1].gene->pTo;
        createGenesSortedIndex(chr);
    }

    if (line) free(line);
}

/**
 * Get an Chromosome structure from a name
 * 
 * @param self a chromosome factory
 * @param name Chromosome name
 * @return an Chromosome structure or NULL
 */
Chromosome_l getChromosomeByName(void *self, char *name) {
    _CHECK_SELF_P(self);
    Chromosome_l chr = NewChromosome(name);
    BtreeRecord_t *res = BTreeFind(((Chromosome_f *) self)->chromosomes, chr, ChromosomeCMP);
    FreeChromosome(chr);
    if (res) return res->value;
    return NULL;
}

/**
 * Get an EntityType structure from a name
 * 
 * @param self a chromosome factory
 * @param name EntityType name
 * @return an EntityType structure or NULL
 */
EntityType_l getEntityTypeByName(void *self, char *name) {
    _CHECK_SELF_P(self);
    int i;
    for (i = 0; i < ((Chromosome_f *) self)->entityTypesLen; i++) {
        if (strncmp(((Chromosome_f *) self)->entityTypes[i]->type, name, (strlen(((Chromosome_f *) self)->entityTypes[i]->type) > strlen(name) ? strlen(((Chromosome_f *) self)->entityTypes[i]->type) : strlen(name))) == 0) {
            return ((Chromosome_f *) self)->entityTypes[i];
        }
    }
    return NULL;
}

/**
 * Print all data for the chromosome factory
 * 
 * @param self a chromosome factory
 * @param outFile OUT file
 * @param entFile ENT file
 */
void ChromosomePrint(void *self, FILE *outFile, FILE *entFile) {
    _CHECK_SELF_P(self);
    void **chrs;
    ChrEntity_l ent;
    int i, j, k, chrLen = 0;

    fprintf(outFile, "Gene_Id\tTranscript_Id\tChr\tLength\tCount_Reads\tTPM\tExon_Length\tExon_Count_Reads\tExon_TPM\tIntron_Length\tIntron_Count_Reads\tIntron_TPM\n");
    fprintf(entFile, "Gene_Id\tTranscript_Id\tChr\tType\tType_Number\tLength\tCount_Reads\tTPM\n");

    if (!((Chromosome_f *) self)->chrArray) {
        BtreeRecordsToArray(&(((Chromosome_f *) self)->chrArray), &(((Chromosome_f *) self)->chrArrayLen), ((Chromosome_f *) self)->chromosomes);
    }
    chrs = ((Chromosome_f *) self)->chrArray;
    for (i = 0; i < ((Chromosome_f *) self)->chrArrayLen; i++) {
        for (j = 0; j < ((Chromosome_l) chrs[i])->genesLen; j++) {
            if (((Chromosome_l) chrs[i])->genesSorted[j].count != 0) {
                fprintf(outFile, "%s\t%s\t%s\t%d\t%d\t%.4f\t%.0f\t%.0f\t%.4f\t%.0f\t%.0f\t%.4f\n",
                        ((Chromosome_l) chrs[i])->genesSorted[j].gene->geneId,
                        ((Chromosome_l) chrs[i])->genesSorted[j].gene->transcriptId,
                        ((Chromosome_l) chrs[i])->name,
                        ((Chromosome_l) chrs[i])->genesSorted[j].gene->pTo - ((Chromosome_l) chrs[i])->genesSorted[j].gene->pFrom,
                        ((Chromosome_l) chrs[i])->genesSorted[j].count,
                        ((Chromosome_l) chrs[i])->genesSorted[j].gene->TPM,
                        ((Chromosome_l) chrs[i])->genesSorted[j].gene->exonLength,
                        ((Chromosome_l) chrs[i])->genesSorted[j].gene->exonCount,
                        ((Chromosome_l) chrs[i])->genesSorted[j].gene->TPMExon,
                        ((Chromosome_l) chrs[i])->genesSorted[j].gene->intronLength,
                        ((Chromosome_l) chrs[i])->genesSorted[j].gene->intronCount,
                        ((Chromosome_l) chrs[i])->genesSorted[j].gene->TPMIntron);

                for (k = 0; k < ((Chromosome_l) chrs[i])->genesSorted[j].gene->chrEntityLen; k++) {
                    ent = ((Chromosome_l) chrs[i])->genesSorted[j].gene->chrEntity[k];
                    fprintf(entFile, "%s\t%s\t%s\t%s\t%d\t%.0f\t%.0f\t%.4f\n",
                            ((Chromosome_l) chrs[i])->genesSorted[j].gene->geneId,
                            ((Chromosome_l) chrs[i])->genesSorted[j].gene->transcriptId,
                            ((Chromosome_l) chrs[i])->name,
                            ent->type->type,
                            k + 1,
                            ent->length,
                            ent->count,
                            ent->TMP);
                }
            }
        }
    }
}

/**
 * Calculate the TPMs for all genes and entities in the chromosome factory
 * 
 * @param chrFactory a chromosome factory
 */
void TPM(Chromosome_f *chrFactory) {
    int i, j, k;
    Chromosome_l chr;
    ChrEntity_l ent;
    double T;

    T = 0.0;

    for (i = 0; i < chrFactory->chrArrayLen; i++) {
        chr = chrFactory->chrArray[i];
        for (j = 0; j < chr->genesLen; j++) {
            if (chr->genesSorted[j].gene->length <= 0) {
                chr->genesSorted[j].gene->length = chr->genesSorted[j].gene->pTo - chr->genesSorted[j].gene->pFrom;
            }
            if (chr->genesSorted[j].gene->length > 0) {
                chr->genesSorted[j].gene->intronCount = chr->genesSorted[j].gene->exonCount = 0;
                chr->genesSorted[j].gene->intronLength = chr->genesSorted[j].gene->exonLength = 0;
                for (k = 0; k < chr->genesSorted[j].gene->chrEntityLen; k++) {
                    ent = chr->genesSorted[j].gene->chrEntity[k];
                    T += (ent->count / ent->length);
                    if (strcmp(ent->type->type, "exon") == 0) {
                        chr->genesSorted[j].gene->exonCount += ent->count;
                        chr->genesSorted[j].gene->exonLength += ent->length;
                    } else if (strcmp(ent->type->type, "intron") == 0) {
                        chr->genesSorted[j].gene->intronCount += ent->count;
                        chr->genesSorted[j].gene->intronLength += ent->length;
                    }
                }
            }
        }
    }

    for (i = 0; i < chrFactory->chrArrayLen; i++) {
        chr = chrFactory->chrArray[i];
        for (j = 0; j < chr->genesLen; j++) {
            chr->genesSorted[j].gene->TPM = ((chr->genesSorted[j].count * 1.0E6) / (chr->genesSorted[j].gene->length * T));
            chr->genesSorted[j].gene->TPMExon = ((chr->genesSorted[j].gene->exonCount * 1.0E6) / (chr->genesSorted[j].gene->exonLength * T));
            chr->genesSorted[j].gene->TPMIntron = ((chr->genesSorted[j].gene->intronCount * 1.0E6) / (chr->genesSorted[j].gene->intronLength * T));
            for (k = 0; k < chr->genesSorted[j].gene->chrEntityLen; k++) {
                ent = chr->genesSorted[j].gene->chrEntity[k];
                ent->TMP = ((ent->count * 1.0E6) / (ent->length * T));
            }
        }
    }
}

/**
 * Perform all calculus developed for the Chromosome factory
 *  
 * @param self a Chromosome factory
 */
void ChrCalculus(void *self) {
    _CHECK_SELF_P(self);
    if (((Chromosome_f *) self)->chrArray == NULL) {
        BtreeRecordsToArray(&(((Chromosome_f *) self)->chrArray), &(((Chromosome_f *) self)->chrArrayLen), ((Chromosome_f *) self)->chromosomes);
    }

    TPM(((Chromosome_f *) self));
}

/**
 * Assign the read delimited by the coordinates (rFrom, rTo) to the genes and 
 * entities
 * 
 * @param chr a chromosome structure
 * @param rFrom starting position
 * @param rTo ending position
 * @param extragenic counter for extragenic regions
 */
void ReadAssignment(Chromosome_l chr, int rFrom, int rTo, int *extragenic) {
    int j, k, l, m;
    bool flag = false;

    k = chr->FindSortedIndex(chr, rTo) - 1;
    if (k < 0) k = 0;

    for (j = k; j < chr->genesLen; j++) {
        l = chr->genesSorted[j].gene->isInside(chr->genesSorted[j].gene, rFrom, rTo);
        switch (l) {
            case 1:
                flag = true;
                chr->genesSorted[j].count++;
                for (m = 0; m < chr->genesSorted[j].gene->chrEntityLen; m++) {
                    ChrEntity_l ent = chr->genesSorted[j].gene->chrEntity[m];
                    if (ent->isInside(ent, rFrom, rTo)) {
                        ent->count++;
                    }
                }
                break;
            case -1:
                j = chr->genesLen;
                break;
        }
    }
    if (!flag) {
        (*extragenic)++;
    }
}

/**
 * Creates a new Chromosome factory
 * 
 * @param gftFile GTF file descriptor
 * @return a Chromosome factory
 */
Chromosome_f* NewChromosomeFactory(FILE *gftFile) {
    Chromosome_f *factory = allocate(sizeof (Chromosome_f), __FILE__, __LINE__);

    factory->chromosomes = NULL;
    factory->chrArray = NULL;
    factory->chrArrayLen = 0;

    factory->entityTypes = allocate(sizeof (EntityType_l) * 2, __FILE__, __LINE__);
    factory->entityTypes[0] = NewEntityType("exon");
    factory->entityTypes[1] = NewEntityType("intron");
    factory->entityTypesLen = 2;

    factory->getChromosomeByName = &getChromosomeByName;
    factory->getEntityTypeByName = &getEntityTypeByName;
    factory->print = &ChromosomePrint;
    factory->calculus = &ChrCalculus;
    factory->readAssignment = &ReadAssignment;

    ChromosomeParse(factory, gftFile);

    return factory;
}

/**
 * Free a Chromosome factory
 * 
 * @param factory a Chromosome factory
 */
void FreeChromosomeFactory(Chromosome_f **factory) {
    int i;
    _CHECK_SELF_P(factory);
    if ((*factory)->chrArray) free((*factory)->chrArray);
    if ((*factory)->chromosomes) {
        BTreeFree((*factory)->chromosomes, FreeChromosome, NULL);
    }
    if ((*factory)->entityTypes) {
        for (i = 0; i < (*factory)->entityTypesLen; i++) {
            FreeEntityType((*factory)->entityTypes[i]);
        }
        free((*factory)->entityTypes);
    }
    free(*factory);
}

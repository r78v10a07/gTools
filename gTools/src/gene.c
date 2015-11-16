#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "btree.h"
#include "chrentity.h"
#include "gene.h"


/**
 * Check if the coordinates (pFrom, pTo) are inside the gene
 * 
 * @param self pointer to a gene structure
 * @param pFrom starting position
 * @param pTo ending position
 * 
 * @return 1 if the coordinates (pFrom, pTo) are inside the gene, 0 if not and 
 * -1 if the pTo position is less than the gene pFrom position  
 */
int isGeneInside(void *self, int pFrom, int pTo) {
    Gene_l gene = (Gene_l) self;
    if ((pFrom >= gene->pFrom && pFrom < gene->pTo) ||
            (pTo > gene->pFrom && pTo <= gene->pTo) ||
            (pFrom <= gene->pFrom && pTo >= gene->pTo))return 1;
    else {
        if (pTo < gene->pFrom) return -1;
    }
    return 0;
}

/**
 * Creates a new gene structure
 * 
 * @param geneString the field number 9 in the GTF format
 * 
 * @return a new gene structure
 */
Gene_l NewGene(char *geneString) {
    int i;
    Gene_l gene = NULL;
    char **fields = NULL;
    size_t fieldsSize = 0;

    fieldsSize = splitString(&fields, geneString, " ");
    if (fieldsSize > 0) {
        gene = allocate(sizeof (Gene_t), __FILE__, __LINE__);
        gene->geneId = NULL;
        gene->geneName = NULL;
        gene->transcriptId = NULL;
        gene->tssId = NULL;
        gene->chrEntity = NULL;
        gene->chrEntityLen = 0;
        gene->pFrom = -1;
        gene->pTo = -1;
        gene->length = 0;
        gene->TPM = 0.0;
        gene->TPMExon = 0.0;
        gene->TPMIntron = 0.0;
        gene->intronCount = 0;
        gene->intronLength = 0;
        gene->exonCount = 0;
        gene->exonLength = 0;
        gene->isInside = &isGeneInside;
        for (i = 0; i < fieldsSize; i += 2) {
            if (strcmp(fields[i], "gene_id") == 0) {
                gene->geneId = strdupdel(fields[i + 1], '"');
            }
            if (strcmp(fields[i], "gene_name") == 0) {
                gene->geneName = strdupdel(fields[i + 1], '"');
            }
            if (strcmp(fields[i], "transcript_id") == 0) {
                gene->transcriptId = strdupdel(fields[i + 1], '"');
            }
            if (strcmp(fields[i], "tss_id") == 0) {
                gene->tssId = strdupdel(fields[i + 1], '"');
            }
        }
        if (!gene->geneId && !gene->geneName && !gene->transcriptId && !gene->tssId) {
            free(gene);
            gene = NULL;
        }

    }
    freeArrayofPointers((void **) fields, fieldsSize);

    return gene;
}

/**
 * Free a gene structure
 * 
 * @param self a gene structure
 */
void FreeGenes(void * self) {
    int i;
    _CHECK_SELF_P(self);
    if (((Gene_l) self)->geneId) free(((Gene_l) self)->geneId);
    if (((Gene_l) self)->geneName) free(((Gene_l) self)->geneName);
    if (((Gene_l) self)->transcriptId) free(((Gene_l) self)->transcriptId);
    if (((Gene_l) self)->tssId) free(((Gene_l) self)->tssId);
    if (((Gene_l) self)->chrEntity) {
        for (i = 0; i < ((Gene_l) self)->chrEntityLen; i++) {
            FreeChrEntity(((Gene_l) self)->chrEntity[i]);
        }
        free(((Gene_l) self)->chrEntity);
    }
    free(((ChrEntity_l) self));
}

/**
 * Creates a GeneSorted_t array from a BTree with the genes
 * 
 * @param genes a GeneSorted_t array to be created
 * @param size number of elements in the GeneSorted_t array * 
 * @param root BTree root for the genes
 */
void BtreeGenesToArray(GeneSorted_t **genes, int *size, BtreeNode_t * root) {
    int i;

    if (root == NULL) {
        return;
    }
    if (!root->is_leaf) {
        for (i = 0; i <= root->num_keys; i++) {
            BtreeGenesToArray(genes, size, root->pointers[i]);
        }
    } else {
        *genes = (GeneSorted_t *) realloc(*genes, sizeof (GeneSorted_t) * (*size + root->num_keys));
        for (i = 0; i < root->num_keys; i++) {
            (*genes)[*size + i].gene = ((BtreeRecord_t *) root->pointers[i])->value;
            (*genes)[*size + i].count = 0;
        }
        *size += root->num_keys;
    }
}

/**
 * Function to compare two gene structures
 * 
 * @param key1 gene structure 1
 * @param key2 gene structure 2
 * 
 * @return 0 if geneId and geneName and transcriptId and tssId are equals in 
 * this order. Otherwise, it return an integer greater than, or less 
 * than 0, according as the string is greater than, equal to, or less than 
 * the string.
 */
int GeneKeyCMP(void *key1, void *key2) {
    int res;
    Gene_l geneKey1 = (Gene_l) key1;
    Gene_l geneKey2 = (Gene_l) key2;

    if ((res = strcmp(geneKey1->geneId, geneKey2->geneId)) == 0) {
        if ((res = strcmp(geneKey1->geneName, geneKey2->geneName)) == 0) {
            if ((res = strcmp(geneKey1->transcriptId, geneKey2->transcriptId)) == 0) {
                res = strcmp(geneKey1->tssId, geneKey2->tssId);
            }
        }
    }

    return res;
}



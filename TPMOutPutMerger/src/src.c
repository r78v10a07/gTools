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
#include <time.h>
#include <math.h>
#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "btree.h"
#include "main.h"

void FreeMEntity(void * self) {
    _CHECK_SELF_P(self);
    if (((MEntity_t *) self)->type) free(((MEntity_t *) self)->type);
    if (((MEntity_t *) self)->samples) free(((MEntity_t *) self)->samples);
}

int GeneKeyCMP(void *key1, void *key2) {
    int res;
    MGene_l geneKey1 = (MGene_l) key1;
    MGene_l geneKey2 = (MGene_l) key2;

    if ((res = strcmp(geneKey1->geneId, geneKey2->geneId)) == 0) {
        return strcmp(geneKey1->transcriptId, geneKey2->transcriptId);
    }
    return res;
}

void FreeMGene(void * self) {
    int i;
    _CHECK_SELF_P(self);
    if (((MGene_l) self)->geneId) free(((MGene_l) self)->geneId);
    if (((MGene_l) self)->chr) free(((MGene_l) self)->chr);
    if (((MGene_l) self)->transcriptId) free(((MGene_l) self)->transcriptId);
    if (((MGene_l) self)->entities) {
        for (i = 0; i < ((MGene_l) self)->entitiesLen; i++) {
            FreeMEntity(&(((MGene_l) self)->entities[i]));
        }
        free(((MGene_l) self)->entities);
    }
    if (((MGene_l) self)->CLLIntronExonRatio) free(((MGene_l) self)->CLLIntronExonRatio);
    if (((MGene_l) self)->NBCIntronExonRatio) free(((MGene_l) self)->NBCIntronExonRatio);
    if (((MGene_l) self)->repExonCount) free(((MGene_l) self)->repExonCount);
    if (((MGene_l) self)->repExonTPM) free(((MGene_l) self)->repExonTPM);
    if (((MGene_l) self)->repIntronCount) free(((MGene_l) self)->repIntronCount);
    if (((MGene_l) self)->repIntronTPM) free(((MGene_l) self)->repIntronTPM);
    free(((MGene_l) self));
}

MGene_l NewMGene(char **fields) {
    MGene_l gene = allocate(sizeof (MGene_t), __FILE__, __LINE__);

    gene->geneId = strdup(fields[0]);
    gene->transcriptId = strdup(fields[1]);
    gene->chr = strdup(fields[2]);
    gene->entities = NULL;
    gene->entitiesLen = 0;
    gene->length = -1;
    gene->CLLTPM = 0.0;
    gene->CLLTPMExon = 0.0;
    gene->CLLTPMIntron = 0.0;
    gene->CLLTMPCount = 0;
    gene->CLLexonCount = 0;
    gene->CLLintronCount = 0;
    gene->NBCTPM = 0.0;
    gene->NBCTPMExon = 0.0;
    gene->NBCTPMIntron = 0.0;
    gene->NBCTMPCount = 0;
    gene->NBCexonCount = 0;
    gene->NBCintronCount = 0;
    gene->CLLIntronExonRatio = NULL;
    gene->NBCIntronExonRatio = NULL;
    gene->repExonCount = NULL;
    gene->repExonTPM = NULL;
    gene->repIntronCount = NULL;
    gene->repIntronTPM = NULL;

    return gene;
}

double log2(const double x){
    return  log(x) * M_LOG2E;
}

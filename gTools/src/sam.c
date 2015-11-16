#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "sam.h"

Cigar_l CigarParse(int *cigarLen, char *c) {
    Cigar_l cigar = NULL;
    *cigarLen = 0;
    char *str = c;
    char *tmp = NULL;
    int tmpIndex = 0;
    int tmpLen = 0;
    
    while (*str != '\0' && *str != '\n' && *str != '\r') {
        if (isdigit(*str)) {
            if (tmpIndex >= tmpLen) {
                tmp = realloc(tmp, sizeof (char) * (tmpLen + 1));
                tmpLen++;
            }
            tmp[tmpIndex] = *str;
            tmpIndex++;
        } else {
            if (tmpIndex >= tmpLen) {
                tmp = realloc(tmp, sizeof (char) * (tmpLen + 1));
                tmpLen++;
            }
            tmp[tmpIndex] = '\0';
            tmpIndex = 0;
            cigar = realloc(cigar, sizeof (CigarValue_t) * (*cigarLen + 1));
            cigar[*cigarLen].code = *str;
            cigar[*cigarLen].value = atoi(tmp);
            (*cigarLen)++;
        }
        str++;
    }

    if (tmp) free(tmp);
    return cigar;
}

void FreeRName(void *self) {
    _CHECK_SELF_P(self);

    if (((RName_l) self)->rName) free(((RName_l) self)->rName);

    free(self);
}

RName_l NewRName(char *rName) {
    RName_l prName = allocate(sizeof (RName_t), __FILE__, __LINE__);

    prName->rName = strdup(rName);
    return prName;
}

void FreeSAM(void * self) {
    if (self) {
        if (((SAM_l) self)->qName) free(((SAM_l) self)->qName);
        if (((SAM_l) self)->cigarStr) free(((SAM_l) self)->cigarStr);
        if (((SAM_l) self)->cigar) free(((SAM_l) self)->cigar);
        if (((SAM_l) self)->rNext) free(((SAM_l) self)->rNext);
        if (((SAM_l) self)->seq) free(((SAM_l) self)->seq);
        if (((SAM_l) self)->qual) free(((SAM_l) self)->qual);

        free(self);
    }
}

SAM_l NewSAM(char **fields, RName_l rName) {
    SAM_l sam = allocate(sizeof (SAM_t), __FILE__, __LINE__);

    sam->qName = strdup(fields[0]);
    sam->flag = atoi(fields[1]);
    sam->rName = rName;
    sam->pos = atoi(fields[3]);
    sam->mapq = atoi(fields[4]);
    sam->cigarStr = strdup(fields[5]);
    sam->cigar = CigarParse(&(sam->cigarLen),fields[5]);
    sam->rNext = strdup(fields[6]);
    sam->pNext = atoi(fields[7]);
    sam->tLen = atoi(fields[8]);
    sam->seq = strdup(fields[9]);
    sam->qual = strdup(fields[10]);

    return sam;
}

/**
 * Parse SAM file
 * 
 * @param self pointer to a SAM factory
 * @param samFile SAM file to be parsed
 * @param max Number of maximun alignments to load on memory
 */
void SAMParse(void *self, FILE *samFile, int max) {
    char *line = NULL;
    size_t len = 0;
    ssize_t read = 0;
    char **fields = NULL;
    size_t fieldsSize = 0;
    RName_l rName = NULL;

    ((SAM_f *) self)->CleanSAMs(self);

    while ((read = getline(&line, &len, samFile)) != -1) {
        fieldsSize = splitString(&fields, line, "\t");
        if (fieldsSize < 11) {
            fprintf(stderr, "\n\n%s\n\n", line);
            printLog(stderr, "Bad SAM format", __FILE__, __LINE__, -1);
        }
        rName = ((SAM_f *) self)->getRNameByName(self, fields[2]);
        if (!rName) {
            ((SAM_f *) self)->rNames = reallocate(((SAM_f *) self)->rNames, sizeof (RName_t) * (((SAM_f *) self)->rNamesLen + 1), __FILE__, __LINE__);
            ((SAM_f *) self)->rNames[((SAM_f *) self)->rNamesLen] = NewRName(fields[2]);
            rName = ((SAM_f *) self)->rNames[((SAM_f *) self)->rNamesLen];
            ((SAM_f *) self)->rNamesLen++;
        }

        ((SAM_f *) self)->sams = reallocate(((SAM_f *) self)->sams, sizeof (SAM_t *) * (((SAM_f *) self)->samsLen + 1), __FILE__, __LINE__);
        ((SAM_f *) self)->sams[((SAM_f *) self)->samsLen] = NewSAM(fields, rName);
        ((SAM_f *) self)->samsLen++;

        freeArrayofPointers((void **) fields, fieldsSize);
        if (((SAM_f *) self)->samsLen == max) break;
    }

    if (line) free(line);
}

/**
 * Get SAM RName by string name
 * 
 * @param self pointer to a SAM factory
 * @param name RName to search
 * 
 * @return RName structure or NULL if it does not exist 
 */
RName_l getRNameByName(void *self, char *name) {
    _CHECK_SELF_P(self);
    int i;
    for (i = 0; i < ((SAM_f *) self)->rNamesLen; i++) {
        if (strncmp(((SAM_f *) self)->rNames[i]->rName, name, (strlen(((SAM_f *) self)->rNames[i]->rName) > strlen(name) ? strlen(((SAM_f *) self)->rNames[i]->rName) : strlen(name))) == 0) {
            return ((SAM_f *) self)->rNames[i];
        }
    }
    return NULL;
}

/**
 * Clean SAM array into the SAM factory
 * 
 * @param self pointer to a SAM factory
 */
void CleanSAMs(void *self) {
    int i;
    _CHECK_SELF_P(self);
    if (((SAM_f *) self)->sams) {
        for (i = 0; i < ((SAM_f *) self)->samsLen; i++) {
            FreeSAM(((SAM_f *) self)->sams[i]);
        }
        free(((SAM_f *) self)->sams);
    }
    ((SAM_f *) self)->sams = NULL;
    ((SAM_f *) self)->samsLen = 0;
}

/**
 * Creates a new SAM factory
 * 
 * @return a new SAM factory
 */
SAM_f *NewSAMFactory(void) {
    SAM_f *factory = allocate(sizeof (SAM_f), __FILE__, __LINE__);

    factory->parse = &SAMParse;
    factory->CleanSAMs = &CleanSAMs;
    factory->getRNameByName = &getRNameByName;
    factory->rNames = NULL;
    factory->rNamesLen = 0;
    factory->samsLen = 0;
    factory->sams = NULL;
    return factory;
}

/**
 * Free a SAM factory
 * 
 * @param self a SAM factory
 */
void FreeSAMFactory(SAM_f **self) {
    int i;
    _CHECK_SELF_P(*self);
    (*self)->CleanSAMs(*self);
    if ((*self)->rNames) {
        for (i = 0; i < (*self)->rNamesLen; i++) {
            FreeRName((*self)->rNames[i]);
        }
        free((*self)->rNames);
    }
    free(*self);
}

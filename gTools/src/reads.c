/* 
 * File:   reads.c
 * Author: veraalva
 *
 * Created on December 10, 2015, 2:06 PM
 */

#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include "berror.h"
#include "bmemory.h"
#include "btree.h"
#include "chrentity.h"
#include "gene.h"
#include "chromosome.h"
#include "sam.h"
#include "reads.h"

void processReadFromSAM(void *self, FILE *samFile, int verbose) {
    int i, p, rFrom, rTo, rLen;
    time_t now, later;
    double seconds;
    Chromosome_l chr = NULL;
    bool flag = false;
    
    time(&now);
    while (1) {
        ((Reads_f *) self)->samFactory->parse(((Reads_f *) self)->samFactory, samFile, 3000);
        if (((Reads_f *) self)->samFactory->sams == NULL) break;
        ((Reads_f *) self)->total += ((Reads_f *) self)->samFactory->samsLen;
        for (i = 0; i < ((Reads_f *) self)->samFactory->samsLen; i++) {
            if (!(chr && strncmp(chr->name, ((Reads_f *) self)->samFactory->sams[i]->rName->rName, (strlen(chr->name) >= strlen(((Reads_f *) self)->samFactory->sams[i]->rName->rName)) ? strlen(chr->name) : strlen(((Reads_f *) self)->samFactory->sams[i]->rName->rName)) == 0)) {
                chr = ((Reads_f *) self)->chrFactory->getChromosomeByName(((Reads_f *) self)->chrFactory, ((Reads_f *) self)->samFactory->sams[i]->rName->rName);
            }
            if (!chr) {
                ((Reads_f *) self)->errors++;
                fprintf(((Reads_f *) self)->errFile, "RNAME\t%s\tquery\t%s\tChromosome not in GTF file\n", ((Reads_f *) self)->samFactory->sams[i]->rName->rName, ((Reads_f *) self)->samFactory->sams[i]->qName);
            } else {
                if (((Reads_f *) self)->samFactory->sams[i]->cigarLen == 1) {
                    if (((Reads_f *) self)->samFactory->sams[i]->cigar[0].code == 'M') {
                        ((Reads_f *) self)->chrFactory->readAssignment(chr, ((Reads_f *) self)->samFactory->sams[i]->pos, ((Reads_f *) self)->samFactory->sams[i]->pos + ((Reads_f *) self)->samFactory->sams[i]->cigar[0].value - 1, &(((Reads_f *) self)->intergenic));
                    } else {
                        ((Reads_f *) self)->errors++;
                        fprintf(((Reads_f *) self)->errFile, "CIGAR\t%c\t%s\tquery\t%s\n", ((Reads_f *) self)->samFactory->sams[i]->cigar[0].code, ((Reads_f *) self)->samFactory->sams[i]->cigarStr, ((Reads_f *) self)->samFactory->sams[i]->qName);
                    }
                } else if (((Reads_f *) self)->samFactory->sams[i]->cigarLen > 1) {
                    flag = false;
                    for (p = 0; p < ((Reads_f *) self)->samFactory->sams[i]->cigarLen; p++) {
                        if (((Reads_f *) self)->samFactory->sams[i]->cigar[p].code == 'N') {
                            flag = true;
                            break;
                        }
                    }
                    if (flag) {
                        rLen = 0;
                        rFrom = ((Reads_f *) self)->samFactory->sams[i]->pos;
                        for (p = 0; p < ((Reads_f *) self)->samFactory->sams[i]->cigarLen; p++) {
                            rLen += ((Reads_f *) self)->samFactory->sams[i]->cigar[p].value;
                            if (((Reads_f *) self)->samFactory->sams[i]->cigar[p].code == 'N') {
                                ((Reads_f *) self)->chrFactory->readAssignment(chr, rFrom, rTo, &(((Reads_f *) self)->intergenic));
                                rFrom = rTo + 1 + ((Reads_f *) self)->samFactory->sams[i]->cigar[p].value;
                                rLen = 0;
                            } else {
                                rTo = rFrom + rLen - 1;
                            }
                        }
                        ((Reads_f *) self)->chrFactory->readAssignment(chr, rFrom, rTo, &(((Reads_f *) self)->intergenic));
                    } else {
                        ((Reads_f *) self)->chrFactory->readAssignment(chr, ((Reads_f *) self)->samFactory->sams[i]->pos, ((Reads_f *) self)->samFactory->sams[i]->pos + (int) strlen(((Reads_f *) self)->samFactory->sams[i]->seq) - 1, &(((Reads_f *) self)->intergenic));
                    }
                } else {
                    ((Reads_f *) self)->errors++;
                    fprintf(((Reads_f *) self)->errFile, "NOCIGAR\t%s\tquery\t%s\n", ((Reads_f *) self)->samFactory->sams[i]->cigarStr, ((Reads_f *) self)->samFactory->sams[i]->qName);
                }
            }
        }
        time(&later);
        seconds = difftime(later, now);
        if (verbose)printf("\tTotal reads: %10d, intergenic reads: %10d, rName errors: %5d,  elapse time: %4.0f s, %5.0f reads/s\r", ((Reads_f *) self)->total, ((Reads_f *) self)->intergenic, ((Reads_f *) self)->errors, seconds, ((Reads_f *) self)->total / seconds);
    }
    
    if (verbose)printf("\tTotal reads: %10d, intergenic reads: %10d, rName errors: %5d,  elapse time: %4.0f s, %5.0f reads/s\n", ((Reads_f *) self)->total, ((Reads_f *) self)->intergenic, ((Reads_f *) self)->errors, seconds, ((Reads_f *) self)->total / seconds);
}

Reads_f *NewReadsFactory(SAM_f *samFactory, Chromosome_f *chrFactory, FILE *errFile) {
    Reads_f *readsFactory = allocate(sizeof (Reads_f), __FILE__, __LINE__);

    readsFactory->samFactory = samFactory;
    readsFactory->chrFactory = chrFactory;
    readsFactory->errFile = errFile;
    
    readsFactory->total = 0;
    readsFactory->errors = 0;
    readsFactory->intergenic = 0;

    readsFactory->processReadFromSAM = &processReadFromSAM;

    return readsFactory;
}

void FreeReadsFactory(void *self) {
    _CHECK_SELF_P(self);
    free(((Reads_f *) self));
}
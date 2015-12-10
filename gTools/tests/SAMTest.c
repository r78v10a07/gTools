/*
 * File:   SAMTest.c
 * Author: veraalva
 *
 * Created on Dec 10, 2015, 10:37:51 AM
 */
#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <CUnit/Basic.h>
#include <stdbool.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <unistd.h>
#include "berror.h"
#include "bmemory.h"
#include "btree.h"
#include "chrentity.h"
#include "gene.h"
#include "chromosome.h"
#include "sam.h"
#include "reads.h"

/*
 * CUnit Test Suite
 */

int init_suite(void) {
    return 0;
}

int clean_suite(void) {
    return 0;
}

void test1() {
    int i, j, k, p, max = 3000;
    bool flag = false;
    Chromosome_l chr = NULL;
    void **chrs;
    ChrEntity_l ent;
    int chrLen = 0;
    char *message = allocate(sizeof (char) * 10000, __FILE__, __LINE__);
    int rFrom, rTo, rLen, errors = 0;
    FILE *gtfFile = fopen("resources/test.gtf", "r");
    FILE *samFile = fopen("resources/test.sam", "r");
    CU_TEST_FATAL(gtfFile != NULL);
    CU_TEST_FATAL(samFile != NULL);

    printf("\nReading chromosomes from GTF file\n");
    Chromosome_f *chrFactory = NewChromosomeFactory(gtfFile);
    CU_TEST_FATAL(chrFactory != NULL);

    BtreeRecordsToArray(&(chrFactory->chrArray), &(chrFactory->chrArrayLen), chrFactory->chromosomes);
    CU_TEST_FATAL(chrFactory->chrArray != NULL);
    CU_TEST_FATAL(chrFactory->chrArrayLen == 1);
    printf("Chromosomes %d\n", chrFactory->chrArrayLen);

    printf("Parsing sam file\n");
    SAM_f *samFactory = NewSAMFactory();
    CU_TEST_FATAL(samFactory != NULL);
    
    Reads_f *readsFactory = NewReadsFactory(samFactory, chrFactory, stderr);
    CU_TEST_FATAL(readsFactory != NULL);
    readsFactory->processReadFromSAM(readsFactory, samFile, 1);
    printf("Reads %d parsed\n", readsFactory->total);
    CU_TEST_FATAL(readsFactory->total == 12);
    CU_TEST_FATAL(readsFactory->errors == 1);

    chrs = chrFactory->chrArray;
    for (i = 0; i < chrFactory->chrArrayLen; i++) {
        for (j = 0; j < ((Chromosome_l) chrs[i])->genesLen; j++) {
            if (((Chromosome_l) chrs[i])->genesSorted[j].count != 0) {
                for (k = 0; k < ((Chromosome_l) chrs[i])->genesSorted[j].gene->chrEntityLen; k++) {
                    ent = ((Chromosome_l) chrs[i])->genesSorted[j].gene->chrEntity[k];
                    printf("%s\t%s\t%s\t%s\t%d\t%d\t%d\n",
                            ((Chromosome_l) chrs[i])->genesSorted[j].gene->geneId,
                            ((Chromosome_l) chrs[i])->genesSorted[j].gene->transcriptId,
                            ((Chromosome_l) chrs[i])->name,
                            ent->type->type,
                            k + 1,
                            ent->length,
                            ent->count);
                    if (strcmp(((Chromosome_l) chrs[i])->genesSorted[j].gene->geneId, "GEN1") == 0) {                        
                        if (strcmp(((Chromosome_l) chrs[i])->genesSorted[j].gene->transcriptId, "TRANS1") == 0) {
                            CU_TEST_FATAL(((Chromosome_l) chrs[i])->genesSorted[j].gene->length == 2845);
                            CU_TEST_FATAL(((Chromosome_l) chrs[i])->genesSorted[j].gene->chrEntityLen == 7);
                            switch (k) {
                                case 0: // EXON
                                    CU_TEST_FATAL(ent->count == 4);
                                    CU_TEST_FATAL(ent->length == 353);
                                    break;
                                case 1: // INTRON
                                    CU_TEST_FATAL(ent->count == 2);
                                    CU_TEST_FATAL(ent->length == 385);
                                    break;
                                case 2: // EXON
                                    CU_TEST_FATAL(ent->count == 2);
                                    CU_TEST_FATAL(ent->length == 108);
                                    break;
                                case 3: // INTRON
                                    CU_TEST_FATAL(ent->count == 0);
                                    CU_TEST_FATAL(ent->length == 499);
                                    break;
                                case 4: // EXON
                                    CU_TEST_FATAL(ent->count == 3);
                                    CU_TEST_FATAL(ent->length == 1188);
                                    break;
                                case 5: // INTRON
                                    CU_TEST_FATAL(ent->count == 1);
                                    CU_TEST_FATAL(ent->length == 199);
                                    break;
                                case 6: // EXON
                                    CU_TEST_FATAL(ent->count == 2);
                                    CU_TEST_FATAL(ent->length == 110);
                                    break;
                            }
                        }
                        if (strcmp(((Chromosome_l) chrs[i])->genesSorted[j].gene->transcriptId, "TRANS2") == 0) {
                            CU_TEST_FATAL(((Chromosome_l) chrs[i])->genesSorted[j].gene->length == 2106);
                            CU_TEST_FATAL(((Chromosome_l) chrs[i])->genesSorted[j].gene->chrEntityLen == 5);
                            switch (k) {
                                case 0: // EXON
                                    CU_TEST_FATAL(ent->count == 2);
                                    CU_TEST_FATAL(ent->length == 108);
                                    break;
                                case 1: // INTRON
                                    CU_TEST_FATAL(ent->count == 0);
                                    CU_TEST_FATAL(ent->length == 499);
                                    break;
                                case 2: // EXON
                                    CU_TEST_FATAL(ent->count == 3);
                                    CU_TEST_FATAL(ent->length == 1188);
                                    break;
                                case 3: // INTRON
                                    CU_TEST_FATAL(ent->count == 1);
                                    CU_TEST_FATAL(ent->length == 199);
                                    break;
                                case 4: // EXON
                                    CU_TEST_FATAL(ent->count == 2);
                                    CU_TEST_FATAL(ent->length == 110);
                                    break;
                            }
                        }
                        if (strcmp(((Chromosome_l) chrs[i])->genesSorted[j].gene->transcriptId, "TRANS3") == 0) {
                            CU_TEST_FATAL(((Chromosome_l) chrs[i])->genesSorted[j].gene->length == 1796);
                            CU_TEST_FATAL(((Chromosome_l) chrs[i])->genesSorted[j].gene->chrEntityLen == 3);
                            switch (k) {
                                case 0: // EXON
                                    CU_TEST_FATAL(ent->count == 2);
                                    CU_TEST_FATAL(ent->length == 108);
                                    break;
                                case 1: // INTRON
                                    CU_TEST_FATAL(ent->count == 0);
                                    CU_TEST_FATAL(ent->length == 499);
                                    break;
                                case 2: // EXON
                                    CU_TEST_FATAL(ent->count == 3);
                                    CU_TEST_FATAL(ent->length == 1188);
                                    break;
                            }
                        }
                    }
                }
            }
        }
    }

    CU_TEST_FATAL(readsFactory->intergenic == 2);
    printf("Total reads: %10d, intergenic reads: %10d, errors: %5d (testing errors, it should be 1)\n\n", readsFactory->total, readsFactory->intergenic, readsFactory->errors);

    free(message);
    FreeSAMFactory(&samFactory);
    FreeChromosomeFactory(&chrFactory);
    fclose(gtfFile);
    fclose(samFile);
}

int main() {
    CU_pSuite pSuite = NULL;

    /* Initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    /* Add a suite to the registry */
    pSuite = CU_add_suite("SAMTest", init_suite, clean_suite);
    if (NULL == pSuite) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    /* Add the tests to the suite */
    if (NULL == CU_add_test(pSuite, "test1", test1)) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    /* Run all tests using the CUnit Basic interface */
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    return CU_get_error();
}

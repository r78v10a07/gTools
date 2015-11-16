/* 
 * File:   main.c
 * Author: veraalva
 *
 * Created on October 28, 2015, 12:35 PM
 */
#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
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

char *program_name;

void print_usage(FILE *stream, int exit_code) {
    fprintf(stream, "\n********************************************************************************\n");
    fprintf(stream, "\nUsage: %s \n", program_name);
    fprintf(stream, "\n\n%s options:\n\n", program_name);
    fprintf(stream, "-v,   --verbose                     Print info\n");
    fprintf(stream, "-h,   --help                        Display this usage information.\n");
    fprintf(stream, "-g,   --gtf                         GTF file\n");
    fprintf(stream, "-i,   --input                       Input file. Use - for stdin\n");
    fprintf(stream, "-o,   --output                      Output prefix\n");
    fprintf(stream, "********************************************************************************\n");
    fprintf(stream, "\n            Roberto Vera Alvarez (e-mail: r78v10a07@gmail.com)\n\n");
    fprintf(stream, "********************************************************************************\n");
    exit(0);
}

/*
 * 
 */
int main(int argc, char** argv) {
    int i, p;    
    int next_option, verbose;
    const char* const short_options = "vhg:i:o:";
    char *outputName, *output = NULL, *input = NULL;
    FILE *gtfFile = NULL;
    FILE *inFile = NULL;
    FILE *outFile = NULL;
    FILE *errFile = NULL;
    FILE *entFile = NULL;
    int max = 3000;
    Chromosome_l chr = NULL;
    int rFrom, rTo, rLen, total = 0, errors = 0;
    time_t now, later;
    double seconds;
    int extragenic = 0;
    bool flag = false;
    
    program_name = argv[0];
    
    const struct option long_options[] = {
        { "help", 0, NULL, 'h'},
        { "verbose", 0, NULL, 'v'},
        { "gtf", 1, NULL, 'g'},
        { "output", 1, NULL, 'o'},
        { "input", 1, NULL, 'i'},
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

            case 'o':
                output = strdup(optarg);
                break;
                
            case 'i':
                if (strlen(optarg) == 1 && optarg[0] == '-'){
                    inFile = stdin;
                }else{
                    inFile = checkPointerError(fopen(optarg, "r"), "Can't open INPUT file", __FILE__, __LINE__, -1);
                }
                break;

            case 'g':
                gtfFile = checkPointerError(fopen(optarg, "r"), "Can't open GTF file", __FILE__, __LINE__, -1);
                break;
        }
    } while (next_option != -1);

    if (!gtfFile || !output || ! inFile) {
        print_usage(stderr, -1);
    }
    outputName = allocate(sizeof(char) * (strlen(output) + 10), __FILE__, __LINE__);
    sprintf(outputName,"%s.out", output);
    outFile = checkPointerError(fopen(outputName, "w"), "Can't open OUT file", __FILE__, __LINE__, -1);
    sprintf(outputName,"%s.err", output);
    errFile = checkPointerError(fopen(outputName, "w"), "Can't open ERR file", __FILE__, __LINE__, -1);
    sprintf(outputName,"%s.ent", output);
    entFile = checkPointerError(fopen(outputName, "w"), "Can't open ENT file", __FILE__, __LINE__, -1);
    free(outputName);
    free(output);

    time(&now);
    if (verbose) printf("Reading chromosomes from GTF file\n");
    Chromosome_f *chrFactory = NewChromosomeFactory(gtfFile);
    time(&later);
    if (verbose) printf("Chromosomes loaded in %.0f s\n", difftime(later, now));

    time(&now);
    if (verbose) printf("Parsing sam file\n");
    SAM_f *samFactory = NewSAMFactory();
    while (1) {
        samFactory->parse(samFactory, stdin, max);
        if (samFactory->sams == NULL) break;
        total += samFactory->samsLen;
        for (i = 0; i < samFactory->samsLen; i++) {
            if (!(chr && strncmp(chr->name, samFactory->sams[i]->rName->rName, (strlen(chr->name) >= strlen(samFactory->sams[i]->rName->rName)) ? strlen(chr->name) : strlen(samFactory->sams[i]->rName->rName)) == 0)) {
                chr = chrFactory->getChromosomeByName(chrFactory, samFactory->sams[i]->rName->rName);
            }
            if (!chr) {
                errors++;
                fprintf(errFile, "RNAME\t%s\tquery\t%s\n", samFactory->sams[i]->rName->rName, samFactory->sams[i]->qName);
            } else {
                if (samFactory->sams[i]->cigarLen == 1) {
                    if (samFactory->sams[i]->cigar[0].code == 'M') {
                        chrFactory->readAssignment(chr, samFactory->sams[i]->pos, samFactory->sams[i]->pos + samFactory->sams[i]->cigar[0].value - 1, &extragenic);
                    } else {
                        errors++;
                        fprintf(errFile, "CIGAR\t%c\t%s\tquery\t%s\n", samFactory->sams[i]->cigar[0].code, samFactory->sams[i]->cigarStr, samFactory->sams[i]->qName);
                    }
                } else if (samFactory->sams[i]->cigarLen > 1) {
                    flag = false;
                    for (p = 0; p < samFactory->sams[i]->cigarLen; p++) {
                        if (samFactory->sams[i]->cigar[p].code == 'N') {
                            flag = true;
                            break;
                        }
                    }
                    if (flag) {
                        rLen = 0;
                        rFrom = samFactory->sams[i]->pos;
                        for (p = 0; p < samFactory->sams[i]->cigarLen; p++) { 
                            rLen += samFactory->sams[i]->cigar[p].value;
                            if (samFactory->sams[i]->cigar[p].code == 'N') {
                                chrFactory->readAssignment(chr, rFrom, rTo, &extragenic);
                                rFrom = rTo + 1 + samFactory->sams[i]->cigar[p].value;   
                                rLen = 0;
                            }else{
                                rTo = rFrom + rLen - 1;
                            }
                        }
                        chrFactory->readAssignment(chr, rFrom, rTo, &extragenic);
                    }else{
                        chrFactory->readAssignment(chr, samFactory->sams[i]->pos, samFactory->sams[i]->pos + (int) strlen(samFactory->sams[i]->seq) - 1, &extragenic);
                    }
                } else {
                    errors++;
                    fprintf(errFile, "NOCIGAR\t%s\tquery\t%s\n", samFactory->sams[i]->cigarStr, samFactory->sams[i]->qName);
                }
            }
        }
        time(&later);
        seconds = difftime(later, now);
        if (verbose)printf("\tTotal reads: %10d, extragenic reads: %10d, rName errors: %5d,  elapse time: %4.0f s, %5.0f reads/s\r", total, extragenic, errors, seconds, total / seconds);
    }
    time(&later);
    seconds = difftime(later, now);
    printf("\tTotal reads: %10d, extragenic reads: %10d, rName errors: %5d,  elapse time: %4.0f s, %5.0f reads/s\n\n", total, extragenic, errors, seconds, total / seconds);

    chrFactory->calculus(chrFactory);
    chrFactory->print(chrFactory, outFile, entFile);

    FreeSAMFactory(&samFactory);
    FreeChromosomeFactory(&chrFactory);
    fclose(gtfFile);
    fclose(outFile);
    fclose(errFile);
    fclose(entFile);
    return (EXIT_SUCCESS);
}
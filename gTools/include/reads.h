/* 
 * File:   reads.h
 * Author: veraalva
 *
 * Created on December 10, 2015, 2:06 PM
 */

#ifndef READS_H
#define READS_H

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct {
        SAM_f *samFactory;
        Chromosome_f *chrFactory;
        FILE *errFile;
        int total;
        int errors;
        int intergenic;

        void (*processReadFromSAM)(void *self, FILE *samFile, int verbose);
    } Reads_f;

    extern Reads_f *NewReadsFactory(SAM_f *samFactory, Chromosome_f *chrFactory, FILE *errFile);

    extern void FreeReadsFactory(void *self);


#ifdef __cplusplus
}
#endif

#endif /* READS_H */


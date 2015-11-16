/* 
 * File:   main.h
 * Author: veraalva
 *
 * Created on November 16, 2015, 2:17 PM
 */

#ifndef MAIN_H
#define	MAIN_H

#ifdef	__cplusplus
extern "C" {
#endif

    typedef struct {
        char *sample;
        int count;
        double TPM;
    } MEntitySample_t;

    typedef struct {
        char *type;
        int length;
        MEntitySample_t *samples;
    } MEntity_t;

    typedef struct {
        char *geneId;
        char *transcriptId;
        char *chr;
        MEntity_t *entities;
        int entitiesLen;
    } MGene_t;

    typedef MGene_t *MGene_l;

    extern void FreeMEntity(void * self);
    extern int GeneKeyCMP(void *key1, void *key2);
    extern void FreeMGene(void * self);
    extern MGene_l NewMGene(char **fields);


#ifdef	__cplusplus
}
#endif

#endif	/* MAIN_H */


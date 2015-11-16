/* 
 * File:   chromosome.h
 * Author: veraalva
 *
 * Created on October 28, 2015, 12:36 PM
 */

#ifndef CHROMOSOME_H
#define	CHROMOSOME_H

#ifdef	__cplusplus
extern "C" {
#endif  

    typedef struct {
        char *name;
        BtreeNode_t *genesTree;
        GeneSorted_t *genesSorted;
        BtreeNode_t *index;
        int genesLen;
        int chromosomeLen;

        int (*FindSortedIndex)(void *self, int key);
    } Chromosome_t;

    typedef Chromosome_t *Chromosome_l;

    /*
     * Chromosome Factory to encapsulate all functions
     */
    typedef struct {
        BtreeNode_t *chromosomes;
        void **chrArray;
        int chrArrayLen;
        EntityType_t **entityTypes;
        int entityTypesLen;

        /**
         * Get an Chromosome structure from a name
         * 
         * @param self a chromosome factory
         * @param name Chromosome name
         * @return an Chromosome structure or NULL
         */
        Chromosome_l(*getChromosomeByName)(void *self, char *name);

        /**
         * Get an EntityType structure from a name
         * 
         * @param self a chromosome factory
         * @param name EntityType name
         * @return an EntityType structure or NULL
         */
        EntityType_l(*getEntityTypeByName)(void *self, char *name);

        /**
         * Print all data for the chromosome factory
         * 
         * @param self a chromosome factory
         * @param outFile OUT file
         * @param entFile ENT file
         */
        void (*print)(void *self, FILE *outFile, FILE *entFile);

        /**
         * Perform all calculus developed for the Chromosome factory
         *  
         * @param self a Chromosome factory
         */
        void (*calculus)(void *self);

        /**
         * Assign the read delimited by the coordinates (rFrom, rTo) to the genes and 
         * entities
         * 
         * @param chr a chromosome structure
         * @param rFrom starting position
         * @param rTo ending position
         * @param extragenic counter for extragenic regions
         */
        void (*readAssignment)(Chromosome_l chr, int rFrom, int rTo, int *extragenic);
    } Chromosome_f;

    /**
     * Creates a new Chromosome factory
     * 
     * @param gftFile GTF file descriptor
     * @return a Chromosome factory
     */
    extern Chromosome_f* NewChromosomeFactory(FILE *gftFile);

    /**
     * Free a Chromosome factory
     * 
     * @param factory a Chromosome factory
     */
    extern void FreeChromosomeFactory(Chromosome_f **factory);


#ifdef	__cplusplus
}
#endif

#endif	/* CHROMOSOME_H */


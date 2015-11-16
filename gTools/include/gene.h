/* 
 * File:   gene.h
 * Author: veraalva
 *
 * Created on November 3, 2015, 12:26 PM
 */

#ifndef GENE_H
#define	GENE_H

#ifdef	__cplusplus
extern "C" {
#endif

    /*
     * Structure to storage the GTF info
     */
    typedef struct {
        char *geneId;
        char *geneName;
        char *transcriptId;
        char *tssId;
        ChrEntity_t **chrEntity;
        int chrEntityLen;
        int pFrom;
        int pTo;
        double length;
        double exonLength;
        double exonCount;
        double intronCount;
        double intronLength;
        double TPM;
        double TPMExon;
        double TPMIntron;

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
        int (*isInside)(void *self, int pFrom, int pTo);
    } Gene_t;

    typedef Gene_t *Gene_l;

    /*
     * Structure used to order the genes
     */
    typedef struct {
        Gene_l gene;
        int count;
    } GeneSorted_t;

    /**
     * Creates a new gene structure
     * 
     * @param geneString the field number 9 in the GTF format
     * 
     * @return a new gene structure
     */
    extern Gene_l NewGene(char *geneString);

    /**
     * Free a gene structure
     * 
     * @param self a gene structure
     */
    extern void FreeGenes(void * self);

    /**
     * Creates a GeneSorted_t array from a BTree with the genes
     * 
     * @param genes a GeneSorted_t array to be created
     * @param size number of elements in the GeneSorted_t array * 
     * @param root BTree root for the genes
     */
    extern void BtreeGenesToArray(GeneSorted_t **index, int *size, BtreeNode_t * root);

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
    extern int GeneKeyCMP(void *key1, void *key2);


#ifdef	__cplusplus
}
#endif

#endif	/* GENE_H */


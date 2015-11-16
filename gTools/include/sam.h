/* 
 * File:   sam.h
 * Author: veraalva
 *
 * Created on November 3, 2015, 1:15 PM
 */

#ifndef SAM_H
#define	SAM_H

#ifdef	__cplusplus
extern "C" {
#endif

    /*
     * Cigar structure     
     */
    typedef struct {
        char code;
        int value;
    } CigarValue_t;

    typedef CigarValue_t *Cigar_l;

    typedef struct {
        char *rName;
    } RName_t;

    typedef RName_t *RName_l;

    /*
     * SAM Structure
     */
    typedef struct {
        char *qName;
        int flag;
        RName_l rName;
        int pos;
        int mapq;
        char *cigarStr;
        Cigar_l cigar;
        int cigarLen;
        char *rNext;
        int pNext;
        int tLen;
        char *seq;
        char *qual;
    } SAM_t;

    typedef SAM_t *SAM_l;

    /**
     * SAM factory
     */
    typedef struct {
        SAM_t **sams;
        int samsLen;
        RName_t **rNames;
        int rNamesLen;

        /**
         * Parse SAM file
         * 
         * @param self pointer to a SAM factory
         * @param samFile SAM file to be parsed
         * @param max Number of maximun alignments to load on memory
         */
        void (*parse)(void *self, FILE *samFile, int max);

        /**
         * Get SAM RName by string name
         * 
         * @param self pointer to a SAM factory
         * @param name RName to search
         * 
         * @return RName structure or NULL if it does not exist 
         */
        RName_l(*getRNameByName)(void *self, char *name);

        /**
         * Clean SAM array into the SAM factory
         * 
         * @param self pointer to a SAM factory
         */
        void (*CleanSAMs)(void *self);

    } SAM_f;

    /**
     * Creates a new SAM factory
     * 
     * @return a new SAM factory
     */
    extern SAM_f *NewSAMFactory(void);

    /**
     * Free a SAM factory
     * 
     * @param self a SAM factory
     */
    extern void FreeSAMFactory(SAM_f **self);


#ifdef	__cplusplus
}
#endif

#endif	/* SAM_H */


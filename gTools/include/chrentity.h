/* 
 * File:   chrentity.h
 * Author: veraalva
 *
 * Created on November 3, 2015, 12:27 PM
 */

#ifndef CHRENTITY_H
#define	CHRENTITY_H

#ifdef	__cplusplus
extern "C" {
#endif

    typedef struct {
        char *type;
    } EntityType_t;

    typedef EntityType_t *EntityType_l;

    /**
     * Create a new Entity type
     * 
     * @param type string with the type to be created
     * @return a new Entity type
     */
    extern EntityType_l NewEntityType(char *type);

    /**
     * Free an Entity type
     * 
     * @param self an Entity type
     */
    extern void FreeEntityType(void *self);

    typedef struct ChrEntity_s {
        int pFrom;
        int pTo;
        double count;
        double length;
        EntityType_t *type;
        char **fields;
        int fieldsLen;
        double TMP;

        /**
         * Check if the coordinates (pFrom, pTo) are inside the entity
         * 
         * @param self pointer to a entity structure
         * @param pFrom starting position
         * @param pTo ending position
         * 
         * @return true if the coordinates (pFrom, pTo) are inside the entity, false if 
         * not and  
         */
        bool (*isInside)(void *self, int pFrom, int pTo);
    } ChrEntity_t;

    typedef ChrEntity_t *ChrEntity_l;

    /**
     * Create a new Entity structure
     * 
     * @param entType EntityType to be used
     * @param fields array of string with the entity data
     * 
     * @return a new Entity structure
     */
    extern ChrEntity_l NewChrEntity(EntityType_l entType, char **fields);

    /**
     * Free an Entity structure
     * 
     * @param self an Entity structure
     */
    extern void FreeChrEntity(void *self);


#ifdef	__cplusplus
}
#endif

#endif	/* CHRENTITY_H */


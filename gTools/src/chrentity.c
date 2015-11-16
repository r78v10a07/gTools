#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "btree.h"
#include "chrentity.h"

/********************************************
 * Entity Functions
 ********************************************/

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
bool isChrInside(void *self, int pFrom, int pTo) {
    ChrEntity_l chrEntity = (ChrEntity_l) self;
    if ((pFrom >= chrEntity->pFrom && pFrom < chrEntity->pTo) ||
            (pTo > chrEntity->pFrom && pTo <= chrEntity->pTo) ||
            (pFrom <= chrEntity->pFrom && pTo >= chrEntity->pTo))return true;

    return false;
}

/**
 * Create a new Entity structure
 * 
 * @param entType EntityType to be used
 * @param fields array of string with the entity data
 * 
 * @return a new Entity structure
 */
ChrEntity_l NewChrEntity(EntityType_l entType, char **fields) {
    ChrEntity_l chrEntity = allocate(sizeof (ChrEntity_t), __FILE__, __LINE__);

    chrEntity->type = entType;
    chrEntity->isInside = &isChrInside;
    if (fields) {
        chrEntity->fields = allocate(sizeof (char *) * 3, __FILE__, __LINE__);
        chrEntity->fields[0] = strdup(fields[5]);
        chrEntity->fields[1] = strdup(fields[6]);
        chrEntity->fields[2] = strdup(fields[7]);
        chrEntity->fieldsLen = 3;
        chrEntity->pFrom = atof(fields[3]);
        chrEntity->pTo = atof(fields[4]);
        chrEntity->length = chrEntity->pTo - chrEntity->pFrom;
    } else {
        chrEntity->fields = NULL;
        chrEntity->fieldsLen = 0;
        chrEntity->pFrom = 0;
        chrEntity->pTo = 0;
    }
    chrEntity->count = 0;
    chrEntity->TMP = 0.0;
    return chrEntity;
}

/**
 * Free an Entity structure
 * 
 * @param self an Entity structure
 */
void FreeChrEntity(void *self) {
    _CHECK_SELF_P(self);
    if (((ChrEntity_l) self)->fields) freeArrayofPointers((void **) ((ChrEntity_l) self)->fields, ((ChrEntity_l) self)->fieldsLen);
    free(((ChrEntity_l) self));
}

/********************************************
 * Entity Type Functions
 ********************************************/

/**
 * Create a new Entity type
 * 
 * @param type string with the type to be created
 * @return a new Entity type
 */
EntityType_l NewEntityType(char *type) {
    EntityType_l entType = allocate(sizeof (EntityType_t), __FILE__, __LINE__);

    entType->type = strdup(type);

    return entType;
}

/**
 * Free an Entity type
 * 
 * @param self an Entity type
 */
void FreeEntityType(void *self) {
    _CHECK_SELF_P(self);
    if (((EntityType_l) self)->type) free(((EntityType_l) self)->type);
    free(((EntityType_l) self));
}

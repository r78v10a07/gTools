
#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "btree.h"
#include "main.h"

double HausdorffDistance(double *s1, int s1Len, double *s2, int s2Len) {
    int i, j;
    double infimum, supremum_s1, supremum_s2;
    double **matrix = allocate(sizeof (double **) * s1Len, __FILE__, __LINE__);

    supremum_s1 = NAN;
    for (i = 0; i < s1Len; i++) {
        matrix[i] = allocate(sizeof (double *) * s2Len, __FILE__, __LINE__);
        infimum = NAN;
        for (j = 0; j < s2Len; j++) {
            matrix[i][j] = fabs(s1[i] - s2[j]);
            if ((!isnan(matrix[i][j]) && !isinf(matrix[i][j])) && (isnan(infimum) || matrix[i][j] < infimum)) infimum = matrix[i][j];
        }
        if (!isnan(infimum) && (isnan(supremum_s1) || supremum_s1 < infimum)) supremum_s1 = infimum;
    }

    supremum_s2 = NAN;
    for (j = 0; j < s2Len; j++) {
        infimum = NAN;
        for (i = 0; i < s1Len; i++) {
            if ((!isnan(matrix[i][j]) && !isinf(matrix[i][j])) && (isnan(infimum) || matrix[i][j] < infimum)) infimum = matrix[i][j];
        }
        if (!isnan(infimum) && (isnan(supremum_s2) || supremum_s2 < infimum)) supremum_s2 = infimum;
    }

    for (i = 0; i < s1Len; i++) free(matrix[i]);
    free(matrix);

    return (supremum_s1 > supremum_s2) ? supremum_s1 : supremum_s2;
}

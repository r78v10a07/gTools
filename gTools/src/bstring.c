/* 
 * File:   bstring.c
 * Author: roberto
 *
 * Created on April 15, 2014, 3:05 PM
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "berror.h"


/**
 * Split a string into a array of string using the delimiter
 * 
 * @param dest the destination pointer
 * @param src the source pointer
 * @param delimiter the delimiter
 * @return a array with the splited strings
 */
size_t splitString(char ***dest, char *src, char *delimiter) {
    size_t count = 0;
    char *token;
    char *srccpy;

    srccpy = strdup(src);
    *dest = NULL;
    token = strtok(srccpy, delimiter);
    while (token) {
        *dest = (char **) checkPointerError(realloc(*dest, sizeof (char **) * (count + 1)), "Can't allocate memory", __FILE__, __LINE__, -1);
        (*dest)[count] = strdup(token);
        token = strtok(NULL, delimiter);
        count++;
    }
    if (srccpy)free(srccpy);
    return count;
}

/**
 * Return 0 if the string haystack ends with the string needle
 * 
 * @param haystack the string to be analyzed
 * @param needle the suffix string
 * @return 0 if the string haystack ends with the string needle, 1 if not
 */
int strbcmp(const char *haystack, const char *needle) {
    int length;
    char *sub;
    if (haystack && needle && strlen(haystack) >= (length = strlen(needle)) && (sub = strstr(haystack, needle)) && strlen(sub) == length) return 0;
    return 1;
}

/**
 * Duplicate the src string until the delimiter or the newline or \0 are found
 * 
 * @param src string to duplicate
 * @param delimiter delimiter to stop duplication
 * @return string
 */
char *strdupdel(char *src, char delimiter) {
    char *id = src;
    char *val, *val_src;
    char *result = NULL;
    if (id) {
        while(*(id++) != delimiter);
        val = val_src = (char *) checkPointerError(malloc(sizeof (char) * (strlen(id) + 1)), "Can't allocate memory", __FILE__, __LINE__, -1);
        
        while (*id != delimiter && *id != '\0' && *id != '\n' && *id != '\r') {
            *(val++) = *(id++);
        }
        *val = '\0';
        result = strdup(val_src);
        free(val_src);
    }
    return result;
}
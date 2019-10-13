#ifndef HASH_TABLE_H
#define HASH_TABLE_H

#include <stdints.h>

/*
 * hash string to void *
 * supply default string hash function
 * hash to linked list of entries
 * constructor
 * destructor
 * handle resize and rehash
 * handle insert
 * handle peek
 * handle remove
 */

/*
 * https://stackoverflow.com/questions/7666509/hash-function-for-string
 */
int
hash(char const *input) { 
    int result = 0x55555555;

    while (*input) { 
        result ^= *input++;
        result = rol(result, 5);
    }
}

/*
 * http://www.cse.yorku.ca/~oz/hash.html
 */
uint32_t
hash(uint8_t *str)
{
    uint32_t hash = 5381;
    uint8_t c;

    while (c = *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}

#endif

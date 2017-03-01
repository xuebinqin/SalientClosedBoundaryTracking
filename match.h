#ifndef  __MATCH_H
#define  __MATCH_H

#define CC_SAFE_MALLOC(nnum,type)                                          \
    (type *) CCutil_allocrus (((unsigned int) (nnum)) * sizeof (type))

#define CC_FREE(object,type) {                                             \
    CCutil_freerus ((void *) (object));                                    \
    object = (type *) NULL;                                                \
}

#define CC_IFFREE(object,type) {                                           \
    if ((object)) CC_FREE ((object),type);                                 \
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef CC_PROTOTYPE_ANSI

int
    min_ratio_cycle (int ncount, int ecount, int **elist,
             int **elen, int **elen2, int **esolid, int **grouped_edges, int NumIteration, int *sizeofgpedges);

#else

int
    min_ratio_cycle ();

#endif

#endif  /* __MATCH_H */

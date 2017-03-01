#define CC_PROTOTYPE_ANSI
extern "C" {
#include "match.h" 
}
#include <cstdio>
#include <vector>
#include <string>

using namespace std;

//static void usage (char *name);
//static int parseargs (int ac, char **av);

int CCutil_getedgelist_n2 (vector<double*> *graph, int ncount, int ecount, int **elist, int **elen, int **elen2, int **esolid, int **grouped_edges);


//static char match_file[256];
//static char edgefilename[256];
//static int nnodes_want = 0;


int RC(vector<double*> *graph, int ncount, int ecount, int NumIteration, vector<double*> &groupedEdges)
{
    int *elist = NULL;
    int *elen = NULL;
    int *elen2 = NULL;
    int *esolid = NULL;

    int **grouped_edges = new (int*);
    int sizeofgpedges;
    //snprintf(edgefilename, sizeof(edgefilename), "%s.w", fname.c_str());
    //snprintf(match_file, sizeof(match_file), "%s", fname.c_str()); /* .cycle is added in match.c */


    if (CCutil_getedgelist_n2 (graph, ncount, ecount,
                               &elist, &elen, &elen2, &esolid, grouped_edges)) {
       fprintf (stderr, "getedgelist_n2 failed\n");
       goto CLEANUP;
     }

    if (min_ratio_cycle (ncount, ecount, &elist, 
                         &elen, &elen2, &esolid, grouped_edges, NumIteration, &sizeofgpedges)){
        fprintf (stderr, "ratio-cycle finding failed\n");
        goto CLEANUP;
    }

    double* one_edge;
    for(int i = 0; i < sizeofgpedges;){
        one_edge = new double[6];
        one_edge[0] = (*grouped_edges)[i++];
        one_edge[1] = (*grouped_edges)[i++];
        one_edge[2] = (*grouped_edges)[i++];
        one_edge[3] = (*grouped_edges)[i++];
        one_edge[4] = (*grouped_edges)[i++];
        one_edge[5] = (*grouped_edges)[i++];
        groupedEdges.push_back(one_edge);
    }

CLEANUP:

   return 0;
}

/* function to read the edge file. elen -- 
   (first) edge weight; elen2 -- second edge weight;
   esolid -- 1 means solid edge and 0 means dashed edge */

int CCutil_getedgelist_n2(vector<double*> *graph, int ncount, int ecount, 
                                 int **elist, int **elen, int **elen2, int **esolid, int **grouped_edges)
{

    int i, k, t;
    double *edge;

    *elist = (int*)malloc(sizeof(int)*2*ecount);
    *elen = (int*)malloc(sizeof(int)*ecount);
    *elen2 = (int*)malloc(sizeof(int)*ecount);
    *esolid = (int*)malloc(sizeof(int)*ecount);

    *grouped_edges = (int*)malloc(sizeof(int)*6*ecount);
//new int[ecount];


    for (i = 0, k = 0, t=0; i < ecount; i++) {
        edge = (*graph)[i];
        (*elist)[k++] = (int)edge[0];
        (*elist)[k++] = (int)edge[1];
        (*elen)[i] = (int)edge[2]; 
        (*elen2)[i] = (int)edge[3];
        (*esolid)[i] = (int)edge[4];
        (*grouped_edges)[t++] = 0;
        (*grouped_edges)[t++] = 0;
        (*grouped_edges)[t++] = 0;
        (*grouped_edges)[t++] = 0;
        (*grouped_edges)[t++] = 0;
        (*grouped_edges)[t++] = 0;
    }

    return 0;
}

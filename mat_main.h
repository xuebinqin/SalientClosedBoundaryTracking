#ifndef MAT_MAIN_H
#define MAT_MAIN_H

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

int RC(vector<double*> *graph, int ncount, int ecount, int NumIteration, vector<double*> &groupedEdges);

#endif // MAT_MAIN_H

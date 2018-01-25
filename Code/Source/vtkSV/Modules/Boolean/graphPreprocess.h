#ifndef GRAPHPREPROCESS_H
#define GRAPHPREPROCESS_H

#include "graph.h"

#ifdef __cplusplus
extern "C" {
#endif

int    gp_CreateDFSTree(graphP theGraph);
int    gp_SortVertices(graphP theGraph);
void   gp_LowpointAndLeastAncestor(graphP theGraph);

#ifdef __cplusplus
}
#endif

#endif

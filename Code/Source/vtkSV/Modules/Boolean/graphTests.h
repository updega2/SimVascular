#ifndef GRAPHTESTS_H
#define GRAPHTESTS_H

#include "graph.h"
#include "stack.h"

/* Public functions for graphs */
#ifdef __cplusplus
extern "C" {
#endif

int    gp_CheckEmbeddingIntegrity(graphP theGraph);
int    gp_CheckKuratowskiSubgraphIntegrity(graphP theGraph);

/* Private function declarations */

int  _TestPath(graphP theGraph, int U, int V);
int  _TryPath(graphP theGraph, int J, int V);
void _MarkPath(graphP theGraph, int J);
int  _TestSubgraph(graphP theSubgraph, graphP theGraph);

#ifdef __cplusplus
}
#endif

#endif

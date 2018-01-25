#ifndef GRAPHEMBED_H
#define GRAPHEMBED_H

#include <stdlib.h>

#include "graph.h"

/* Imported functions */

#ifdef __cplusplus
extern "C" {
#endif

int    gp_Embed(graphP theGraph, int embedFlags);

extern void _InitGraphNode(graphP theGraph, int I);
extern void _FillVisitedFlags(graphP, int);

extern int _IsolateKuratowskiSubgraph(graphP theEmbedding, int I);

/* Private functions (some are exported to system only) */

void _CreateSortedSeparatedDFSChildLists(graphP theEmbedding);
void _CreateFwdArcLists(graphP theGraph);
void _CreateDFSTreeEmbedding(graphP theGraph);

void _EmbedBackEdgeToDescendant(graphP theEmbedding, int RootSide, int RootVertex, int W, int WPrevLink);

int  _GetNextVertexOnExternalFace(graphP theEmbedding, int curVertex, int *pPrevLink);

void _InvertVertex(graphP theEmbedding, int V);
void _MergeVertex(graphP theEmbedding, int W, int WPrevLink, int R);
void _MergeBicomps(graphP theEmbedding);

void _RecordPertinentChildBicomp(graphP theEmbedding, int I, int RootVertex);
#ifndef SPEED_MACROS
int  _GetPertinentChildBicomp(graphP theEmbedding, int W);
#endif

void _WalkUp(graphP theEmbedding, int I, int W);
void _WalkDown(graphP theEmbedding, int I, int RootVertex);

void _OrientVerticesInEmbedding(graphP theEmbedding);
void _OrientVerticesInBicomp(graphP theEmbedding, int BicompRoot, int PreserveSigns);
int  _JoinBicomps(graphP theEmbedding);

#ifdef __cplusplus
}
#endif

#endif

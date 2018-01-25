#ifndef GRAPHNONPLANAR_H
#define GRAPHNONPLANAR_H

#include "graph.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void _ClearIsolatorContext(graphP theGraph);
extern void _FillVisitedFlags(graphP, int);
extern void _FillVisitedFlagsInBicomp(graphP theGraph, int BicompRoot, int FillValue);
extern void _SetVertexTypeInBicomp(graphP theGraph, int BicompRoot, int theType);

extern int  _GetNextVertexOnExternalFace(graphP theEmbedding, int curVertex, int *pPrevLink);
extern int  _GetPertinentChildBicomp(graphP theEmbedding, int W);
extern void _WalkDown(graphP theEmbedding, int I, int RootVertex);
extern void _OrientVerticesInEmbedding(graphP theEmbedding);
extern void _OrientVerticesInBicomp(graphP theEmbedding, int BicompRoot, int PreserveSigns);

/* Private functions (exported to system) */

int  _ChooseTypeOfNonplanarityMinor(graphP theEmbedding, int I, int R);
int  _InitializeNonplanarityContext(graphP theEmbedding, int I, int R);

int  _FindNonplanarityBicompRoot(graphP theEmbedding);
void _FindActiveVertices(graphP theEmbedding, int R, int *pX, int *pY);
int  _FindPertinentVertex(graphP theEmbedding);

void _PopAndUnmarkVerticesAndEdges(graphP theEmbedding, int Z);

int  _MarkHighestXYPath(graphP theEmbedding);
int  _MarkZtoRPath(graphP theEmbedding);
int  _FindExtActivityBelowXYPath(graphP theEmbedding);

#ifdef __cplusplus
}
#endif

#endif

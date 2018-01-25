#ifndef GRAPHSTRUCTURE_H
#define GRAPHSTRUCTURE_H

#include <stdlib.h>

#include "graph.h"

#ifdef __cplusplus
extern "C" {
#endif

graphP gp_New(void);
int    gp_InitGraph(graphP theGraph, int N);
void   gp_ReinitializeGraph(graphP theGraph);
int    gp_CopyGraph(graphP dstGraph, graphP srcGraph);
graphP gp_DupGraph(graphP theGraph);

int    gp_CreateRandomGraph(graphP theGraph);
int    gp_CreateRandomGraphEx(graphP theGraph, int numEdges);

void   gp_Free(graphP *pGraph);

int    gp_Read(graphP theGraph, char *FileName);
int    gp_Write(graphP theGraph, char *FileName, int Mode);

int    gp_IsNeighbor(graphP theGraph, int u, int v);
int    gp_GetVertexDegree(graphP theGraph, int v);

int    gp_AddEdge(graphP theGraph, int u, int ulink, int v, int vlink);
void   gp_HideEdge(graphP theGraph, int arcPos);
void   gp_RestoreEdge(graphP theGraph, int arcPos);
int    gp_DeleteEdge(graphP theGraph, int J, int nextLink);

void _InitGraphNode(graphP theGraph, int I);
void _ClearIsolatorContext(graphP theGraph);
void _FillVisitedFlags(graphP theGraph, int FillValue);
void _FillVisitedFlagsInBicomp(graphP theGraph, int BicompRoot, int FillValue);
void _FillVisitedFlagsInOtherBicomps(graphP theGraph, int BicompRoot, int FillValue);
void _SetVertexTypeInBicomp(graphP theGraph, int BicompRoot, int theType);

void _ClearGraph(graphP theGraph);

void _InitVertexRec(graphP theGraph, int I);

int  _GetRandomNumber(int NMin, int NMax);

void _AddArc(graphP theGraph, int u, int v, int arcPos, int link);
void _HideArc(graphP theGraph, int arcPos);

int    gp_GetTwinArc(graphP theGraph, int Arc);

/********************************************************************
 int  gp_GetTwinArc(graphP theGraph, int Arc);
 This macro function returns the calculated twin arc of a given arc.
 If the arc location is even, then the successor is the twin.
 If the arc node is odd, then the predecessor is the twin.

 Logically, we return (Arc & 1) ? Arc-1 : Arc+1
 ********************************************************************/

#ifdef __cplusplus
}
#endif

#endif

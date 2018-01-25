#ifndef GRAPHISOLATOR_H
#define GRAPHISOLATOR_H

#include "graph.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void _FillVisitedFlags(graphP, int);

extern int  _GetNextVertexOnExternalFace(graphP theEmbedding, int curVertex, int *pPrevLink);
extern int  _JoinBicomps(graphP theEmbedding);
extern void _RecordPertinentChildBicomp(graphP theEmbedding, int I, int RootVertex);
extern int  _GetPertinentChildBicomp(graphP theEmbedding, int W);

extern int _ChooseTypeOfNonplanarityMinor(graphP theEmbedding, int I, int R);

/* Private function declarations (exported within system) */

int _IsolateKuratowskiSubgraph(graphP theEmbedding, int I);

int  _FindUnembeddedEdgeToAncestor(graphP theEmbedding, int cutVertex,
                                   int *pAncestor, int *pDescendant);
int  _FindUnembeddedEdgeToCurVertex(graphP theEmbedding, int cutVertex,
                                    int *pDescendant);
int  _FindUnembeddedEdgeToSubtree(graphP theEmbedding, int ancestor,
                                  int SubtreeRoot, int *pDescendant);

int  _MarkDFSPath(graphP theEmbedding, int ancestor, int descendant);
int  _MarkPathAlongBicompExtFace(graphP theEmbedding, int startVert, int endVert);

int  _AddAndMarkEdge(graphP theEmbedding, int ancestor, int descendant);
void _AddBackEdge(graphP theEmbedding, int ancestor, int descendant);
int  _DeleteUnmarkedVerticesAndEdges(graphP theEmbedding);

int  _InitializeIsolatorContext(graphP theEmbedding);

int  _IsolateMinorA(graphP theEmbedding);
int  _IsolateMinorB(graphP theEmbedding);
int  _IsolateMinorC(graphP theEmbedding);
int  _IsolateMinorD(graphP theEmbedding);
int  _IsolateMinorE(graphP theEmbedding);

int  _IsolateMinorE1(graphP theEmbedding);
int  _IsolateMinorE2(graphP theEmbedding);
int  _IsolateMinorE3(graphP theEmbedding);
int  _IsolateMinorE4(graphP theEmbedding);

int  _GetLeastAncestorConnection(graphP theEmbedding, int cutVertex);
int  _MarkDFSPathsToDescendants(graphP theEmbedding);
int  _AddAndMarkUnembeddedEdges(graphP theEmbedding);

#ifdef __cplusplus
}
#endif

#endif

#ifndef GRAPHIO_H
#define GRAPHIO_H

#include <stdlib.h>
#include <string.h>

#include "graph.h"

#ifdef __cplusplus
extern "C" {
#endif

int  _ReadAdjMatrix(graphP theGraph, FILE *Infile);
int  _ReadAdjList(graphP theGraph, FILE *Infile);
int  _WriteAdjList(graphP theGraph, FILE *Outfile);
int  _WriteAdjMatrix(graphP theGraph, FILE *Outfile);
int  _WriteDebugInfo(graphP theGraph, FILE *Outfile);

#ifdef __cplusplus
}
#endif

#endif


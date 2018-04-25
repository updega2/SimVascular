/*=========================================================================
 *
 * Copyright (c) 2014-2015 The Regents of the University of California.
 * All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *=========================================================================*/

/**
 *  \class vtkSVSurfaceCenterlineGrouper
 *  \brief Using a polydata centerlines, separate the polydata into regions
 *  based on the centerlines
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#ifndef vtkSVSurfaceCenterlineGrouper_h
#define vtkSVSurfaceCenterlineGrouper_h

#include "vtkSVSegmentationModule.h" // For export

#include "vtkIdList.h"
#include "vtkMatrix4x4.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyData.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"

#include "vtkSVPolyBallLine.h"
#include "vtkSVGlobals.h"

class VTKSVSEGMENTATION_EXPORT vtkSVSurfaceCenterlineGrouper : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkSVSurfaceCenterlineGrouper,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  static vtkSVSurfaceCenterlineGrouper *New();

  //@{
  /// \brief Get/Set macro for merged centerlines
  vtkSetObjectMacro(MergedCenterlines,vtkPolyData);
  vtkGetObjectMacro(MergedCenterlines,vtkPolyData);
  //@}

  //@{
  /// \brief Get/Set macro for surface polycube
  vtkSetObjectMacro(PolycubePd,vtkPolyData);
  vtkGetObjectMacro(PolycubePd,vtkPolyData);
  //@}

  //@{
  /// \brief Get/Set macro for array name used by the filter. Must
  //  be present on the centerlines.
  vtkSetStringMacro(CenterlineGroupIdsArrayName);
  vtkGetStringMacro(CenterlineGroupIdsArrayName);
  vtkSetStringMacro(CenterlineRadiusArrayName);
  vtkGetStringMacro(CenterlineRadiusArrayName);
  vtkSetStringMacro(GroupIdsArrayName);
  vtkGetStringMacro(GroupIdsArrayName);
  vtkSetStringMacro(BlankingArrayName);
  vtkGetStringMacro(BlankingArrayName);
  vtkSetStringMacro(CenterlineIdsArrayName);
  vtkGetStringMacro(CenterlineIdsArrayName);
  vtkSetStringMacro(TractIdsArrayName);
  vtkGetStringMacro(TractIdsArrayName);
  vtkSetStringMacro(PatchIdsArrayName);
  vtkGetStringMacro(PatchIdsArrayName);
  vtkSetStringMacro(SlicePointsArrayName);
  vtkGetStringMacro(SlicePointsArrayName);
  //@}

  //@{
  /// \brief Get/Set the radius information
  vtkSetMacro(UseRadiusInformation,int);
  vtkGetMacro(UseRadiusInformation,int);
  vtkBooleanMacro(UseRadiusInformation,int);
  //@}

  //@{
  /// \brief Get/Set the radius information
  vtkSetMacro(EnforcePolycubeConnectivity, int);
  vtkGetMacro(EnforcePolycubeConnectivity, int);
  vtkBooleanMacro(EnforcePolycubeConnectivity, int);
  //@}

  //@{
  /// \brief Get/Set the radius information
  vtkSetMacro(EnforceCenterlinesConnectivity, int);
  vtkGetMacro(EnforceCenterlinesConnectivity, int);
  vtkBooleanMacro(EnforceCenterlinesConnectivity, int);
  //@}

  //@{
  /// \brief Get/Set the radius information
  vtkSetMacro(GroupSurface, int);
  vtkGetMacro(GroupSurface, int);
  vtkBooleanMacro(GroupSurface, int);
  //@}


  /** \brief Correct cells on the boundary by updating val if they have
   *  multiple neighboring cells of the same value */
  static int CorrectCellBoundaries(vtkPolyData *pd, std::string cellArrayName);

  /** \brief Naive implementation to get most reoccuring number in list. Okay
   *  because list size is small. */
  static void GetMostOccuringVal(vtkIdList *idList, int &output, int &max_count);

  /** \brief Run basic edge weithing cvt with pd */
  static int RunEdgeWeightedCVT(vtkPolyData *pd, vtkPolyData *generatorPd);

  /** \brief From three vectors, compute transformation from global to local */
  static int ComputeRotationMatrix(const double vx[3], const double vy[3],
                                   const double vz[3], double rotMatrix[9]);

  static int SmoothBoundaries(vtkPolyData *pd, std::string arrayName);

  static int GetPointEdgeCells(vtkPolyData *pd, std::string arrayName,
                               const int cellId, const int pointId,
                               vtkIdList *sameCells);

  static int GetRegions(vtkPolyData *pd, std::string arrayName,
                        std::vector<Region> &allRegions);

  static int CurveFitBoundaries(vtkPolyData *pd, std::string arrayName,
                                std::vector<Region> allRegions);

  static int SplitBoundary(vtkPolyData *pd, std::vector<int> boundary,
                           int numDivs, int groupId, std::vector<int> &newSlicePoints, std::string slicePointsArrayName);

  static int GetCCWPoint(vtkPolyData *pd, const int pointId, const int cellId);
  static int GetCWPoint(vtkPolyData *pd, const int pointId, const int cellId);

  static int CheckBoundaryEdge(vtkPolyData *pd, std::string arrayName, const int cellId, const int pointId0, const int pointId1);

  static int CheckCellValuesEdge(vtkPolyData *pd, std::string arrayName, const int cellId, const int pointId0, const int pointId1);
  static void SplineKnots(std::vector<int> &u, int n, int t);

  static void SplineCurve(const std::vector<XYZ> &inp, int n, const std::vector<int> &knots, int t, std::vector<XYZ> &outp, int res);

  static void SplinePoint(const std::vector<int> &u, int n, int t, double v, const std::vector<XYZ> &control, XYZ &output);

  static double SplineBlend(int k, int t, const std::vector<int> &u, double v);

  static int FindPointMatchingValues(vtkPointSet *ps, std::string arrayName, vtkIdList *matchingVals, int &returnPtId);

  static int RotateGroupToGlobalAxis(vtkPolyData *pd,
                                     const int thresholdId,
                                     std::string arrayName,
                                     vtkPolyData *rotPd,
                                     vtkMatrix4x4 *rotMatrix0,
                                     vtkMatrix4x4 *rotMatrix1);
  static int InterpolateMapOntoTarget(vtkPolyData *sourceBasePd,
                                      vtkPolyData *targetPd,
                                      vtkPolyData *targetBasePd,
                                      vtkPolyData *mappedPd,
                                      std::string dataMatchingArrayName);
  static int GetCellRingNeighbors(vtkPolyData *pd,
                                  vtkIdList *cellIds,
                                  int ringNumber,
                                  int totNumberOfRings,
                                  std::vector<std::vector<int> > &neighbors);


  const static double GlobalCoords[3][3];


protected:
  vtkSVSurfaceCenterlineGrouper();
  ~vtkSVSurfaceCenterlineGrouper();

  // Usual data generation method
  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *) override;

  int PrepFilter(); // Prep work.
  int RunFilter(); // Run filter operations.

  int MergeCenterlines();
  int CheckPolycubeEnforcePossible();
  int MatchSurfaceToPolycube();
  int CheckSlicePoints();
  int SplitCellsAroundPoint(vtkPolyData *pd, int ptId);
  int SplitEdge(vtkPolyData *pd, int cellId, int ptId0, int ptId1,
                vtkCellArray *newCells, std::vector<std::vector<int> >  &splitCells);
  int FixMultipleGroups(vtkPolyData *pd, vtkPolyData *polycubePd,
                        std::vector<Region> surfaceGroups,
                        std::vector<Region> polycubeGroups);

  int CheckGroups(vtkPolyData *pd);
  int CheckGroups2();
  int FixEdges(vtkPolyData *pd, vtkPolyData *origPd, std::string arrayName,
               const Region region, std::vector<int> allEdges,
               std::vector<int> fixEdges, vtkIdList *critPts);
  int FixPlanarTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, std::string arrayName,
                            const Region region, std::vector<int> allEdges,
                            std::vector<int> badEdges, vtkIdList *critPts);
  int FixPerpenTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, std::string arrayName,
                            const Region region, std::vector<int> allEdges,
                            std::vector<int> badEdges, vtkIdList *critPts);
  int FixCornerTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, std::string arrayName,
                            const Region region, std::vector<int> allEdges,
                            std::vector<int> badEdges, vtkIdList *critPts);
  int FixOffsetTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, vtkPolyData *polyPd,
                            std::string arrayName,
                            const Region region, const Region polyRegion,
                            std::vector<int> allEdges,
                            std::vector<int> badEdges, vtkIdList *critPts);
  int FixFilledTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, vtkPolyData *polyPd,
                            std::string arrayName,
                            const Region region, const Region polyRegion,
                            std::vector<int> allEdges,
                            std::vector<int> badEdges, vtkIdList *critPts);
  int FixSplitsTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, vtkPolyData *polyPd,
                            std::string arrayName,
                            const Region region, const Region polyRegion,
                            std::vector<int> allEdges,
                            std::vector<int> badEdges, vtkIdList *critPts);
  int FixCloseGroup(vtkPolyData *pd, vtkPolyData *origPd, vtkPolyData *polyPd,
                    std::string arrayName,
                    const Region region, const Region polyRegion,
                    std::vector<int> allEdges,
                    std::vector<int> badEdges, vtkIdList *critPts);
  int FixGroupsWithPolycube();
  int FixGroupsWithCenterlines(int fixIters);

  int GetConnectedEdges(std::vector<std::vector<int> > inputEdges,
                        std::vector<std::vector<int> > &connectedCornerPts);

  int RemoveNegativeGroups(vtkPolyData *pd, std::string arrayName);
  int RemoveDuplicateGroups(vtkPolyData *pd, std::string arrayName);
  int FixRegions(vtkPolyData *pd, std::string arrayName,
                        std::vector<Region> &allRegions,
                        std::vector<int> badRegions,
                        const int currentValue,
                        const int onlyFixIslands = 0);

  char *CenterlineGroupIdsArrayName;
  char *CenterlineRadiusArrayName;
  char *CenterlineIdsArrayName;
  char *GroupIdsArrayName;
  char *BlankingArrayName;
  char *TractIdsArrayName;
  char *PatchIdsArrayName;
  char *SlicePointsArrayName;

  vtkPolyData *WorkPd;
  vtkPolyData *MergedCenterlines;
  vtkPolyData *PolycubePd;

  int EnforcePolycubeConnectivity;
  int EnforceCenterlinesConnectivity;
  int GroupSurface;
  int UseRadiusInformation;

private:
  vtkSVSurfaceCenterlineGrouper(const vtkSVSurfaceCenterlineGrouper&);  // Not implemented.
  void operator=(const vtkSVSurfaceCenterlineGrouper&);  // Not implemented.
};

#endif

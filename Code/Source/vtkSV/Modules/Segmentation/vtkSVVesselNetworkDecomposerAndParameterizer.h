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
 *  \class vtkSVVesselNetworkDecomposerAndParameterizer
 *  \brief Using a polydata centerlines, separate the polydata into regions
 *  based on the centerlines
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#ifndef vtkSVVesselNetworkDecomposerAndParameterizer_h
#define vtkSVVesselNetworkDecomposerAndParameterizer_h

#include "vtkSVSegmentationModule.h" // For export

#include "vtkIdList.h"
#include "vtkMatrix4x4.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyData.h"
#include "vtkSVPolyBallLine.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"

#include "vtkSVGlobals.h"

class VTKSVSEGMENTATION_EXPORT vtkSVVesselNetworkDecomposerAndParameterizer : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkSVVesselNetworkDecomposerAndParameterizer,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  static vtkSVVesselNetworkDecomposerAndParameterizer *New();

  //@{
  /// \brief Get/Set macro for the object's centerlines
  vtkSetObjectMacro(Centerlines,vtkPolyData);
  vtkGetObjectMacro(Centerlines,vtkPolyData);
  //@}

  //@{
  /// \brief Get the graph for the model
  vtkSetObjectMacro(GraphPd,vtkPolyData);
  vtkGetObjectMacro(GraphPd,vtkPolyData);
  //@}

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
  /// \brief Get/Set macro for surface polycube
  vtkSetObjectMacro(PolycubeUg,vtkUnstructuredGrid);
  vtkGetObjectMacro(PolycubeUg,vtkUnstructuredGrid);
  //@}

  //@{
  /// \brief Get/Set macro for surface polycube
  vtkSetObjectMacro(FinalHexMesh,vtkUnstructuredGrid);
  vtkGetObjectMacro(FinalHexMesh,vtkUnstructuredGrid);
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
  //@}

  //@{
  /// \brief Get/Set the cutoff radius factor for clipping of the surface
  //  distance functions
  vtkSetMacro(CutoffRadiusFactor,double);
  vtkGetMacro(CutoffRadiusFactor,double);
  //@}

  //@{
  /// \brief Get/Set the factor for enforcing of the boundary directions. Approximately represents the number of centerline points to enforce per branch. Default is 1, and typically a fairly low value works well. The larger the value, the larger the portion of the vessel is set explicitly, and sometimes this can cause large problems.
  vtkSetMacro(BoundaryEnforceFactor,int);
  vtkGetMacro(BoundaryEnforceFactor,int);
  //@}


  //@{
  /// \brief Get/Set the clip value for clipping of the surface distance functions.
  vtkSetMacro(ClipValue,double);
  vtkGetMacro(ClipValue,double);
  //@}

  //@{
  /// \brief Get/Set the radius information
  vtkSetMacro(UseRadiusInformation,int);
  vtkGetMacro(UseRadiusInformation,int);
  vtkBooleanMacro(UseRadiusInformation,int);
  //@}

  //@{
  /// \brief Get/Set the initial group clipper to use
  vtkSetMacro(UseVmtkClipping,int);
  vtkGetMacro(UseVmtkClipping,int);
  vtkBooleanMacro(UseVmtkClipping,int);
  //@}

  //@{
  /// \brief Get/Set whether the boundary at separating patches should be more
  //  strictly enforced.
  vtkSetMacro(EnforceBoundaryDirections,int);
  vtkGetMacro(EnforceBoundaryDirections,int);
  vtkBooleanMacro(EnforceBoundaryDirections,int);
  //@}

  //@{
  /// \brief Get/Set the number of divisions to use along width and height of polycube
  vtkSetMacro(PolycubeDivisions,int);
  vtkGetMacro(PolycubeDivisions,int);
  //@}

  //@{
  /// \brief Get/Set the unit length for each division of the polycube
  vtkSetMacro(PolycubeUnitLength,double);
  vtkGetMacro(PolycubeUnitLength,double);
  //@}

  //@{
  /// \brief Get/Set the scalar determing how much influence to put on the normal
  // of the cell and how much influence to put on the position of the cell for
  // the cube patch clustering.
  vtkSetMacro(NormalsWeighting,double);
  vtkGetMacro(NormalsWeighting,double);
  //@}

  //@{
  /// \brief Get/Set whether the model is a vascular model with artificial truncated
  //  boundaries
  vtkSetMacro(IsVasculature,int);
  vtkGetMacro(IsVasculature,int);
  vtkBooleanMacro(IsVasculature,int);
  //@}

  //@{
  /// \brief Get/Set If model is not vasculature, indicate how many centerline
  //  points to remove from the ends
  vtkSetMacro(NumberOfCenterlineRemovePts,int);
  vtkGetMacro(NumberOfCenterlineRemovePts,int);
  //@}
  //
  //@{
  /// \brief Get/Set to use absolute distance
  vtkSetMacro(UseAbsoluteMergeDistance, int);
  vtkGetMacro(UseAbsoluteMergeDistance, int);
  vtkBooleanMacro(UseAbsoluteMergeDistance, int);
  //@}

  //@{
  /// \brief Get/Set to use absolute distance
  vtkSetMacro(RadiusMergeRatio, double);
  vtkGetMacro(RadiusMergeRatio, double);
  //@}

  //@{
  /// \brief Get/Set to use absolute distance
  vtkSetMacro(MergeDistance, double);
  vtkGetMacro(MergeDistance, double);
  //@}

protected:
  vtkSVVesselNetworkDecomposerAndParameterizer();
  ~vtkSVVesselNetworkDecomposerAndParameterizer();

  // Usual data generation method
  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *) override;

  int PrepFilter(); // Prep work.
  int RunFilter(); // Run filter operations.

  int MergeCenterlines();

  char *CenterlineGroupIdsArrayName;
  char *CenterlineRadiusArrayName;
  char *GroupIdsArrayName;
  char *BlankingArrayName;

  vtkPolyData *WorkPd;
  vtkPolyData *GraphPd;
  vtkPolyData *Centerlines;
  vtkPolyData *MergedCenterlines;
  vtkPolyData *PolycubePd;

  vtkUnstructuredGrid *PolycubeUg;
  vtkUnstructuredGrid *FinalHexMesh;

  int UseRadiusInformation;
  int UseVmtkClipping;
  int EnforceBoundaryDirections;
  int IsVasculature;
  int NumberOfCenterlineRemovePts;
  int PolycubeDivisions;
  int BoundaryEnforceFactor;
  int UseAbsoluteMergeDistance;

  double CutoffRadiusFactor;
  double ClipValue;
  double PolycubeUnitLength;
  double NormalsWeighting;
  double RadiusMergeRatio;
  double MergeDistance;

private:
  vtkSVVesselNetworkDecomposerAndParameterizer(const vtkSVVesselNetworkDecomposerAndParameterizer&);  // Not implemented.
  void operator=(const vtkSVVesselNetworkDecomposerAndParameterizer&);  // Not implemented.
};

#endif

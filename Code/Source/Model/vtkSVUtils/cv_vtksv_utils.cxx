/* Copyright (c) Stanford University, The Regents of the University of
 *               California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
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
 */

/** @file cv_vtksv_utils.cxx
 *  @brief These functions are utilities that implement vtkvtksv classes
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "cv_vtksv_utils.h"

#include "SimVascular.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkDecimatePro.h"
#include "vtkDoubleArray.h"
#include "vtkGeometryFilter.h"
#include "vtkIntArray.h"
#include "vtkMath.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkSortDataArray.h"
#include "vtkThreshold.h"
#include "vtkTriangleFilter.h"

#include "vtkSVCenterlineParallelTransportVectors.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVNURBSUtils.h"
#include "vtkSVLoftNURBSSurface.h"
#include "vtkSVPolycubeGenerator.h"
#include "vtkSVParameterizeSurfaceOnPolycube.h"
#include "vtkSVParameterizeVolumeOnPolycube.h"
#include "vtkSVPassDataArray.h"
#include "vtkSVSurfaceCenterlineAttributesPasser.h"
#include "vtkSVSurfaceCenterlineGrouper.h"
#include "vtkSVSurfaceCuboidPatcher.h"

#include "cv_occtsolid_utils.h"

#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

int VTKSVUtils_DecomposePolyData( cvPolyData *polydata, cvPolyData *mergedCenterlines, cvPolyData **decomposedPolyData, std::vector<cvOCCTSolidModel*> &loftedSurfs)
{
  vtkPolyData *geom = polydata->GetVtkPolyData();
  vtkPolyData *centerlines = mergedCenterlines->GetVtkPolyData();
  cvPolyData *result1 = NULL;
  *decomposedPolyData = NULL;

  try {

    std::cout<<"Polycubing..."<<endl;

    vtkNew(vtkSVPolycubeGenerator, polycuber);
    polycuber->SetInputData(centerlines);
    polycuber->SetCenterlineGroupIdsArrayName("GroupIds");
    polycuber->SetCenterlineRadiusArrayName("MaximumInscribedSphereRadius");
    polycuber->SetGridIdsArrayName("GridIds");
    polycuber->SetPolycubeDivisions(5);
    polycuber->Update();

    if (polycuber->GetErrorCode() != 0)
    {
      std::cerr<<"Error creating polycube" << endl;
      return SV_ERROR;
    }

    vtkNew(vtkPolyData, polycubePd);
    vtkNew(vtkUnstructuredGrid, polycubeUg);
    vtkNew(vtkPolyData, graphPd);

    polycubePd->DeepCopy(polycuber->GetOutput());
    polycubeUg->DeepCopy(polycuber->GetVolumePolycubeUg());
    graphPd->DeepCopy(polycuber->GetGraphPd());

    // Generate normals just in case they don't exist
    vtkNew(vtkPolyDataNormals, normaler);
    normaler->SetInputData(geom);
    normaler->ComputePointNormalsOff();
    normaler->ComputeCellNormalsOn();
    normaler->SplittingOff();
    normaler->Update();

    geom->DeepCopy(normaler->GetOutput());
    geom->BuildLinks();
    vtkDataArray *normalsArray =
      geom->GetCellData()->GetArray("Normals");

    int stopCellNumber = ceil(geom->GetNumberOfCells()*0.0001);

    vtkNew(vtkSVSurfaceCenterlineGrouper, grouper);
    grouper->SetInputData(geom);
    grouper->SetPolycubePd(polycubePd);
    grouper->SetMergedCenterlines(centerlines);
    grouper->SetUseRadiusInformation(1);
    grouper->SetCenterlineRadiusArrayName("MaximumInscribedSphereRadius");
    grouper->SetCenterlineGroupIdsArrayName("GroupIds");
    grouper->SetCenterlineIdsArrayName("CenterlineIds");
    grouper->SetGroupIdsArrayName("GroupIds");
    grouper->SetTractIdsArrayName("TractIds");
    grouper->SetPatchIdsArrayName("PatchIds");
    grouper->SetSlicePointsArrayName("SlicePoints");
    grouper->GroupSurfaceOn();
    grouper->EnforceCenterlinesConnectivityOn();
    grouper->EnforcePolycubeConnectivityOn();
    grouper->DebugOn();
    grouper->Update();

    if (grouper->GetErrorCode() != 0)
    {
      std::cerr<<"Error grouping polydata around centerlines" << endl;
      return SV_ERROR;
    }

    geom->DeepCopy(grouper->GetOutput());

    vtkNew(vtkSVCenterlineParallelTransportVectors, parallelTransport);
    parallelTransport->SetInputData(centerlines);
    parallelTransport->SetGroupIdsArrayName("GroupIds");
    parallelTransport->SetParallelTransportVectorArrayName("ParallelTransportVector");
    parallelTransport->Update();

    if (parallelTransport->GetErrorCode() != 0)
    {
      std::cerr<<"Error creating a parallel transport frame along the centerlines" << endl;
      return SV_ERROR;
    }

    centerlines->DeepCopy(parallelTransport->GetOutput());

    vtkNew(vtkSVSurfaceCuboidPatcher, patcher);
    patcher->SetInputData(geom);
    patcher->SetPolycubePd(polycubePd);
    patcher->SetMergedCenterlines(centerlines);
    patcher->SetCenterlineRadiusArrayName("MaximumInscribedSphereRadius");
    patcher->SetCenterlineGroupIdsArrayName("GroupIds");
    patcher->SetCenterlineIdsArrayName("CenterlineIds");
    patcher->SetGroupIdsArrayName("GroupIds");
    patcher->SetTractIdsArrayName("TractIds");
    patcher->SetPatchIdsArrayName("PatchIds");
    patcher->SetSlicePointsArrayName("SlicePoints");
    patcher->SetClusteringVectorArrayName("NormalsTransformedToCenterlines");
    patcher->SetParallelTransportVectorArrayName("ParallelTransportVector");
    patcher->SetIsVasculature(1);
    patcher->SetNormalsWeighting(0.6);
    patcher->EnforcePolycubeConnectivityOn();
    patcher->Update();

    if (patcher->GetErrorCode() != 0)
    {
      std::cerr<<"Error creating nurbs patches on polydata" << endl;
      return SV_ERROR;
    }

    geom->DeepCopy(patcher->GetOutput());

    vtkNew(vtkSVParameterizeSurfaceOnPolycube, surfParameterizer);
    surfParameterizer->SetInputData(geom);
    surfParameterizer->SetPolycubePd(polycubePd);
    surfParameterizer->SetPolycubeUg(polycubeUg);
    surfParameterizer->SetGroupIdsArrayName("GroupIds");
    surfParameterizer->SetGridIdsArrayName("GridIds");
    surfParameterizer->Update();

    if (surfParameterizer->GetErrorCode() != 0)
    {
      std::cerr<<"Error parameterizing a final nurbs surface" << endl;
      return SV_ERROR;
    }

    // Get all the groups on the surface
    vtkNew(vtkIdList, groupIds);
    for (int i=0; i<polycubePd->GetNumberOfCells(); i++)
    {
      int groupVal = polycubePd->GetCellData()->GetArray(
          "GroupIds")->GetTuple1(i);
      groupIds->InsertUniqueId(groupVal);
    }
    vtkSortDataArray::Sort(groupIds);
    int numGroups = groupIds->GetNumberOfIds();

    // Get the divisions on the polycube
    vtkDataArray *polycubeDivisions = polycubeUg->GetFieldData()->GetArray("PolycubeDivisions");
    if (polycubeDivisions == NULL)
    {
      std::cerr<< "Array with name PolycubeDivivisions needs to be present on volume polycube" << endl;
      return SV_ERROR;
    }
    if (polycubeDivisions->GetNumberOfTuples() != numGroups ||
        polycubeDivisions->GetNumberOfComponents() != 4)
    {
      std::cerr<< "PolycubeDivisions array has " << polycubeDivisions->GetNumberOfTuples() << " tuples and  " << polycubeDivisions->GetNumberOfComponents() << ". Expected " << numGroups << " tuples, and 4 components" << endl;
      return SV_ERROR;
    }

    if (vtkSVGeneralUtils::CheckArrayExists(polycubeUg, 1, "GroupIds") != SV_OK)
    {
      std::cerr << "Group Ids Array with name GroupIds does not exist on volume polycube" << endl;
      return SV_ERROR;
    }

    if (vtkSVGeneralUtils::CheckArrayExists(polycubeUg, 0, "GridIds") != SV_OK)
    {
      std::cerr << "Grid point ids array with name GridIds does not exist on volume polycube" << endl;
      return SV_ERROR;
    }

    // Do a little smoothing
    vtkNew(vtkPolyData, helperPd);
    vtkSVGeneralUtils::ThresholdPd(surfParameterizer->GetPolycubeOnSurfacePd(), 0, 0, 0, "InteriorPoints", helperPd);

    vtkNew(vtkCleanPolyData, cleaner);
    cleaner->SetInputData(helperPd);
    cleaner->Update();

    vtkNew(vtkTriangleFilter, triangulator);
    triangulator->SetInputData(cleaner->GetOutput());
    triangulator->Update();

    helperPd->DeepCopy(triangulator->GetOutput());
    helperPd->BuildLinks();

    // Get regions of group ids, so that we can get critical points splitting groups
    std::vector<Region> groupRegions;
    if (vtkSVGeneralUtils::GetRegions(helperPd, "GroupIds", groupRegions) != SV_OK)
    {
      std::cerr << "Error gettring groups of input surface" << endl;
      return SV_ERROR;
    }

    vtkNew(vtkIdList, frontNeighbors);
    vtkNew(vtkIdList, backNeighbors);
    vtkNew(vtkIdList, frontGroupNeighbors);
    vtkNew(vtkIdList, backGroupNeighbors);
    for (int i=0; i<numGroups; i++)
    {
      int groupId = groupIds->GetId(i);

      //if (!(groupId == 2 || groupId == 6))
      //{
      //  continue;
      //}

      // Get centerline info
      vtkIdType nlinepts, *linepts;
      int centerlineId = centerlines->GetCellData()->GetArray("GroupIds")->LookupValue(groupId);
      centerlines->GetCellPoints(centerlineId, nlinepts, linepts);
      int isTerminating = 1;

      centerlines->GetPointCells(linepts[0], frontNeighbors);
      frontGroupNeighbors->Reset();
      for (int j=0; j<frontNeighbors->GetNumberOfIds(); j++)
      {
        frontGroupNeighbors->InsertNextId(centerlines->GetCellData()->GetArray(
          "GroupIds")->GetTuple1(frontNeighbors->GetId(j)));
      }

      centerlines->GetPointCells(linepts[nlinepts-1], backNeighbors);
      backGroupNeighbors->Reset();
      for (int j=0; j<backNeighbors->GetNumberOfIds(); j++)
      {
        backGroupNeighbors->InsertNextId(centerlines->GetCellData()->GetArray(
          "GroupIds")->GetTuple1(backNeighbors->GetId(j)));
      }

      if (backNeighbors->GetNumberOfIds() != 1 && frontNeighbors->GetNumberOfIds() != 1)
        isTerminating = 0;

      double whl_divs[4];
      for (int j=0; j<4; j++)
        whl_divs[j] = -1.0;
      for (int j=0; j<polycubeDivisions->GetNumberOfTuples(); j++)
      {
        polycubeDivisions->GetTuple(j, whl_divs);
        if (whl_divs[0] == groupId)
          break;
      }
      if (whl_divs[0] == -1.0)
      {
        std::cerr<< "Field data array PolycubeDivisions did not have divisions for group number " << groupId << endl;
        return SV_ERROR;
      }

      std::cout << "SETTING UP AND LOFTING GROUP: " << groupId << " DIMS: " << whl_divs[1] << " " << whl_divs[2] << " " << whl_divs[3] << endl;

      // Threshold out each group
      vtkNew(vtkPolyData, thresholdPd);
      thresholdPd->DeepCopy(surfParameterizer->GetPolycubeOnSurfacePd());
      vtkSVGeneralUtils::ThresholdPd(thresholdPd, groupId, groupId, 1, "GroupIds");

      vtkDataArray *ptIds = thresholdPd->GetPointData()->GetArray("GridIds");

      int dim[3];
      dim[0] = whl_divs[1]; dim[1] = whl_divs[2]; dim[2] = whl_divs[3];

      // Set up new structured grid for nurbs surface lofting
      int newDims[3];
      newDims[0] = whl_divs[3]; newDims[1] = 2*whl_divs[1] + 2*whl_divs[2] - 3; newDims[2] = 1;
      vtkNew(vtkStructuredGrid, inputGrid);
      vtkNew(vtkPoints, inputGridPoints);
      inputGridPoints->SetNumberOfPoints(newDims[0] * newDims[1]);
      inputGrid->SetPoints(inputGridPoints);
      inputGrid->SetDimensions(newDims);

      // Loop through length and get exterior
      int rowCount = 0;
      int pos[3], newPos[3];
      int ptId, realPtId;
      double pt[3];
      for (int j=0; j<newDims[0]; j++)
      {
        rowCount = 0;
        // Go along bottom edge
        for (int k=0; k<dim[0]; k++)
        {
          pos[0] = k; pos[1] = 0; pos[2] = j;
          ptId = vtkStructuredData::ComputePointId(dim, pos);

          realPtId = ptIds->LookupValue(ptId);

          if (realPtId == -1)
          {
            std::cerr << "BOTTOM EDGE DIDN'T WORK" << endl;
            return SV_ERROR;
          }
          else
          {
            thresholdPd->GetPoint(realPtId, pt);

            newPos[0] = j; newPos[1] = rowCount++; newPos[2] = 0;
            ptId = vtkStructuredData::ComputePointId(newDims, newPos);
            inputGrid->GetPoints()->SetPoint(ptId, pt);
          }
        }
        // Go along right edge
        for (int k=1; k<dim[1]; k++)
        {
          pos[0] = dim[0]-1; pos[1] = k; pos[2] = j;
          ptId = vtkStructuredData::ComputePointId(dim, pos);

          realPtId = ptIds->LookupValue(ptId);

          if (realPtId == -1)
          {
            std::cerr << "RIGHT EDGE DIDN'T WORK" << endl;
            return SV_ERROR;
          }
          else
          {
            thresholdPd->GetPoint(realPtId, pt);

            newPos[0] = j; newPos[1] = rowCount++; newPos[2] = 0;
            ptId = vtkStructuredData::ComputePointId(newDims, newPos);
            inputGrid->GetPoints()->SetPoint(ptId, pt);
          }
        }
        // Go along top edge
        for (int k=1; k<dim[0]; k++)
        {
          pos[0] = dim[0]-k-1; pos[1] = dim[1]-1; pos[2] = j;
          ptId = vtkStructuredData::ComputePointId(dim, pos);

          realPtId = ptIds->LookupValue(ptId);

          if (realPtId == -1)
          {
            std::cerr << "TOP EDGE DIDN'T WORK" << endl;
            return SV_ERROR;
          }
          else
          {
            thresholdPd->GetPoint(realPtId, pt);

            newPos[0] = j; newPos[1] = rowCount++; newPos[2] = 0;
            ptId = vtkStructuredData::ComputePointId(newDims, newPos);
            inputGrid->GetPoints()->SetPoint(ptId, pt);
          }
        }

        // Go along left edge
        for (int k=1; k<dim[1]; k++)
        {
          pos[0] = 0; pos[1] = dim[1]-k-1; pos[2] = j;
          ptId = vtkStructuredData::ComputePointId(dim, pos);

          realPtId = ptIds->LookupValue(ptId);

          if (realPtId == -1)
          {
            std::cerr << "LEFT EDGE DIDN'T WORK" << endl;
            return SV_ERROR;
          }
          else
          {
            thresholdPd->GetPoint(realPtId, pt);

            newPos[0] = j; newPos[1] = rowCount++; newPos[2] = 0;
            ptId = vtkStructuredData::ComputePointId(newDims, newPos);
            inputGrid->GetPoints()->SetPoint(ptId, pt);
          }
        }
      }

      //// Now loft each surface
      //int uDegree = 3;
      //int vDegree = 3;
      //vtkNew(vtkSVLoftNURBSSurface,lofter);
      //lofter->SetInputData(inputGrid);
      //lofter->SetUDegree(uDegree);
      //lofter->SetVDegree(vDegree);
      //lofter->SetPolyDataUSpacing(0.1);
      //lofter->SetPolyDataVSpacing(0.1);
      //lofter->SetUKnotSpanType("average");
      //lofter->SetVKnotSpanType("average");
      //lofter->SetUParametricSpanType("chord");
      //lofter->SetVParametricSpanType("chord");
      //lofter->Update();

      //vtkNew(vtkSVNURBSSurface, NURBSSurface);
      //NURBSSurface->DeepCopy(lofter->GetSurface());

      // See what looks like by using as control grid

      int uDegree = 3;
      int vDegree = 3;
      std::string putype = "chord";
      std::string pvtype = "chord";
      std::string kutype = "average";
      std::string kvtype = "average";

      // Set the temporary control points
      vtkNew(vtkPoints, tmpUPoints);
      tmpUPoints->SetNumberOfPoints(newDims[0]);
      for (int k=0; k<newDims[0]; k++)
      {
        int pos[3]; pos[0] = k; pos[1] = 0; pos[2] = 0;
        int ptId = vtkStructuredData::ComputePointId(newDims, pos);
        tmpUPoints->SetPoint(k, inputGrid->GetPoint(ptId));
      }

      // Get the input point set u representation
      vtkNew(vtkDoubleArray, U);
      if (vtkSVNURBSUtils::GetUs(tmpUPoints, putype, U) != SV_OK)
      {
        return SV_ERROR;
      }

      // Get the knots in the u direction
      vtkNew(vtkDoubleArray, uKnots);
      if (vtkSVNURBSUtils::GetKnots(U, uDegree, kutype, uKnots) != SV_OK)
      {
        std::cerr<<"Error getting knots"<<endl;
        return SV_ERROR;
      }
      //
      vtkNew(vtkPoints, tmpVPoints);
      tmpVPoints->SetNumberOfPoints(newDims[1]);
      for (int k=0; k<newDims[1]; k++)
      {
        int pos[3]; pos[0] = 0; pos[1] = k; pos[2] = 0;
        int ptId = vtkStructuredData::ComputePointId(dim, pos);
        tmpVPoints->SetPoint(k, inputGrid->GetPoint(ptId));
      }
      // Get the input point set v representation
      vtkNew(vtkDoubleArray, V);
      if (vtkSVNURBSUtils::GetUs(tmpVPoints, pvtype, V) != SV_OK)
      {
        return SV_ERROR;
      }

      // Get the knots in the v direction
      vtkNew(vtkDoubleArray, vKnots);
      if (vtkSVNURBSUtils::GetKnots(V, vDegree, kvtype, vKnots) != SV_OK)
      {
        std::cerr<<"Error getting knots"<<endl;
        return SV_ERROR;
      }

      vtkNew(vtkSVNURBSSurface, NURBSSurface);
      NURBSSurface->SetKnotVector(uKnots, 0);
      NURBSSurface->SetKnotVector(vKnots, 1);
      NURBSSurface->SetControlPoints(inputGrid);
      NURBSSurface->SetUDegree(uDegree);
      NURBSSurface->SetVDegree(vDegree);

      // Get multiplicities, then convert everything to double arrays
      vtkNew(vtkDoubleArray, USingleKnotArray);
      vtkNew(vtkIntArray,    UMultArray);
      NURBSSurface->GetUMultiplicity(UMultArray, USingleKnotArray);
      vtkNew(vtkDoubleArray, VSingleKnotArray);
      vtkNew(vtkIntArray,    VMultArray);
      NURBSSurface->GetVMultiplicity(VMultArray, VSingleKnotArray);
      vtkSVControlGrid *controlPointGrid = NURBSSurface->GetControlPointGrid();

      // Get all information needed by creation of bspline surface
      int dims[3];
      controlPointGrid->GetDimensions(dims);
      int Xlen1 = dims[0];
      int Xlen2 = dims[1];
      double **Xarr = new double*[Xlen1];
      double **Yarr = new double*[Xlen1];
      double **Zarr = new double*[Xlen1];
      for (int j=0; j<Xlen1; j++)
      {
        Xarr[j] = new double[Xlen2];
        Yarr[j] = new double[Xlen2];
        Zarr[j] = new double[Xlen2];
        for (int k=0; k<Xlen2; k++)
        {
          double pt[3];
          double w;
          controlPointGrid->GetControlPoint(j, k, 0, pt, w);
          Xarr[j][k] = pt[0];
          Yarr[j][k] = pt[1];
          Zarr[j][k] = pt[2];
        }
      }

      // Get mult information as arrays
      int uKlen = USingleKnotArray->GetNumberOfTuples();
      double *uKarr = new double[uKlen];
      for (int j=0; j<uKlen; j++)
        uKarr[j] = USingleKnotArray->GetTuple1(j);

      int vKlen = VSingleKnotArray->GetNumberOfTuples();
      double *vKarr = new double[vKlen];
      for (int j=0; j<vKlen; j++)
        vKarr[j] = VSingleKnotArray->GetTuple1(j);

      int uMlen = UMultArray->GetNumberOfTuples();
      double *uMarr = new double[uMlen];
      for (int j=0; j<uMlen; j++)
        uMarr[j] = UMultArray->GetTuple1(j);

      int vMlen = VMultArray->GetNumberOfTuples();
      double *vMarr = new double[vMlen];
      for (int j=0; j<vMlen; j++)
        vMarr[j] = VMultArray->GetTuple1(j);

      // Flipping order!
      cvOCCTSolidModel* surf=new cvOCCTSolidModel();
      if (surf->CreateBSplineSurface(Xarr, Yarr, Zarr,
                                     Xlen1, Xlen2,
                                     vKarr, vKlen,
                                     uKarr, uKlen,
                                     vMarr, vMlen,
                                     uMarr, uMlen,
                                     vDegree, uDegree) != SV_OK )
      {
        std::cerr << "poly manipulation error ";
          delete surf;
          //Clean up
          for (int j=0;j<Xlen1;j++)
          {
            delete [] Xarr[j];
            delete [] Yarr[j];
            delete [] Zarr[j];
          }
          delete [] Xarr;
          delete [] Yarr;
          delete [] Zarr;

          delete [] uKarr;
          delete [] vKarr;
          delete [] uMarr;
          delete [] vMarr;
          return NULL;
      }
      //Clean up
      for (int j=0;j<Xlen1;j++)
      {
        delete [] Xarr[j];
        delete [] Yarr[j];
        delete [] Zarr[j];
      }
      delete [] Xarr;
      delete [] Yarr;
      delete [] Zarr;

      delete [] uKarr;
      delete [] vKarr;
      delete [] uMarr;
      delete [] vMarr;

      // Now we need to get the points to split at
      for (int surfEnd = 0; surfEnd < 2; surfEnd++)
      {
        vtkNew(vtkPoints, cartesianPoints);
        for (int j=0; j<groupRegions.size(); j++)
        {
          if (groupRegions[j].IndexCluster != groupId)
          {
            continue;
          }

          fprintf(stdout,"NUM CORNER POINTS: %d\n", groupRegions[j].CornerPoints.size());
          int cornerPtId;
          double criticalPt[3];
          vtkNew(vtkIdList, pointGroupIds);
          vtkNew(vtkIdList, intersectList);
          for (int k=0; k<groupRegions[j].CornerPoints.size(); k++)
          {
            cornerPtId = groupRegions[j].CornerPoints[k];
            helperPd->GetPoint(cornerPtId, criticalPt);
            //geom->GetPoint(cornerPtId, criticalPt);

            pointGroupIds->Reset();
            vtkSVGeneralUtils::GetPointCellsValues(helperPd, "GroupIds", cornerPtId, pointGroupIds);
            //vtkSVGeneralUtils::GetPointCellsValues(geom, "GroupIds", cornerPtId, pointGroupIds);

            intersectList->Reset();
            intersectList->DeepCopy(pointGroupIds);

            if (surfEnd == 0)

            {
              intersectList->IntersectWith(frontGroupNeighbors);

              if (frontGroupNeighbors->GetNumberOfIds() == intersectList->GetNumberOfIds() &&
                  pointGroupIds->GetNumberOfIds() == intersectList->GetNumberOfIds())
              {
                cartesianPoints->InsertNextPoint(criticalPt);
              }
            }
            else
            {
              intersectList->IntersectWith(backGroupNeighbors);
              if (backGroupNeighbors->GetNumberOfIds() == intersectList->GetNumberOfIds() &&
                  pointGroupIds->GetNumberOfIds() == intersectList->GetNumberOfIds())
              {
                cartesianPoints->InsertNextPoint(criticalPt);
              }
            }
          }
        }

        int edgeNumber;
        if (groupId == 0)
        {
          edgeNumber = surfEnd;
        }
        else
        {
          edgeNumber = (surfEnd+1)%2;
        }

        // Must split edges at common points so that we can actually merge shapes
        fprintf(stdout,"SPLITTING NUM POINTS: %d\n", cartesianPoints->GetNumberOfPoints());
        if (cartesianPoints->GetNumberOfPoints() >= 2 && cartesianPoints->GetNumberOfPoints() <= 6)
        {
          if (OCCTUtils_SplitEdgeOnShapeNearPoints(*(surf->geom_), edgeNumber, cartesianPoints) != SV_OK)
          {
            std::cerr << "Could not split edge on surface" << endl;
            return NULL;
          }
        }
      }

      loftedSurfs.push_back(surf);
    }

    // Now going to do caps as well
    for (int i=0; i<numGroups; i++)
    {
      int groupId = groupIds->GetId(i);

      double whl_divs[4];
      for (int j=0; j<4; j++)
        whl_divs[j] = -1.0;
      for (int j=0; j<polycubeDivisions->GetNumberOfTuples(); j++)
      {
        polycubeDivisions->GetTuple(j, whl_divs);
        if (whl_divs[0] == groupId)
          break;
      }
      if (whl_divs[0] == -1.0)
      {
        std::cerr<< "Field data array PolycubeDivisions did not have divisions for group number " << groupId << endl;
        return SV_ERROR;
      }

      std::cout << "SETTING UP AND LOFTING GROUP: " << groupId << " DIMS: " << whl_divs[1] << " " << whl_divs[2] << " " << whl_divs[3] << endl;

      // Threshold out each group
      vtkNew(vtkPolyData, thresholdPd);
      thresholdPd->DeepCopy(surfParameterizer->GetPolycubeOnSurfacePd());
      vtkSVGeneralUtils::ThresholdPd(thresholdPd, groupId, groupId, 1, "GroupIds");

      vtkNew(vtkPolyData, groupPolycubePd);
      groupPolycubePd->DeepCopy(polycubePd);
      vtkSVGeneralUtils::ThresholdPd(groupPolycubePd, groupId, groupId, 1, "GroupIds");
      if (groupPolycubePd->GetNumberOfCells() == 4)
      {
        continue;
      }

      vtkDataArray *ptIds = thresholdPd->GetPointData()->GetArray("GridIds");

      int dim[3];
      dim[0] = whl_divs[1]; dim[1] = whl_divs[2]; dim[2] = whl_divs[3];

      // Set up new structured grid for nurbs surface lofting
      int newDims[3];
      newDims[0] = whl_divs[1]; newDims[1] = whl_divs[2]; newDims[2] = 1;

      vtkNew(vtkStructuredGrid, inputGrid0);
      vtkNew(vtkPoints, inputGridPoints0);
      inputGridPoints0->SetNumberOfPoints(newDims[0] * newDims[1]);
      inputGrid0->SetPoints(inputGridPoints0);
      inputGrid0->SetDimensions(newDims);

      vtkNew(vtkStructuredGrid, inputGrid1);
      vtkNew(vtkPoints, inputGridPoints1);
      inputGridPoints1->SetNumberOfPoints(newDims[0] * newDims[1]);
      inputGrid1->SetPoints(inputGridPoints1);
      inputGrid1->SetDimensions(newDims);

      // Loop through length and get exterior
      int pos[3], newPos[3];
      int ptId, realPtId;
      double pt[3];
      for (int j=0; j<newDims[0]; j++)
      {
        // bottom cap
        for (int k=0; k<newDims[1]; k++)
        {
          pos[0] = j; pos[1] = k; pos[2] = 0;
          ptId = vtkStructuredData::ComputePointId(dim, pos);

          realPtId = ptIds->LookupValue(ptId);

          if (realPtId == -1)
          {
            std::cerr << "BOTTOM CAP DIDN'T WORK" << endl;
            return SV_ERROR;
          }
          else
          {
            thresholdPd->GetPoint(realPtId, pt);

            newPos[0] = j; newPos[1] = k; newPos[2] = 0;
            ptId = vtkStructuredData::ComputePointId(newDims, newPos);
            inputGrid0->GetPoints()->SetPoint(ptId, pt);
          }
        }

        // top cap
        for (int k=0; k<newDims[1]; k++)
        {
          pos[0] = j; pos[1] = k; pos[2] = dim[2]-1;
          ptId = vtkStructuredData::ComputePointId(dim, pos);

          realPtId = ptIds->LookupValue(ptId);

          if (realPtId == -1)
          {
            std::cerr << "TOP CAP DIDN'T WORK" << endl;
            return SV_ERROR;
          }
          else
          {
            thresholdPd->GetPoint(realPtId, pt);

            newPos[0] = j; newPos[1] = k; newPos[2] = 0;
            ptId = vtkStructuredData::ComputePointId(newDims, newPos);
            inputGrid1->GetPoints()->SetPoint(ptId, pt);
          }
        }
      }

      //// Now loft each surface
      //int uDegree = 2;
      //int vDegree = 2;
      //vtkNew(vtkSVLoftNURBSSurface,lofter0);
      //lofter0->SetInputData(inputGrid0);
      //lofter0->SetUDegree(uDegree);
      //lofter0->SetVDegree(vDegree);
      //lofter0->SetPolyDataUSpacing(0.1);
      //lofter0->SetPolyDataVSpacing(0.1);
      //lofter0->SetUKnotSpanType("average");
      //lofter0->SetVKnotSpanType("average");
      //lofter0->SetUParametricSpanType("chord");
      //lofter0->SetVParametricSpanType("chord");
      //lofter0->Update();
      //
      //vtkNew(vtkSVNURBSSurface, NURBSSurface0);
      //NURBSSurface0->DeepCopy(lofter0->GetSurface());

      // Try as control grid
      int uDegree = 2;
      int vDegree = 2;
      std::string putype = "chord";
      std::string pvtype = "chord";
      std::string kutype = "average";
      std::string kvtype = "average";

      // Set the temporary control points
      vtkNew(vtkPoints, tmpUPoints0);
      tmpUPoints0->SetNumberOfPoints(newDims[0]);
      for (int k=0; k<newDims[0]; k++)
      {
        int pos[3]; pos[0] = k; pos[1] = 0; pos[2] = 0;
        int ptId = vtkStructuredData::ComputePointId(newDims, pos);
        tmpUPoints0->SetPoint(k, inputGrid0->GetPoint(ptId));
      }

      // Get the input point set u representation
      vtkNew(vtkDoubleArray, U0);
      if (vtkSVNURBSUtils::GetUs(tmpUPoints0, putype, U0) != SV_OK)
      {
        return SV_ERROR;
      }

      // Get the knots in the u direction
      vtkNew(vtkDoubleArray, uKnots0);
      if (vtkSVNURBSUtils::GetKnots(U0, uDegree, kutype, uKnots0) != SV_OK)
      {
        std::cerr<<"Error getting knots"<<endl;
        return SV_ERROR;
      }
      //
      vtkNew(vtkPoints, tmpVPoints0);
      tmpVPoints0->SetNumberOfPoints(newDims[1]);
      for (int k=0; k<newDims[1]; k++)
      {
        int pos[3]; pos[0] = 0; pos[1] = k; pos[2] = 0;
        int ptId = vtkStructuredData::ComputePointId(dim, pos);
        tmpVPoints0->SetPoint(k, inputGrid0->GetPoint(ptId));
      }
      // Get the input point set v representation
      vtkNew(vtkDoubleArray, V0);
      if (vtkSVNURBSUtils::GetUs(tmpVPoints0, pvtype, V0) != SV_OK)
      {
        return SV_ERROR;
      }

      // Get the knots in the v direction
      vtkNew(vtkDoubleArray, vKnots0);
      if (vtkSVNURBSUtils::GetKnots(V0, vDegree, kvtype, vKnots0) != SV_OK)
      {
        std::cerr<<"Error getting knots"<<endl;
        return SV_ERROR;
      }

      vtkNew(vtkSVNURBSSurface, NURBSSurface0);
      NURBSSurface0->SetKnotVector(uKnots0, 0);
      NURBSSurface0->SetKnotVector(vKnots0, 1);
      NURBSSurface0->SetControlPoints(inputGrid0);
      NURBSSurface0->SetUDegree(uDegree);
      NURBSSurface0->SetVDegree(vDegree);

      //vtkNew(vtkSVLoftNURBSSurface,lofter1);
      //lofter1->SetInputData(inputGrid1);
      //lofter1->SetUDegree(uDegree);
      //lofter1->SetVDegree(vDegree);
      //lofter1->SetPolyDataUSpacing(0.1);
      //lofter1->SetPolyDataVSpacing(0.1);
      //lofter1->SetUKnotSpanType("average");
      //lofter1->SetVKnotSpanType("average");
      //lofter1->SetUParametricSpanType("chord");
      //lofter1->SetVParametricSpanType("chord");
      //lofter1->Update();

      //vtkNew(vtkSVNURBSSurface, NURBSSurface1);
      //NURBSSurface1->DeepCopy(lofter1->GetSurface());

      // Set the temporary control points
      vtkNew(vtkPoints, tmpUPoints1);
      tmpUPoints1->SetNumberOfPoints(newDims[0]);
      for (int k=0; k<newDims[0]; k++)
      {
        int pos[3]; pos[0] = k; pos[1] = 0; pos[2] = 0;
        int ptId = vtkStructuredData::ComputePointId(newDims, pos);
        tmpUPoints1->SetPoint(k, inputGrid1->GetPoint(ptId));
      }

      // Get the input point set u representation
      vtkNew(vtkDoubleArray, U1);
      if (vtkSVNURBSUtils::GetUs(tmpUPoints1, putype, U1) != SV_OK)
      {
        return SV_ERROR;
      }

      // Get the knots in the u direction
      vtkNew(vtkDoubleArray, uKnots1);
      if (vtkSVNURBSUtils::GetKnots(U1, uDegree, kutype, uKnots1) != SV_OK)
      {
        std::cerr<<"Error getting knots"<<endl;
        return SV_ERROR;
      }
      //
      vtkNew(vtkPoints, tmpVPoints1);
      tmpVPoints1->SetNumberOfPoints(newDims[1]);
      for (int k=0; k<newDims[1]; k++)
      {
        int pos[3]; pos[0] = 0; pos[1] = k; pos[2] = 0;
        int ptId = vtkStructuredData::ComputePointId(dim, pos);
        tmpVPoints1->SetPoint(k, inputGrid1->GetPoint(ptId));
      }
      // Get the input point set v representation
      vtkNew(vtkDoubleArray, V1);
      if (vtkSVNURBSUtils::GetUs(tmpVPoints1, pvtype, V1) != SV_OK)
      {
        return SV_ERROR;
      }

      // Get the knots in the v direction
      vtkNew(vtkDoubleArray, vKnots1);
      if (vtkSVNURBSUtils::GetKnots(V1, vDegree, kvtype, vKnots1) != SV_OK)
      {
        std::cerr<<"Error getting knots"<<endl;
        return SV_ERROR;
      }

      std::cout <<"U KNOTS"<<endl;
      vtkSVNURBSUtils::PrintArray(uKnots1);
      std::cout <<"V KNOTS"<<endl;
      vtkSVNURBSUtils::PrintArray(vKnots1);
      std::cout<<"CONTROL GRID" << endl;
      vtkSVNURBSUtils::PrintStructuredGrid(inputGrid1);
      std::cout<<"U DEGREE " << uDegree << endl;
      std::cout<<"V DEGREE " << vDegree << endl;

      vtkNew(vtkSVNURBSSurface, NURBSSurface1);
      NURBSSurface1->SetKnotVector(uKnots1, 0);
      NURBSSurface1->SetKnotVector(vKnots1, 1);
      NURBSSurface1->SetControlPoints(inputGrid1);
      NURBSSurface1->SetUDegree(uDegree);
      NURBSSurface1->SetVDegree(vDegree);

      std::cout <<"AFTER!!" << endl;
      std::cout <<"U KNOTS"<<endl;
      vtkSVNURBSUtils::PrintArray(NURBSSurface1->GetUKnotVector());
      std::cout <<"V KNOTS"<<endl;
      vtkSVNURBSUtils::PrintArray(NURBSSurface1->GetVKnotVector());
      std::cout<<"CONTROL GRID" << endl;
      vtkSVNURBSUtils::PrintStructuredGrid(NURBSSurface1->GetControlPointGrid());
      std::cout<<"U DEGREE " << NURBSSurface1->GetUDegree() << endl;
      std::cout<<"V DEGREE " << NURBSSurface1->GetVDegree() << endl;

      if (groupId == 0)
      {
        NURBSSurface1->GeneratePolyDataRepresentation(0.1, 0.1);
        vtkSVIOUtils::WriteVTPFile("/Users/adamupdegrove/Desktop/tmp/CAp.vtp", NURBSSurface1->GetSurfaceRepresentation());
        cvOCCTSolidModel* surf1=new cvOCCTSolidModel();
        if (VTKSVUtils_NURBSSurfaceToOCCTBSpline(NURBSSurface1, surf1) != SV_OK)
        {
          std::cerr << "Error converting cap to occt model" << endl;
          return SV_ERROR;
        }

        loftedSurfs.push_back(surf1);

        if (groupPolycubePd->GetNumberOfCells() == 6)
        {

          cvOCCTSolidModel* surf0=new cvOCCTSolidModel();
          if (VTKSVUtils_NURBSSurfaceToOCCTBSpline(NURBSSurface0, surf0) != SV_OK)
          {
            std::cerr << "Error converting cap to occt model" << endl;
            return SV_ERROR;
          }

          loftedSurfs.push_back(surf0);
        }
      }
      else
      {
        cvOCCTSolidModel* surf0=new cvOCCTSolidModel();
        if (VTKSVUtils_NURBSSurfaceToOCCTBSpline(NURBSSurface0, surf0) != SV_OK)
        {
          std::cerr << "Error converting cap to occt model" << endl;
          return SV_ERROR;
        }

        loftedSurfs.push_back(surf0);
      }
    }

    //vtkNew(vtkSVParameterizeVolumeOnPolycube, volParameterizer);
    //volParameterizer->SetInputData(geom);
    //volParameterizer->SetPolycubeUg(polycubeUg);
    //volParameterizer->SetSurfaceOnPolycubePd(surfParameterizer->GetOutput());
    //volParameterizer->SetGroupIdsArrayName("GroupIds");
    //volParameterizer->Update();

    result1 = new cvPolyData( patcher->GetOutput() );
    *decomposedPolyData = result1;
  }
  catch (...) {
    fprintf(stderr,"ERROR in surface creation.\n");
    fflush(stderr);
    return SV_ERROR;
  }

  return SV_OK;
}

int VTKSVUtils_NURBSSurfaceToOCCTBSpline(vtkSVNURBSSurface *NURBSSurface,
                                         cvOCCTSolidModel *solidModel)
{
  // Get multiplicities, then convert everything to double arrays
  vtkNew(vtkDoubleArray, USingleKnotArray);
  vtkNew(vtkIntArray,    UMultArray);
  NURBSSurface->GetUMultiplicity(UMultArray, USingleKnotArray);
  vtkNew(vtkDoubleArray, VSingleKnotArray);
  vtkNew(vtkIntArray,    VMultArray);
  NURBSSurface->GetVMultiplicity(VMultArray, VSingleKnotArray);
  vtkSVControlGrid *controlPointGrid = NURBSSurface->GetControlPointGrid();
  int uDegree = NURBSSurface->GetUDegree();
  int vDegree = NURBSSurface->GetVDegree();

  std::cout <<"AFTER!!" << endl;
  std::cout <<"U SINGLE KNOTS"<<endl;
  vtkSVNURBSUtils::PrintArray(USingleKnotArray);
  std::cout <<"U MULTS"<<endl;
  vtkSVNURBSUtils::PrintArray(UMultArray);
  std::cout <<"V SINGLE KNOTS"<<endl;
  vtkSVNURBSUtils::PrintArray(VSingleKnotArray);
  std::cout <<"V MULTS"<<endl;
  vtkSVNURBSUtils::PrintArray(VMultArray);
  std::cout<<"CONTROL GRID" << endl;
  vtkSVNURBSUtils::PrintStructuredGrid(controlPointGrid);
  std::cout<<"U DEGREE " << uDegree << endl;
  std::cout<<"V DEGREE " << uDegree << endl;


  // Get all information needed by creation of bspline surface
  int dims[3];
  controlPointGrid->GetDimensions(dims);
  int Xlen1 = dims[0];
  int Xlen2 = dims[1];
  double **Xarr = new double*[Xlen1];
  double **Yarr = new double*[Xlen1];
  double **Zarr = new double*[Xlen1];
  for (int j=0; j<Xlen1; j++)
  {
    Xarr[j] = new double[Xlen2];
    Yarr[j] = new double[Xlen2];
    Zarr[j] = new double[Xlen2];
    for (int k=0; k<Xlen2; k++)
    {
      double pt[3];
      double w;
      controlPointGrid->GetControlPoint(j, k, 0, pt, w);
      Xarr[j][k] = pt[0];
      Yarr[j][k] = pt[1];
      Zarr[j][k] = pt[2];
    }
  }

  // Get mult information as arrays
  int uKlen = USingleKnotArray->GetNumberOfTuples();
  double *uKarr = new double[uKlen];
  for (int j=0; j<uKlen; j++)
    uKarr[j] = USingleKnotArray->GetTuple1(j);

  int vKlen = VSingleKnotArray->GetNumberOfTuples();
  double *vKarr = new double[vKlen];
  for (int j=0; j<vKlen; j++)
    vKarr[j] = VSingleKnotArray->GetTuple1(j);

  int uMlen = UMultArray->GetNumberOfTuples();
  double *uMarr = new double[uMlen];
  for (int j=0; j<uMlen; j++)
    uMarr[j] = UMultArray->GetTuple1(j);

  int vMlen = VMultArray->GetNumberOfTuples();
  double *vMarr = new double[vMlen];
  for (int j=0; j<vMlen; j++)
    vMarr[j] = VMultArray->GetTuple1(j);

  // Flipping order!
  if (solidModel->CreateBSplineCap(Xarr, Yarr, Zarr,
                                 Xlen1, Xlen2,
                                 vKarr, vKlen,
                                 uKarr, uKlen,
                                 vMarr, vMlen,
                                 uMarr, uMlen,
                                 vDegree, uDegree) != SV_OK )
  {
      std::cerr << "poly manipulation error " << endl;
      delete solidModel;
      //Clean up
      for (int j=0;j<Xlen1;j++)
      {
        delete [] Xarr[j];
        delete [] Yarr[j];
        delete [] Zarr[j];
      }
      delete [] Xarr;
      delete [] Yarr;
      delete [] Zarr;

      delete [] uKarr;
      delete [] vKarr;
      delete [] uMarr;
      delete [] vMarr;
      return SV_ERROR;
  }
  //Clean up
  for (int j=0;j<Xlen1;j++)
  {
    delete [] Xarr[j];
    delete [] Yarr[j];
    delete [] Zarr[j];
  }
  delete [] Xarr;
  delete [] Yarr;
  delete [] Zarr;

  delete [] uKarr;
  delete [] vKarr;
  delete [] uMarr;
  delete [] vMarr;

  return SV_OK;
}

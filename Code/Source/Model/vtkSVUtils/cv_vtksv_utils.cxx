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
#include "vtkSVLoftNURBSSurface.h"
#include "vtkSVPolycubeGenerator.h"
#include "vtkSVParameterizeSurfaceOnPolycube.h"
#include "vtkSVParameterizeVolumeOnPolycube.h"
#include "vtkSVPassDataArray.h"
#include "vtkSVSurfaceCenterlineAttributesPasser.h"
#include "vtkSVSurfaceCenterlineGrouper.h"
#include "vtkSVSurfaceCuboidPatcher.h"

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
    patcher->EnforceBoundaryDirectionsOn();
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
      return SV_OK;
    }

    if (vtkSVGeneralUtils::CheckArrayExists(polycubeUg, 0, "GridIds") != SV_OK)
    {
      std::cerr << "Grid point ids array with name GridIds does not exist on volume polycube" << endl;
      return SV_OK;
    }

    vtkNew(vtkPoints, savePoints);
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

      vtkDataArray *ptIds = thresholdPd->GetPointData()->GetArray("GridIds");

      int dim[3];
      dim[0] = whl_divs[1]; dim[1] = whl_divs[2]; dim[2] = whl_divs[3];

      // Set up new structured grid for nurbs surface lofting
      int newDims[3];
      newDims[0] = whl_divs[3]; newDims[1] = 2*whl_divs[1] + 2*whl_divs[2] - 3; newDims[2] = 1;
      //if (groupId == 3)
      //{
      //  continue;
      //}
      //else if (groupId == 2)
      //{
      //  newDims[0] = 4; newDims[1] = savePoints->GetNumberOfPoints(); newDims[2] = 1;
      //}
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
            //if (groupId == 0)
            //{
            //  if (j == 0)
            //  {
            //    savePoints->InsertNextPoint(pt);
            //  }
            //}
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
            //if (groupId == 0)
            //{
            //  if (j == 0)
            //  {
            //    savePoints->InsertNextPoint(pt);
            //  }
            //}
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
            //if (groupId == 0)
            //{
            //  if (j == 0)
            //  {
            //    savePoints->InsertNextPoint(pt);
            //  }
            //}
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
            //if (groupId == 0)
            //{
            //  if (j == 0)
            //  {
            //    savePoints->InsertNextPoint(pt);
            //  }
            //}
          }
        }

        //if (groupId == 2)
        //{
        //  for (int k=0; k<savePoints->GetNumberOfPoints(); k++)
        //  {
        //    savePoints->GetPoint(k, pt);

        //    newPos[0] = j; newPos[1] = k; newPos[2] = 0;
        //    ptId = vtkStructuredData::ComputePointId(newDims, newPos);

        //    for (int l=0; l<3; l++)
        //    {
        //      pt[2] = pt[2] - j*10.0;
        //    }
        //    inputGrid->GetPoints()->SetPoint(ptId, pt);
        //  }
        //}
      }

      // Now loft each surface
      int uDegree = 3;
      int vDegree = 3;
      vtkNew(vtkSVLoftNURBSSurface,lofter);
      lofter->SetInputData(inputGrid);
      lofter->SetUDegree(uDegree);
      lofter->SetVDegree(vDegree);
      lofter->SetPolyDataUSpacing(0.1);
      lofter->SetPolyDataVSpacing(0.1);
      lofter->SetUKnotSpanType("average");
      lofter->SetVKnotSpanType("average");
      lofter->SetUParametricSpanType("chord");
      lofter->SetVParametricSpanType("chord");
      lofter->Update();

      vtkNew(vtkSVNURBSSurface, NURBSSurface);
      NURBSSurface->DeepCopy(lofter->GetSurface());

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
      for (int i=0; i<Xlen1; i++)
      {
        Xarr[i] = new double[Xlen2];
        Yarr[i] = new double[Xlen2];
        Zarr[i] = new double[Xlen2];
        for (int j=0; j<Xlen2; j++)
        {
          double pt[3];
          double w;
          controlPointGrid->GetControlPoint(i, j, 0, pt, w);
          Xarr[i][j] = pt[0];
          Yarr[i][j] = pt[1];
          Zarr[i][j] = pt[2];
        }
      }

      // Get mult information as arrays
      int uKlen = USingleKnotArray->GetNumberOfTuples();
      double *uKarr = new double[uKlen];
      for (int i=0; i<uKlen; i++)
        uKarr[i] = USingleKnotArray->GetTuple1(i);

      int vKlen = VSingleKnotArray->GetNumberOfTuples();
      double *vKarr = new double[vKlen];
      for (int i=0; i<vKlen; i++)
        vKarr[i] = VSingleKnotArray->GetTuple1(i);

      int uMlen = UMultArray->GetNumberOfTuples();
      double *uMarr = new double[uMlen];
      for (int i=0; i<uMlen; i++)
        uMarr[i] = UMultArray->GetTuple1(i);

      int vMlen = VMultArray->GetNumberOfTuples();
      double *vMarr = new double[vMlen];
      for (int i=0; i<vMlen; i++)
        vMarr[i] = VMultArray->GetTuple1(i);

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
          for (int i=0;i<Xlen1;i++)
          {
            delete [] Xarr[i];
            delete [] Yarr[i];
            delete [] Zarr[i];
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
      for (int i=0;i<Xlen1;i++)
      {
        delete [] Xarr[i];
        delete [] Yarr[i];
        delete [] Zarr[i];
      }
      delete [] Xarr;
      delete [] Yarr;
      delete [] Zarr;

      delete [] uKarr;
      delete [] vKarr;
      delete [] uMarr;
      delete [] vMarr;

      //cvOCCTSolidModel* surfCapped=new cvOCCTSolidModel();
      //if ( surfCapped->CapSurfToSolid(surf) != SV_OK )
      //{
      //    std::cerr << "error in cap / bound operation ";
      //    delete surf;
      //    delete surfCapped;
      //    return NULL;
      //}
      ////delete surf

      loftedSurfs.push_back(surf);
      //loftedSurfs.push_back(surfCapped);
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
    fprintf(stderr,"ERROR in centerline separation.\n");
    fflush(stderr);
    return SV_ERROR;
  }

  return SV_OK;
}

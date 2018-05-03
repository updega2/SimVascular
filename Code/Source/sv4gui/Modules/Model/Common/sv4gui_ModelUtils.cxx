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

#include "sv4gui_ModelUtils.h"

#include "sv4gui_ModelElementPolyData.h"
#include "sv4gui_PathElement.h"
#include "sv4gui_Math3.h"
#include "sv4gui_MitkSeg3D.h"

#include "SimVascular.h"
#include "cv_sys_geom.h"
#include "cv_vmtk_utils.h"
#include "cv_vtksv_utils.h"
#include "cvPolyData.h"
#include "cv_polydatasolid_utils.h"
#include "cv_occtsolid_utils.h"

#include <vtkCellType.h>
#include <vtkFillHolesFilter.h>
#include <vtkPlaneSource.h>
#include <vtkClipPolyData.h>
#include <vtkImplicitDataSet.h>
#include <vtkThreshold.h>
#include <vtkDataSetSurfaceFilter.h>
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVNURBSSurface.h"

#include "BRepCheck.hxx"
#include "BRepCheck_ListOfStatus.hxx"
#include "BRepCheck_DataMapOfShapeListOfStatus.hxx"
#include "BRepCheck_Analyzer.hxx"
#include "BRepCheck_Face.hxx"
#include "BRepCheck_Shell.hxx"
#include "BRepCheck_ListOfStatus.hxx"
#include "BRepCheck_ListIteratorOfListOfStatus.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepBuilderAPI_MakeVertex.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_NurbsConvert.hxx"
#include "BRepGProp.hxx"
#include "BRep_Builder.hxx"
#include "BRep_Tool.hxx"
#include "BRep_TFace.hxx"
#include "BRepTools_ReShape.hxx"
#include "BRepExtrema_ExtPC.hxx"
#include "BRepExtrema_DistShapeShape.hxx"
#include "BRepLib_FuseEdges.hxx"
#include "BRepLib_FindSurface.hxx"

#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "GProp_GProps.hxx"

#include "GeomAPI_ProjectPointOnCurve.hxx"
#include "GeomConvert_CompCurveToBSplineCurve.hxx"
#include "Geom2d_Line.hxx"

#include "TopoDS.hxx"
#include "TopoDS_Edge.hxx"
#include "TopoDS_Wire.hxx"

#include "ShapeFix_FreeBounds.hxx"
#include "ShapeFix_Shape.hxx"
#include "Standard_Real.hxx"

#include "vtkXMLPolyDataWriter.h"

vtkPolyData* sv4guiModelUtils::CreatePolyData(std::vector<sv4guiContourGroup*> groups, std::vector<vtkPolyData*> vtps, int numSamplingPts, svLoftingParam *param, unsigned int t, int noInterOut, double tol)
{
    int groupNumber=groups.size();
    int vtpNumber=vtps.size();
    cvPolyData **srcs=new cvPolyData* [groupNumber+vtpNumber];
    for(int i=0;i<groupNumber;i++)
    {
      sv4guiContourGroup* group=groups[i];
      vtkPolyData *vtkpd = CreateLoftSurface(group,numSamplingPts,1,param,t);
      if (vtkpd == NULL)
      {
        for (int j=0; j< i-1; j++)
          delete srcs[j];
        delete [] srcs;
        return NULL;
      }
      srcs[i]=new cvPolyData(vtkpd);
      vtkpd->Delete();
    }

    for(int i=0;i<vtpNumber;i++)
    {
        if(vtps[i]==NULL)
        {
            for(int j=0;j<i+groupNumber-1;j++)
                delete srcs[j];
            delete [] srcs;
            return NULL;
        }
        vtkPolyData* newvtp=vtkPolyData::New();
        newvtp->DeepCopy(vtps[i]);

        if(!newvtp->GetCellData()->HasArray("CapID"))
        {
            vtkIntArray* capArray=vtkIntArray::New();
            capArray->SetName("CapID");
            for(int i=0;i<newvtp->GetNumberOfCells();i++)
                capArray->InsertNextValue(-1);

            newvtp->GetCellData()->AddArray(capArray);
        }

        cvPolyData* cvpd=new cvPolyData(newvtp);
        srcs[i+groupNumber]=cvpd;
        newvtp->Delete();
    }

    cvPolyData *dst=NULL;

    int status=sys_geom_all_union(srcs, groupNumber+vtpNumber, noInterOut, tol, &dst);

    for(int i=0;i<groupNumber+vtpNumber;i++)
    {
        delete srcs[i];
    }
    delete [] srcs;

    if(status!=SV_OK)
        return NULL;
    else
        return dst->GetVtkPolyData();
}

sv4guiModelElementPolyData* sv4guiModelUtils::CreateModelElementPolyData(std::vector<mitk::DataNode::Pointer> segNodes, int numSamplingPts, int stats[], svLoftingParam *param, unsigned int t, int noInterOut, double tol)
{
    std::vector<sv4guiContourGroup*> groups;
    std::vector<vtkPolyData*> vtps;
    std::vector<std::string> segNames;

    for(int i=0;i<segNodes.size();i++)
    {
        mitk::DataNode::Pointer segNode=segNodes[i];
        sv4guiContourGroup* group = dynamic_cast<sv4guiContourGroup*>(segNode->GetData());
        if(group!=NULL)
        {
            groups.push_back(group);
            segNames.push_back(segNode->GetName());
        }
    }

    for(int i=0;i<segNodes.size();i++)
    {
        mitk::DataNode::Pointer segNode=segNodes[i];
        sv4guiMitkSeg3D* seg3D=dynamic_cast<sv4guiMitkSeg3D*>(segNode->GetData());
        if(seg3D && seg3D->GetVtkPolyData())
        {
            vtps.push_back(seg3D->GetVtkPolyData());
            segNames.push_back(segNode->GetName());
        }
    }

    vtkPolyData* solidvpd=CreatePolyData(groups,vtps,numSamplingPts,param,t,noInterOut,tol);
    if(solidvpd==NULL) return NULL;

    cvPolyData *src=new cvPolyData(solidvpd);
    cvPolyData *dst = NULL;

    if(stats&&sys_geom_checksurface(src,stats,tol)!=SV_OK)
    {
      solidvpd->Delete();
      return NULL;
    }

    int *doublecaps;
    int numfaces=0;

    if (sys_geom_set_ids_for_caps(src, &dst,  &doublecaps,&numfaces) != SV_OK)
    {
      solidvpd->Delete();
      return NULL;
    }

    vtkSmartPointer<vtkPolyData> forClean =
      vtkSmartPointer<vtkPolyData>::New();
    forClean->DeepCopy(dst->GetVtkPolyData());
    vtkSmartPointer<vtkPolyData> nowClean =
      vtkSmartPointer<vtkPolyData>::New();
    nowClean = sv4guiModelUtils::OrientVtkPolyData(forClean);

    solidvpd->DeepCopy(nowClean);;

    int numSeg=segNames.size();
    int numCap2=0;
    for(int i=numSeg-1;i>-1;i--)
    {
        if(doublecaps[i]!=0)
        {
            numCap2=doublecaps[i];
            break;
        }
    }

    std::string *allNames=new std::string[2*numSeg+numCap2];

    for(int i=0;i<numSeg;i++)
    {
        allNames[i]="wall_"+segNames[i];
        allNames[numSeg+i]="cap_"+segNames[i];
        if(doublecaps[i]!=0)
            allNames[2*numSeg+doublecaps[i]-1]="cap_"+segNames[i]+"_2";
    }

    std::vector<sv4guiModelElement::svFace*> faces;

    for(int i=0;i<2*numSeg+numCap2;i++)
    {
        vtkPolyData *facepd = vtkPolyData::New();
        int faceid=i+1;
        PlyDtaUtils_GetFacePolyData(solidvpd, &faceid, facepd);

        if(facepd==NULL||facepd->GetNumberOfPoints()==0)
            continue;

        sv4guiModelElement::svFace* face =new sv4guiModelElement::svFace;
        face->id=faceid;
        face->name=allNames[i];
        face->vpd=facepd;

        if(face->name.substr(0,5)=="wall_")
            face->type="wall";
        else if(face->name.substr(0,4)=="cap_")
            face->type="cap";

        faces.push_back(face);
    }

    delete[] allNames;

    sv4guiModelElementPolyData* modelElement=new sv4guiModelElementPolyData();
    modelElement->SetSegNames(segNames);
    modelElement->SetFaces(faces);
    modelElement->SetWholeVtkPolyData(solidvpd);
    modelElement->SetNumSampling(numSamplingPts);


    bool ok = false;
    if(modelElement->MarkCellsByFaces(modelElement->GetCapFaceIDs()))
    {
      int numDivs = 1;
      ok=modelElement->LinearSubdivideLocal(numDivs);
    }
    if(!ok)
    {
      MITK_ERROR << "Failed to subdivide caps of created PolyData";
      return NULL;
    }

    return modelElement;
}

vtkPolyData* sv4guiModelUtils::CreatePolyDataByBlend(vtkPolyData* vpdsrc, int faceID1, int faceID2, double radius, sv4guiModelElement::svBlendParam* param)
{
    if(vpdsrc==NULL)
        return NULL;

    cvPolyData *src=new cvPolyData(vpdsrc);
    cvPolyData *dst = NULL;

    int vals[2];
    vals[0]=faceID1;
    vals[1]=faceID2;

    if ( sys_geom_set_array_for_local_op_face_blend(src,&dst, "ModelFaceID", vals, 2, radius, "ActiveCells", 1)
         != SV_OK )
    {
        MITK_ERROR << "poly blend (using radius) error ";
        return NULL;
    }

    cvPolyData *dst2 = NULL;

    if ( sys_geom_local_blend( dst, &dst2, param->numblenditers,
                               param->numsubblenditers, param->numsubdivisioniters,
                               param->numcgsmoothiters, param->numlapsmoothiters,
                               param->targetdecimation,
                               NULL, "ActiveCells")

         != SV_OK )
    {
        MITK_ERROR << "poly blend error ";
        return NULL;
    }

    vtkPolyData* vpd=dst2->GetVtkPolyData();
    vpd->GetCellData()->RemoveArray("ActiveCells");

    return vpd;

}

sv4guiModelElementPolyData* sv4guiModelUtils::CreateModelElementPolyDataByBlend(sv4guiModelElementPolyData* mepdsrc, std::vector<sv4guiModelElement::svBlendParamRadius*> blendRadii, sv4guiModelElement::svBlendParam* param)
{

    vtkSmartPointer<vtkPolyData> oldVpd=mepdsrc->GetWholeVtkPolyData();
    if(oldVpd==NULL) return NULL;

    vtkSmartPointer<vtkPolyData> lastVpd=oldVpd;

    for(int i=0;i<blendRadii.size();i++)
    {
        int faceID1=0;
        int faceID2=0;
        double radius=0.0;
        if(blendRadii[i] && blendRadii[i]->radius>0)
        {
            faceID1=blendRadii[i]->faceID1;
            faceID2=blendRadii[i]->faceID2;
            radius=blendRadii[i]->radius;

        }

        lastVpd=sv4guiModelUtils::CreatePolyDataByBlend(lastVpd, faceID1, faceID2, radius, param);

        if(lastVpd==NULL) return NULL;

    }

    sv4guiModelElementPolyData* mepddst =mepdsrc->Clone();
    mepddst->SetWholeVtkPolyData(lastVpd);
    std::vector<sv4guiModelElement::svFace*> faces=mepddst->GetFaces();
    for(int i=0;i<faces.size();i++)
    {
        faces[i]->vpd=mepddst->CreateFaceVtkPolyData(faces[i]->id);
    }

    mepddst->AssignBlendParam(param);

    mepddst->AddBlendRadii(blendRadii);

    return mepddst;
}

vtkPolyData* sv4guiModelUtils::CreateLoftSurface(sv4guiContourGroup* contourGroup, int numSamplingPts, int addCaps, svLoftingParam* param, unsigned int t)
{

    svLoftingParam* usedParam= contourGroup->GetLoftingParam();
    if(param!=NULL) usedParam=param;

    std::vector<sv4guiContour*> contourSet=contourGroup->GetValidContourSet(t);

    return CreateLoftSurface(contourSet,numSamplingPts,usedParam,addCaps);
}

vtkPolyData* sv4guiModelUtils::CreateLoftSurface(std::vector<sv4guiContour*> contourSet, int numSamplingPts, svLoftingParam* param, int addCaps)
{
    int contourNumber=contourSet.size();
    if (contourNumber < 2)
      return NULL;

    if(param==NULL)
        return NULL;

    param->numOutPtsAlongLength=param->samplePerSegment*contourNumber;
    param->numPtsInLinearSampleAlongLength=param->linearMuliplier*param->numOutPtsAlongLength;

    param->numSuperPts=0;
    for(int i=0;i<contourNumber;i++)
    {
        int pointNunumber=contourSet[i]->GetContourPointNumber();
        if(pointNunumber>param->numSuperPts)
            param->numSuperPts=pointNunumber;
    }

    if(param->numOutPtsInSegs>param->numSuperPts)
        param->numSuperPts=param->numOutPtsInSegs;

    int newNumSamplingPts=param->numOutPtsInSegs;

    if(numSamplingPts>3)
    {
        newNumSamplingPts=numSamplingPts;
        if(numSamplingPts>param->numSuperPts)
            param->numSuperPts=numSamplingPts;
    }

    std::vector<cvPolyData*> superSampledContours;
    for(int i=0;i<contourNumber;i++)
    {
        vtkPolyData* vtkpd=vtkPolyData::New();

        vtkpd->DeepCopy(contourSet[i]->CreateVtkPolyDataFromContour(false));
        cvPolyData* cvpd=new cvPolyData(vtkpd);
        vtkpd->Delete();
        cvPolyData* cvpd2=sys_geom_sampleLoop(cvpd,param->numSuperPts);
        if(cvpd2==NULL)
        {
            MITK_ERROR << "Supersampling error ";
            return NULL;
        }
        superSampledContours.push_back(cvpd2);
    }

    std::vector<cvPolyData*> alignedContours;
    for(int i=0;i<contourNumber;i++)
    {
        if(i==0)
        {
            alignedContours.push_back(superSampledContours[0]);
        }
        else
        {
            cvPolyData* cvpd3;
            if(param->vecFlag==1)
                cvpd3=sys_geom_Align(alignedContours[i-1],superSampledContours[i]);
            else
                cvpd3=sys_geom_AlignByDist(alignedContours[i-1],superSampledContours[i]);

            if(cvpd3==NULL)
            {
                MITK_ERROR << "aligning error ";
                // Clean up
                for (int i=0; i<contourNumber; i++)
                  delete superSampledContours[i];

                return NULL;
            }

            alignedContours.push_back(cvpd3);
        }
    }

    cvPolyData **sampledContours=new cvPolyData*[contourNumber];
    for(int i=0;i<contourNumber;i++)
    {
        cvPolyData * cvpd4=sys_geom_sampleLoop(alignedContours[i],newNumSamplingPts);
        if(cvpd4==NULL)
        {
            MITK_ERROR << "sampling error ";
            for (int j=0; j<i; j++)
              delete sampledContours[j];
            delete [] sampledContours;
            return NULL;
        }
        sampledContours[i]=cvpd4;
    }

    cvPolyData *dst;
    vtkPolyData* outpd=NULL;

    if (param->method=="spline")
    {
      if ( sys_geom_loft_solid(sampledContours, contourNumber,param->useLinearSampleAlongLength,param->useFFT,
                               param->numOutPtsAlongLength,newNumSamplingPts,
                               param->numPtsInLinearSampleAlongLength,param->numModes,param->splineType,param->bias,param->tension,param->continuity,
                               &dst )
           != SV_OK )
      {
          MITK_ERROR << "poly manipulation error ";
          outpd=NULL;
      }
      else
      {

          if(addCaps==1)
              outpd=CreateOrientClosedPolySolidVessel(dst->GetVtkPolyData());
          else
              outpd=CreateOrientOpenPolySolidVessel(dst->GetVtkPolyData());
      }
    }
    else if (param->method=="nurbs")
    {
      // Degrees of surface
      int uDegree = param->uDegree;
      int vDegree = param->vDegree;

      // Override to maximum possible degree if too large a degree for given number of inputs!
      if (uDegree >= contourNumber)
        uDegree = contourNumber-1;
      if (vDegree >= newNumSamplingPts)
        vDegree = newNumSamplingPts-1;

      // Set to average knot span and chord length if just two inputs
      if (contourNumber == 2)
      {
        param->uKnotSpanType = "average";
        param->uParametricSpanType = "chord";
      }
      if (newNumSamplingPts == 2)
      {
        param->uKnotSpanType = "average";
        param->uParametricSpanType = "chord";
      }

      // Output spacing function of given input points
      double uSpacing = 1.0/param->numOutPtsAlongLength;
      double vSpacing = 1.0/newNumSamplingPts;

      // span types
      const char *uKnotSpanType       = param->uKnotSpanType.c_str();
      const char *vKnotSpanType       = param->vKnotSpanType.c_str();
      const char *uParametricSpanType = param->uParametricSpanType.c_str();
      const char *vParametricSpanType = param->vParametricSpanType.c_str();
      vtkNew(vtkSVNURBSSurface, NURBSSurface);

      if ( sys_geom_loft_solid_with_nurbs(sampledContours, contourNumber,
                                          uDegree, vDegree, uSpacing,
                                          vSpacing, uKnotSpanType,
                                          vKnotSpanType,
                                          uParametricSpanType,
                                          vParametricSpanType,
                                          NURBSSurface,
                                          &dst )
           != SV_OK )
      {
          MITK_ERROR << "poly manipulation error ";
          outpd=NULL;
      }
      else
      {
          if (PlyDtaUtils_CheckLoftSurface(dst->GetVtkPolyData()) != SV_OK)
          {
            MITK_ERROR << "Error lofting surface";
            outpd=NULL;
          }
          else
          {
            if(addCaps==1)
                outpd=CreateOrientClosedPolySolidVessel(dst->GetVtkPolyData());
            else
                outpd=CreateOrientOpenPolySolidVessel(dst->GetVtkPolyData());
          }

      }
    }

    // Clean up
    for (int i=0; i<contourNumber; i++)
    {
      delete superSampledContours[i];
      delete sampledContours[i];
    }
    delete [] sampledContours;

    if(dst!=NULL) delete dst;

    return outpd;

}

vtkPolyData* sv4guiModelUtils::CreateOrientOpenPolySolidVessel(vtkPolyData* inpd)
{
    int originalCellNumber=inpd->GetNumberOfCells();

    vtkPolyData* tmppd=FillHoles(inpd);

    vtkSmartPointer<vtkPolyDataNormals> nrmls = vtkSmartPointer<vtkPolyDataNormals>::New();
    nrmls->SplittingOff();
    nrmls->ConsistencyOn();
    nrmls->AutoOrientNormalsOn();
    nrmls->ComputeCellNormalsOn();
    nrmls->ComputePointNormalsOff();
    nrmls->SetInputData(tmppd);
    nrmls->Update();

    vtkPolyData* outpd=vtkPolyData::New();
    outpd->DeepCopy(nrmls->GetOutput());
    //    int VTK_TRIANGLE=5;
    for(int i=outpd->GetNumberOfCells()-1;i>=originalCellNumber;i--)
    {
        if(outpd->GetCellType(i)==VTK_TRIANGLE)
            outpd->DeleteCell(i);
    }
    outpd->RemoveDeletedCells();

    tmppd->Delete();

    return outpd;
}

vtkPolyData* sv4guiModelUtils::FillHoles(vtkPolyData* inpd)
{
    vtkSmartPointer<vtkFillHolesFilter> filler = vtkSmartPointer<vtkFillHolesFilter>::New();
    filler->SetHoleSize(filler->GetHoleSizeMaxValue());
    filler->SetInputData(inpd);
    filler->Update();

    return Orient(filler->GetOutput());
}

vtkPolyData* sv4guiModelUtils::Orient(vtkPolyData* inpd)
{
    vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->PointMergingOn();
    cleaner->ConvertLinesToPointsOff();
    cleaner->ConvertPolysToLinesOff();
    cleaner->SetInputData(inpd);
    cleaner->Update();

    vtkSmartPointer<vtkPolyDataNormals> orienter = vtkSmartPointer<vtkPolyDataNormals>::New();
    orienter->SetInputData(cleaner->GetOutput());
    orienter->AutoOrientNormalsOn();
    orienter->ComputePointNormalsOff();
    orienter->FlipNormalsOn();
    orienter->SplittingOff();
    orienter->ComputeCellNormalsOn();
    orienter->ConsistencyOn();
    orienter->NonManifoldTraversalOff();
    orienter->Update();

    vtkPolyData* outpd=vtkPolyData::New();
    outpd->DeepCopy(orienter->GetOutput());

    return outpd;
}

vtkPolyData* sv4guiModelUtils::CreateOrientClosedPolySolidVessel(vtkPolyData* inpd)
{
    int fillID=0;
    int fillType=0;

    vtkPolyData* tmppd=FillHolesWithIDs(inpd,fillID,fillType);

    vtkSmartPointer<vtkPolyDataNormals> nrmls = vtkSmartPointer<vtkPolyDataNormals>::New();
    nrmls->SplittingOff();
    nrmls->ConsistencyOn();
    nrmls->AutoOrientNormalsOn();
    nrmls->ComputeCellNormalsOn();
    nrmls->ComputePointNormalsOff();
    nrmls->SetInputData(tmppd);
    nrmls->Update();

    vtkPolyData* outpd=vtkPolyData::New();
    outpd->DeepCopy(nrmls->GetOutput());

    tmppd->Delete();

    return outpd;
}

vtkPolyData* sv4guiModelUtils::FillHolesWithIDs(vtkPolyData* inpd, int fillID, int fillType)
{
    cvPolyData* cvpd=new cvPolyData(inpd);
    int numFilled=0;
    cvPolyData* tmpcvpd;
    if(VMTKUtils_CapWithIds(cvpd,&tmpcvpd,fillID,numFilled,fillType)!=SV_OK)
        return NULL;

    if(tmpcvpd==NULL)
        return NULL;

    vtkPolyData* outpd=Orient(tmpcvpd->GetVtkPolyData());

    delete tmpcvpd;

    return outpd;
}

bool sv4guiModelUtils::CheckArrayName(vtkDataSet *object,int datatype,std::string arrayname )
{
    vtkIdType i;
    int numArrays;

    if (datatype == 0)
    {
        numArrays = object->GetPointData()->GetNumberOfArrays();
        for (i=0;i<numArrays;i++)
        {
            if (strcmp(object->GetPointData()->GetArrayName(i),arrayname.c_str())==0)
            {
                return true;
            }
        }

        //    if(object->GetPointData()->HasArray(arrayname.c_str()))
        //        return true;

    }
    else
    {
        numArrays = object->GetCellData()->GetNumberOfArrays();
        for (i=0;i<numArrays;i++)
        {
            if (strcmp(object->GetCellData()->GetArrayName(i),arrayname.c_str())==0)
            {
                return true;
            }
        }

        //    if(object->GetCellData()->HasArray(arrayname.c_str()))
        //        return true;
    }

    return false;
}

vtkSmartPointer<vtkPolyData> sv4guiModelUtils::OrientVtkPolyData(vtkSmartPointer<vtkPolyData> inpd)
{
    vtkSmartPointer<vtkCleanPolyData> cleaner=vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->PointMergingOn();
    cleaner->ConvertLinesToPointsOff();
    cleaner->ConvertPolysToLinesOff();
    cleaner->SetInputDataObject(inpd);
    cleaner->Update();

    vtkSmartPointer<vtkPolyDataNormals> orienter=vtkSmartPointer<vtkPolyDataNormals>::New();
    orienter->SetInputDataObject(cleaner->GetOutput());
    orienter->AutoOrientNormalsOn();
    orienter->ComputePointNormalsOff();
    orienter->SplittingOff();
    orienter->ComputeCellNormalsOn();
    orienter->ConsistencyOn();
    orienter->NonManifoldTraversalOff();
    orienter->Update();

    return orienter->GetOutput();
}

vtkSmartPointer<vtkPolyData> sv4guiModelUtils::MarkCells(vtkSmartPointer<vtkPolyData> inpd, std::vector<int> cellIDs)
{
    if(inpd==NULL)
        return NULL;

    cvPolyData *src=new cvPolyData(inpd);
    cvPolyData *dst = NULL;

    int* cellIDArray = &cellIDs[0];

    if ( sys_geom_set_array_for_local_op_cells(src, &dst, cellIDArray, cellIDs.size(), "ActiveCells", 1) != SV_OK )
    {
        MITK_ERROR << "poly marking cells (by cell ids) error ";
        return NULL;
    }

    return dst->GetVtkPolyData();
}


vtkSmartPointer<vtkPolyData> sv4guiModelUtils::MarkCellsBySphere(vtkSmartPointer<vtkPolyData> inpd, double radius, double center[3])
{
    if(inpd==NULL)
        return NULL;

    cvPolyData *src=new cvPolyData(inpd);
    cvPolyData *dst = NULL;

    if ( sys_geom_set_array_for_local_op_sphere(src, &dst, radius, center, "ActiveCells", 1) != SV_OK )
    {
        MITK_ERROR << "poly marking cells (by sphere) error ";
        return NULL;
    }

    return dst->GetVtkPolyData();
}

vtkSmartPointer<vtkPolyData> sv4guiModelUtils::MarkCellsByFaces(vtkSmartPointer<vtkPolyData> inpd, std::vector<int> faceIDs)
{
    if(inpd==NULL)
        return NULL;

    cvPolyData *src=new cvPolyData(inpd);
    cvPolyData *dst = NULL;

    int* faceIDArray = &faceIDs[0];

    if ( sys_geom_set_array_for_local_op_face(src, &dst, "ModelFaceID", faceIDArray, faceIDs.size(), "ActiveCells", 1) != SV_OK )
    {
        MITK_ERROR << "poly marking cells (by face ids) error ";
        return NULL;
    }

    return dst->GetVtkPolyData();
}

vtkSmartPointer<vtkPolyData> sv4guiModelUtils::MarkCellsByFaceJunctions(vtkSmartPointer<vtkPolyData> inpd, std::vector<int> faceIDs, double radius)
{
    if(inpd==NULL)
        return NULL;

    cvPolyData *src=new cvPolyData(inpd);
    cvPolyData *dst = NULL;

    int* faceIDArray = &faceIDs[0];

    if ( sys_geom_set_array_for_local_op_face_blend(src, &dst, "ModelFaceID", faceIDArray, faceIDs.size(), radius, "ActiveCells", 1) != SV_OK )
    {
        MITK_ERROR << "poly marking cells (by face functions) error ";
        return NULL;
    }

    return dst->GetVtkPolyData();
}

vtkSmartPointer<vtkPolyData> sv4guiModelUtils::DecimateLocal(vtkSmartPointer<vtkPolyData> inpd, double targetRate)
{
    if(inpd==NULL)
        return NULL;

    cvPolyData *src=new cvPolyData(inpd);
    cvPolyData *dst = NULL;

    if ( sys_geom_local_quadric_decimation(src, &dst, targetRate, NULL, "ActiveCells") != SV_OK )
    {
        MITK_ERROR << "poly local decimation error ";
        return NULL;
    }

    return dst->GetVtkPolyData();
}

vtkSmartPointer<vtkPolyData> sv4guiModelUtils::LaplacianSmoothLocal(vtkSmartPointer<vtkPolyData> inpd, int numIters, double relaxFactor)
{
    if(inpd==NULL)
        return NULL;

    cvPolyData *src=new cvPolyData(inpd);
    cvPolyData *dst = NULL;

    if ( sys_geom_local_laplacian_smooth(src, &dst, numIters, relaxFactor, NULL, "ActiveCells") != SV_OK )
    {
        MITK_ERROR << "poly local Laplacian smooth error ";
        return NULL;
    }

    return dst->GetVtkPolyData();
}

vtkSmartPointer<vtkPolyData> sv4guiModelUtils::ConstrainSmoothLocal(vtkSmartPointer<vtkPolyData> inpd, int numIters, double constrainFactor, int numCGSolves)
{
    if(inpd==NULL)
        return NULL;

    cvPolyData *src=new cvPolyData(inpd);
    cvPolyData *dst = NULL;

    if ( sys_geom_local_constrain_smooth(src, &dst, numIters, constrainFactor, numCGSolves, NULL, "ActiveCells") != SV_OK )
    {
        MITK_ERROR << "poly local constrain smooth error ";
        return NULL;
    }

    return dst->GetVtkPolyData();
}

vtkSmartPointer<vtkPolyData> sv4guiModelUtils::LinearSubdivideLocal(vtkSmartPointer<vtkPolyData> inpd, int numDivs)
{
    if(inpd==NULL)
        return NULL;

    cvPolyData *src=new cvPolyData(inpd);
    cvPolyData *dst = NULL;

    if ( sys_geom_local_linear_subdivision(src, &dst, numDivs, NULL, "ActiveCells") != SV_OK )
    {
        MITK_ERROR << "poly local linear subdivision error ";
        return NULL;
    }

    return dst->GetVtkPolyData();
}

vtkSmartPointer<vtkPolyData> sv4guiModelUtils::LoopSubdivideLocal(vtkSmartPointer<vtkPolyData> inpd, int numDivs)
{
    if(inpd==NULL)
        return NULL;

    cvPolyData *src=new cvPolyData(inpd);
    cvPolyData *dst = NULL;

    if ( sys_geom_local_loop_subdivision(src, &dst, numDivs, NULL, "ActiveCells") != SV_OK )
    {
        MITK_ERROR << "poly local loop subdivision error ";
        return NULL;
    }

    return dst->GetVtkPolyData();
}

vtkSmartPointer<vtkPolyData> sv4guiModelUtils::CutByPlane(vtkSmartPointer<vtkPolyData> inpd, double origin[3], double point1[3], double point2[3], bool above )
{
    if(inpd==NULL)
        return NULL;

    vtkSmartPointer<vtkPlaneSource> plane= vtkSmartPointer<vtkPlaneSource>::New();
    plane->SetOrigin(origin);
    plane->SetPoint1(point1);
    plane->SetPoint2(point2);
    plane->Update();

    double* nrm=plane->GetNormal();
    nrm[0]=-nrm[0];
    nrm[1]=-nrm[1];
    nrm[2]=-nrm[2];

    vtkSmartPointer<vtkPlane> impPlane=vtkSmartPointer<vtkPlane>::New();
    impPlane->SetOrigin(origin);
    impPlane->SetNormal(nrm);

    vtkSmartPointer<vtkClipPolyData> clipper=vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetInputData(inpd);
    clipper->GenerateClippedOutputOn();
    clipper->SetClipFunction(impPlane);
    clipper->Update();

    vtkSmartPointer<vtkFillHolesFilter> triangulator=vtkSmartPointer<vtkFillHolesFilter>::New();
    if(above)
    {
        triangulator->SetInputData(clipper->GetOutput());
    }
    else
    {
        triangulator->SetInputData(clipper->GetClippedOutput());
    }
    triangulator->Update();

    return triangulator->GetOutput();
}

vtkSmartPointer<vtkPolyData> sv4guiModelUtils::CutByBox(vtkSmartPointer<vtkPolyData> inpd, vtkSmartPointer<vtkPlanes> boxPlanes, bool inside)
{
    if(inpd==NULL)
        return NULL;

    if(boxPlanes==NULL)
        return NULL;

    vtkSmartPointer<vtkClipPolyData> clipper=vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetInputData(inpd);
    clipper->GenerateClippedOutputOn();
    clipper->SetClipFunction(boxPlanes);
    clipper->Update();

    vtkSmartPointer<vtkTriangleFilter> triangulator=vtkSmartPointer<vtkTriangleFilter>::New();
    if(inside)
    {
        triangulator->SetInputData(clipper->GetOutput());
    }
    else
    {
        triangulator->SetInputData(clipper->GetClippedOutput());
    }
    triangulator->Update();

    return triangulator->GetOutput();
}

bool sv4guiModelUtils::DeleteRegions(vtkSmartPointer<vtkPolyData> inpd, std::vector<int> regionIDs)
{
    if(inpd==NULL)
        return false;

    std::string arrayname="ModelFaceID";
    bool existing=false;

    if(inpd->GetCellData()->HasArray(arrayname.c_str()))
        existing=true;

    if(!existing)
        return false;

    for(int i=0;i<regionIDs.size();i++)
    {
        vtkSmartPointer<vtkIntArray> boundaryRegions = vtkSmartPointer<vtkIntArray>::New();
        boundaryRegions = vtkIntArray::SafeDownCast(inpd->GetCellData()-> GetScalars("ModelFaceID"));

        inpd->BuildLinks();

        for (vtkIdType cellId=0; cellId< inpd->GetNumberOfCells(); cellId++)
        {
            if (boundaryRegions->GetValue(cellId) == regionIDs[i])
            {
                inpd->DeleteCell(cellId);
            }
        }

        inpd->RemoveDeletedCells();
    }

    return true;
}

vtkPolyData* sv4guiModelUtils::CreateCenterlines(sv4guiModelElement* modelElement,
                                             vtkIdList *sourceCapIds, int useVmtk)
{
    if(modelElement==NULL || modelElement->GetWholeVtkPolyData()==NULL)
        return NULL;

    vtkSmartPointer<vtkPolyData> inpd=vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> fullpd=vtkSmartPointer<vtkPolyData>::New();
    inpd->DeepCopy(modelElement->GetWholeVtkPolyData());
    fullpd->DeepCopy(modelElement->GetWholeVtkPolyData());

    if(!DeleteRegions(inpd,modelElement->GetCapFaceIDs()))
    {
        return NULL;
    }

    cvPolyData *src=new cvPolyData(inpd);
    cvPolyData *tmpCleaned = NULL;
    cvPolyData *cleaned = NULL;
    cvPolyData *capped  = NULL;
    int numCapCenterIds;
    int *capCenterIds=NULL;

    tmpCleaned = sys_geom_Clean(src);
    delete src;

    // Need to triangulate in between to make sure valid input to capper
    vtkNew(vtkTriangleFilter, triangulator);
    triangulator->SetInputData(tmpCleaned->GetVtkPolyData());
    triangulator->PassLinesOff();
    triangulator->PassVertsOff();
    triangulator->Update();

    cleaned = new cvPolyData(triangulator->GetOutput());

    if ( VMTKUtils_Cap(cleaned, &capped, &numCapCenterIds, &capCenterIds, 1 ) != SV_OK)
    {
      delete cleaned;
      if (capped != NULL)
        delete capped;
      return NULL;
    }
    if (numCapCenterIds < 2)
    {
      delete cleaned;
      if (capped != NULL)
        delete capped;
      return NULL;
    }
    delete cleaned;

    vtkNew(vtkCleanPolyData, cleaner);
    cleaner->SetInputData(capped->GetVtkPolyData());
    cleaner->Update();

    vtkNew(vtkPolyData, tmpPd);
    tmpPd->DeepCopy(cleaner->GetOutput());

    for (int i=0; i<tmpPd->GetNumberOfCells(); i++)
    {
      if (tmpPd->GetCellType(i) != VTK_TRIANGLE)
      {
        tmpPd->DeleteCell(i);
      }
    }

    // Just in case things are removed by cleaner or deleted cells, reset
    // cap center ids
    tmpPd->RemoveDeletedCells();
    tmpPd->BuildLinks();

    vtkNew(vtkPointLocator, pointLocator);
    pointLocator->SetDataSet(tmpPd);
    pointLocator->BuildLocator();

    int currCapCenterId, newCapCenterId;
    double pt[3];
    for (int i=0; i<numCapCenterIds; i++)
    {
      currCapCenterId = capCenterIds[i];
      capped->GetVtkPolyData()->GetPoint(currCapCenterId, pt);

      newCapCenterId = pointLocator->FindClosestPoint(pt);
      capCenterIds[i] = newCapCenterId;
    }

    // Reset capped now
    delete capped;
    capped = new cvPolyData(tmpPd);

    vtkSmartPointer<vtkIdList> sourcePtIds = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> targetPtIds = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> capCenterPtIds = vtkSmartPointer<vtkIdList>::New();

    int capIdsGiven = 0;
    if (sourceCapIds != NULL)
    {
      if (sourceCapIds->GetNumberOfIds() > 0)
        capIdsGiven = 1;
    }

    if (!capIdsGiven)
    {
      //sourcePtIds->InsertNextId(capCenterIds[0]);
      sourcePtIds->InsertNextId(0);
      for (int i=1; i<numCapCenterIds; i++)
      {
        //targetPtIds->InsertNextId(capCenterIds[i]);
        targetPtIds->InsertNextId(i);
      }
    }
    else
    {
      vtkSmartPointer<vtkCellLocator> locator =
        vtkSmartPointer<vtkCellLocator>::New();
      locator->SetDataSet(fullpd);
      locator->BuildLocator();

      int subId;
      double distance;
      double capPt[3];
      double closestPt[3];
      vtkIdType closestCell;
      vtkSmartPointer<vtkGenericCell> genericCell =
        vtkSmartPointer<vtkGenericCell>::New();

      for (int i=0; i<numCapCenterIds; i++)
      {
        int ptId = capCenterIds[i];
        capped->GetVtkPolyData()->GetPoint(ptId, capPt);

        locator->FindClosestPoint(capPt,closestPt,genericCell,closestCell,
          subId,distance);

        int capFaceId = fullpd->GetCellData()->GetArray("ModelFaceID")->GetTuple1(closestCell);

        if (sourceCapIds->IsId(capFaceId) != -1)
        {
          //sourcePtIds->InsertNextId(ptId);
          sourcePtIds->InsertNextId(i);
        }
        else
        {
          //targetPtIds->InsertNextId(ptId);
          targetPtIds->InsertNextId(i);
        }
      }
    }

    for (int i=0; i<numCapCenterIds; i++)
      capCenterPtIds->InsertNextId(capCenterIds[i]);
    delete [] capCenterIds;

    vtkNew(vtkCleanPolyData, testCleaner);
    testCleaner->SetInputData(capped->GetVtkPolyData());
    testCleaner->Update();

    if (testCleaner->GetOutput()->GetNumberOfPoints() != capped->GetVtkPolyData()->GetNumberOfPoints() || testCleaner->GetOutput()->GetNumberOfCells() != capped->GetVtkPolyData()->GetNumberOfCells())
    {
      MITK_ERROR << "Error after cap" << endl;
      return NULL;
    }

    vtkPolyData* centerlines=CreateCenterlines(capped->GetVtkPolyData(),
                                               sourcePtIds, targetPtIds,
                                               capCenterPtIds,
                                               useVmtk);
    delete capped;

    return centerlines;
}

vtkPolyData* sv4guiModelUtils::CreateCenterlines(vtkPolyData* inpd, int useVmtk)
{
  // If given just a polydata, assume it is a wall, cap and get source and
  // target points and then send to centerline extraction

  // Cap the solid to get centerline ids
  cvPolyData *src = new cvPolyData(inpd);
  cvPolyData *cleaned = NULL;
  cvPolyData *capped  = NULL;
  int numCapCenterIds;
  int *capCenterIds=NULL;

  cleaned = sys_geom_Clean(src);

  if ( VMTKUtils_Cap(cleaned, &capped, &numCapCenterIds, &capCenterIds, 1 ) != SV_OK)
  {
    delete cleaned;
    if (capped != NULL)
      delete capped;
    return NULL;
  }
  if (numCapCenterIds < 2)
  {
    delete cleaned;
    if (capped != NULL)
      delete capped;
    return NULL;
  }
  delete cleaned;

  vtkSmartPointer<vtkIdList> sourcePtIds = vtkSmartPointer<vtkIdList>::New();
  //sourcePtIds->InsertNextId(capCenterIds[0]);
  sourcePtIds->InsertNextId(0);
  vtkSmartPointer<vtkIdList> targetPtIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> capCenterPtIds = vtkSmartPointer<vtkIdList>::New();
  for (int i=1; i<numCapCenterIds; i++)
  {
    //targetPtIds->InsertNextId(capCenterIds[i]);
    targetPtIds->InsertNextId(i);
  }

  // capped and got ids

  for (int i=0; i<numCapCenterIds; i++)
    capCenterPtIds->InsertNextId(capCenterIds[i]);
  delete [] capCenterIds;

  return CreateCenterlines(capped->GetVtkPolyData(), sourcePtIds, targetPtIds, capCenterPtIds, useVmtk);

}


vtkPolyData* sv4guiModelUtils::CreateCenterlines(vtkPolyData* inpd,
                                             vtkIdList *sourcePtIds,
                                             vtkIdList *targetPtIds,
                                             vtkIdList *capCenterPtIds,
                                             int useVmtk)
{
    if(inpd==NULL)
        return NULL;

    cvPolyData *src = new cvPolyData(inpd);
    cvPolyData *tempCenterlines = NULL;
    cvPolyData *voronoi = NULL;

    int numSourcePts = sourcePtIds->GetNumberOfIds();
    int *sources=new int[numSourcePts];
    for (int i=0; i<numSourcePts; i++)
      sources[i]=sourcePtIds->GetId(i);

    int numTargetPts = targetPtIds->GetNumberOfIds();
    int *targets=new int[numTargetPts];
    for (int i=0; i<numTargetPts; i++)
      targets[i]=targetPtIds->GetId(i);

    int numCapCenterPts = capCenterPtIds->GetNumberOfIds();
    int *capCenters = new int[numCapCenterPts];
    for (int i=0; i<numCapCenterPts; i++)
      capCenters[i]=capCenterPtIds->GetId(i);

    if ( VMTKUtils_Centerlines(src, sources, numSourcePts, targets, numTargetPts, capCenters, numCapCenterPts, useVmtk, &tempCenterlines, &voronoi) != SV_OK )
    {
        delete src;
        delete [] sources;
        delete [] targets;
        delete [] capCenters;
        return NULL;
    }
    delete src;
    delete voronoi;
    delete [] sources;
    delete [] targets;
    delete [] capCenters;

    cvPolyData *centerlines=NULL;
    if ( VMTKUtils_SeparateCenterlines(tempCenterlines, useVmtk, &centerlines) != SV_OK )
    {
        delete tempCenterlines;
        return NULL;
    }
    delete tempCenterlines;

    return centerlines->GetVtkPolyData();
}

vtkPolyData* sv4guiModelUtils::MergeCenterlines(vtkPolyData* centerlinesPD, int useVmtk)
{
    if(centerlinesPD==NULL)
        return NULL;

    cvPolyData *centerlines =new cvPolyData(centerlinesPD);

    cvPolyData *merged_centerlines=NULL;
    int mergeblanked = 1;
    if (VMTKUtils_MergeCenterlines(centerlines, mergeblanked, useVmtk, &merged_centerlines) != SV_OK )
    {
      delete centerlines;
      return NULL;
    }
    delete centerlines;

    return merged_centerlines->GetVtkPolyData();
}

vtkPolyData* sv4guiModelUtils::CalculateDistanceToCenterlines(vtkPolyData* centerlines, vtkPolyData* original)
{
    if(centerlines==NULL || original==NULL)
        return NULL;

    cvPolyData *src=new cvPolyData(original);
    cvPolyData *lines=new cvPolyData(centerlines);
    cvPolyData *distance = NULL;
    if ( VMTKUtils_DistanceToCenterlines(src, lines, &distance) != SV_OK )
    {
        return NULL;
    }

    return distance->GetVtkPolyData();
}

std::vector<sv4guiPathElement::sv4guiPathPoint> sv4guiModelUtils::ConvertToPathPoints(std::vector<mitk::Point3D> posPoints)
{
    std::vector<sv4guiPathElement::sv4guiPathPoint> pathPoints;

    for(int i=0;i<posPoints.size()-1;i++)
    {
        sv4guiPathElement::sv4guiPathPoint pathPoint;
        pathPoint.pos=posPoints[i];
        pathPoint.tangent=posPoints[i+1]-posPoints[i];
        pathPoint.tangent.Normalize();
        pathPoint.rotation=sv4guiMath3::GetPerpendicularNormalVector(pathPoint.tangent);

        pathPoints.push_back(pathPoint);
    }

    sv4guiPathElement::sv4guiPathPoint lastPathPoint=pathPoints.back();
    lastPathPoint.pos=posPoints.back();
    pathPoints.push_back(lastPathPoint);

    return pathPoints;
}

vtkSmartPointer<vtkPolyData> sv4guiModelUtils::GetThresholdRegion(vtkSmartPointer<vtkPolyData> pd, vtkDataObject::FieldAssociations dataType, std::string arrayName, double minValue, double maxValue )
{
    vtkSmartPointer<vtkThreshold> thresholder=vtkSmartPointer<vtkThreshold>::New();
    thresholder->SetInputData(pd);
    thresholder->SetInputArrayToProcess(0, 0, 0, dataType, arrayName.c_str());
    thresholder->ThresholdBetween(minValue, maxValue);
    thresholder->Update();

    vtkSmartPointer<vtkDataSetSurfaceFilter> surfacer=vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfacer->SetInputData(thresholder->GetOutput());
    surfacer->Update();

    return surfacer->GetOutput();
}

std::vector<sv4guiPathElement*> sv4guiModelUtils::CreatePathElements(sv4guiModelElement* modelElement,
                                                             vtkSmartPointer<vtkPolyData> centerlinesPD)
{
    std::vector<sv4guiPathElement*> pathElements;

    if(centerlinesPD==NULL || !centerlinesPD->GetCellData()->HasArray("CenterlineIds"))
        return pathElements;

    int numCenterlines=centerlinesPD->GetCellData()->GetArray("CenterlineIds")->GetRange()[1]+1;

    for(int i=0;i<numCenterlines;i++)
    {
        vtkSmartPointer<vtkPolyData> polyline=GetThresholdRegion(centerlinesPD,vtkDataObject::FIELD_ASSOCIATION_CELLS, "CenterlineIds", i, i);

        vtkSmartPointer<vtkDataArray> groupArray=polyline->GetCellData()->GetArray("GroupIds");
        int lowerValue=groupArray->GetRange()[0];
        int upperValue=groupArray->GetRange()[1];

        vtkSmartPointer<vtkPolyData> centerline=GetThresholdRegion(polyline,vtkDataObject::FIELD_ASSOCIATION_CELLS,"GroupIds",lowerValue,upperValue);
        std::vector<mitk::Point3D> posPoints;
        for(int j=0;j<centerline->GetNumberOfPoints();j++)
        {
            mitk::Point3D point;
            point[0]=centerline->GetPoint(j)[0];
            point[1]=centerline->GetPoint(j)[1];
            point[2]=centerline->GetPoint(j)[2];

            posPoints.push_back(point);
        }

        sv4guiPathElement* pe=new sv4guiPathElement();
        pe->SetMethod(sv4guiPathElement::CONSTANT_TOTAL_NUMBER);
        pe->SetCalculationNumber(centerline->GetNumberOfPoints());
        //pe->SetControlPoints(controlPoints,false);
        pe->SetPathPoints(ConvertToPathPoints(posPoints));

        pathElements.push_back(pe);
    }

    return pathElements;
}

double sv4guiModelUtils::CalculateVpdArea(vtkPolyData* vpd)
{
    if(vpd==NULL)
        return 0;

    double area=0;
    cvPolyData *src=new cvPolyData(vpd);

    if ( sys_geom_SurfArea(src, &area) != SV_OK )
    {
        return 0;
    }

    return area;
}

bool sv4guiModelUtils::CheckPolyDataSurface(vtkPolyData* pd, std::string &msg)
{
    if(pd==NULL)
    {
      msg = "Polydata is empty\n";
      return SV_ERROR;
    }

    pd->BuildLinks();

    int numPolys  = pd->GetNumberOfCells();
    int numPoints = pd->GetNumberOfPoints();

    bool valid=true;
    int numNotTriangles=0;
    int numOpenEdges=0;
    int numNonManifoldEdges=0;
    for (int i=0; i<numPolys; i++)
    {
      vtkIdType npts, *pts;
      pd->GetCellPoints(i, npts, pts);
      if (npts != 3)
      {
        valid = false;
        numNotTriangles++;
      }
      for (int j=0; j<npts; j++)
      {
        vtkIdType p0, p1;
        p0 = pts[j];
        p1 = pts[(j+1)%npts];

        vtkNew(vtkIdList, edgeNeighbor);
        pd->GetCellEdgeNeighbors(i, p0, p1, edgeNeighbor);

        if (edgeNeighbor->GetNumberOfIds() > 1)
          numNonManifoldEdges++;
        else if (edgeNeighbor->GetNumberOfIds() == 0)
          numOpenEdges++;
      }
    }

    msg = "";
    if (!valid)
    {
      msg = msg + "Surface contains non-triangular elements!\n";
      msg = msg + "  Number of non-triangular elements: "+ std::to_string(numNotTriangles) + "\n";
    }
    msg = msg +  "  Number of elements: " + std::to_string(numPolys) + "\n";
    msg = msg +  "  Number of points: " + std::to_string(numPoints) + "\n";
    msg = msg +  "  Number of  non-manifold edges: " + std::to_string(numNonManifoldEdges) + "\n";
    msg = msg +  "  Number of  open edges: " + std::to_string(numOpenEdges) + "\n";

    return valid;
}

bool sv4guiModelUtils::TriangulateSurface(vtkPolyData* pd)
{
  vtkNew(vtkTriangleFilter, triangulator);
  triangulator->SetInputData(pd);
  triangulator->Update();

  pd->DeepCopy(triangulator->GetOutput());

  return SV_OK;
}

vtkPolyData* sv4guiModelUtils::RunDecomposition(sv4guiModelElement* modelElement,
                                                    vtkPolyData *mergedCenterlines)
{
  vtkSmartPointer<vtkPolyData> wallPd=vtkSmartPointer<vtkPolyData>::New();
  wallPd->DeepCopy(modelElement->GetWholeVtkPolyData());

  if(!DeleteRegions(wallPd,modelElement->GetCapFaceIDs()))
  {
    return NULL;
  }

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(wallPd);
  cleaner->Update();
  wallPd->DeepCopy(cleaner->GetOutput());

  cvOCCTSolidModel *cylinderSolid = new cvOCCTSolidModel();
  double radius = 2.0;
  double length = 2.0;
  double ctr[3]; ctr[0] = 0.0; ctr[1] = 0.0; ctr[2] = 0.0;
  double axis[3]; axis[0] = 0.0; axis[1] = 0.0; axis[2] = 1.0;
  cylinderSolid->MakeCylinder(radius, length, ctr, axis);

  TopExp_Explorer seeHowMany;
  seeHowMany.Init(*(cylinderSolid->geom_), TopAbs_EDGE);
  int cylEdges = 0;
  for (int r=0; seeHowMany.More(); seeHowMany.Next(), r++)
  {
    cylEdges++;
    TopoDS_Edge cylEdge = TopoDS::Edge(seeHowMany.Current());

    TopLoc_Location newLoc;
    Standard_Real newFirst, newLast;
    Handle(Geom_Curve) cylCurve = BRep_Tool::Curve (cylEdge, newLoc, newFirst, newLast);
    vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
    for (int k=0; k<20; k++)
    {
      double newVal = (newLast-newFirst)*(k/20.) + newFirst;
      gp_Pnt newPnt = cylCurve->Value(newVal);
      newPoints->InsertNextPoint(newPnt.X(), newPnt.Y(), newPnt.Z());
    }
    vtkSmartPointer<vtkPolyData> newPointsPd = vtkSmartPointer<vtkPolyData>::New();
    newPointsPd->SetPoints(newPoints);
    std::string newFn = "/Users/adamupdegrove/Desktop/tmp/CYLEDGES_"+std::to_string(r)+".vtp";
    vtkSVIOUtils::WriteVTPFile(newFn, newPointsPd);
  }
  fprintf(stdout,"NUM CYL EDEG: %d\n", cylEdges);

  ////=========================== SOSUDDDDDDDDDDDD =======================0l

  if(modelElement==NULL || modelElement->GetWholeVtkPolyData()==NULL)
      return NULL;

  cvPolyData *modelPolyData = new cvPolyData(modelElement->GetWholeVtkPolyData());
  cvPolyData *mergedPolyData = new cvPolyData(mergedCenterlines);
  cvPolyData *decomposedPolyData = NULL;
  std::vector<cvOCCTSolidModel*> loftedSurfs;

  if ( VTKSVUtils_DecomposePolyData(modelPolyData, mergedPolyData, &decomposedPolyData, loftedSurfs) != SV_OK )
  {
    delete modelPolyData;
    delete mergedPolyData;
    if (decomposedPolyData != NULL)
    {
      delete decomposedPolyData;
    }
    return NULL;
  }

  delete modelPolyData;
  delete mergedPolyData;

  if(loftedSurfs.size()==0)
  {
    fprintf(stderr,"No surfs from decomposition\n");
    return NULL;
  }
  gp_Pnt pt0(1.2062324285507, -5.3187279701233, -76.282005310059);
  gp_Pnt pt1(1.1930881738663, 5.764274597168, -76.197998046875);
  TopoDS_Vertex vertex0 = BRepBuilderAPI_MakeVertex(pt0);
  TopoDS_Vertex vertex1 = BRepBuilderAPI_MakeVertex(pt1);

  // ============================ TEST SPLIT EDGE IN HAIFF ================
  for (int surfer=0; surfer<loftedSurfs.size(); surfer++)
  {
    Standard_Real sewtoler =  1.e-6;
    Standard_Real closetoler =  1.e-4;
    ShapeFix_FreeBounds findFree(*(loftedSurfs[surfer]->geom_),sewtoler,closetoler,
              Standard_False,Standard_False);
    TopoDS_Compound freeWires = findFree.GetClosedWires();
    TopExp_Explorer NewEdgeExp;
    NewEdgeExp.Init(freeWires,TopAbs_EDGE);
    for (int i=0;NewEdgeExp.More();NewEdgeExp.Next(),i++)
    {
      //if (surfer == 0)
      //{
        if (i != 0)
        {
          continue;
        }
      //}
      //else
      //{
      //  if (i != 1)
      //  {
      //    continue;
      //  }
      //}


      TopoDS_Edge tmpEdge = TopoDS::Edge(NewEdgeExp.Current());
      GProp_GProps tmpEdgeProps;
      BRepGProp::LinearProperties(tmpEdge,tmpEdgeProps);
      fprintf(stdout,"FULL EDGE PROPS: %.6f\n", tmpEdgeProps.Mass());

      TopLoc_Location tmpEdgeLoc;
      Standard_Real tmpEdgeFirst, tmpEdgeLast;
      Handle(Geom_Curve) tmpCurve = BRep_Tool::Curve (tmpEdge, tmpEdgeLoc, tmpEdgeFirst, tmpEdgeLast);


      // Get closest point on curve
      BRepExtrema_DistShapeShape closestPointFinder0(tmpEdge, vertex0);
      closestPointFinder0.Perform();

      fprintf(stdout,"ACTUAL CLOSE POINT 0: %.6f %.6f %.6f\n", closestPointFinder0.PointOnShape1(1).X(), closestPointFinder0.PointOnShape1(1).Y(), closestPointFinder0.PointOnShape1(1).Z());

      Standard_Real newParam0;
      closestPointFinder0.ParOnEdgeS1(1, newParam0);
      fprintf(stdout,"ACTUAL CLOSE PARAMETER 0: %.6f\n", newParam0);

      BRepExtrema_DistShapeShape closestPointFinder1(tmpEdge, vertex1);
      closestPointFinder1.Perform();

      fprintf(stdout,"ACTUAL CLOSE POINT 1: %.6f %.6f %.6f\n", closestPointFinder1.PointOnShape1(1).X(), closestPointFinder1.PointOnShape1(1).Y(), closestPointFinder1.PointOnShape1(1).Z());
      Standard_Real newParam1;
      closestPointFinder1.ParOnEdgeS1(1, newParam1);
      fprintf(stdout,"ACTUAL CLOSE PARAMETER 1: %.6f\n", newParam1);

      double midPt = (tmpEdgeLast-tmpEdgeFirst)*0.5;

      if (newParam0 > newParam1)
      {
        Standard_Real tmp = newParam1;
        newParam1 = newParam0;
        newParam0 = tmp;
      }

      //TopoDS_Edge edge0 = BRepBuilderAPI_MakeEdge(tmpCurve, tmpEdgeFirst, midPt);
      //TopoDS_Edge edge0 = BRepBuilderAPI_MakeEdge(tmpCurve, tmpEdgeFirst, newParam0);
      TopoDS_Edge edge0 = BRepBuilderAPI_MakeEdge(tmpCurve, 0.0, 0.5);
      GProp_GProps edge0Props;
      BRepGProp::LinearProperties(edge0,edge0Props);
      fprintf(stdout,"EDGE 0 PROPS: %.6f\n", edge0Props.Mass());

      //TopoDS_Edge edge1 = BRepBuilderAPI_MakeEdge(tmpCurve, midPt, tmpEdgeLast);
      //TopoDS_Edge edge1 = BRepBuilderAPI_MakeEdge(tmpCurve, newParam0, newParam1);
      TopoDS_Edge edge1 = BRepBuilderAPI_MakeEdge(tmpCurve, 0.5, 1.0);
      GProp_GProps edge1Props;
      BRepGProp::LinearProperties(edge1,edge1Props);
      fprintf(stdout,"EDGE 1 PROPS: %.6f\n", edge1Props.Mass());

      //TopoDS_Edge edge2 = BRepBuilderAPI_MakeEdge(tmpCurve, newParam1, tmpEdgeLast);
      TopoDS_Edge edge2 = BRepBuilderAPI_MakeEdge(tmpCurve, 0.25, 0.5);
      GProp_GProps edge2Props;
      BRepGProp::LinearProperties(edge2,edge2Props);
      fprintf(stdout,"EDGE 2 PROPS: %.6f\n", edge2Props.Mass());

      // Pull together edges 0 and 2
      BRepBuilderAPI_MakeWire newWiremaker(edge0, edge2);
      newWiremaker.Build();
      TopoDS_Wire halfWire = newWiremaker.Wire();



      //// Fuse together
      //BRepLib_FuseEdges fuser(halfWire, Standard_False);
      //fprintf(stdout,"NUM POINTS TO BE REMOVERRRED: %d\n", fuser.NbVertices());
      //fuser.Perform();
      //TopoDS_Shape fusedShape = fuser.Shape();

      //TopExp_Explorer RealDumbExp;
      //RealDumbExp.Init(fusedShape,TopAbs_EDGE);
      //TopoDS_Edge newestEdge = TopoDS::Edge(RealDumbExp.Current());
      //int numWireEdges = 0;
      //for (int j=0;RealDumbExp.More();RealDumbExp.Next(),j++)
      //{
      //  numWireEdges++;
      //}
      //fprintf(stdout,"SHOUDL BE ONE FUSED EDGE: %d\n", numWireEdges);

      //Handle(Geom_BSplineCurve) curv0BS = OCCTUtils_EdgeToBSpline(edge0);
      //Handle(Geom_BSplineCurve) curv2BS = OCCTUtils_EdgeToBSpline(edge2);

      //Standard_Real aTolV = Precision::Confusion();
      //aTolV = 1.e-3;
      //GeomConvert_CompCurveToBSplineCurve compBS(curv0BS);
      ////Standard_Boolean didWork = compBS.Add(curv2BS, aTolV, Standard_True, Standard_False, 1);
      //Standard_Boolean didWork = compBS.Add(curv2BS, aTolV);
      //Handle(Geom_BSplineCurve) BS = compBS.BSplineCurve();
      //fprintf(stdout,"DID WORK? %d\n", didWork);

      //BRepBuilderAPI_MakeEdge edgeMaker(BS);
      //edgeMaker.Build();
      //TopoDS_Edge newestEdge = edgeMaker.Edge();

      //TopExp_Explorer RealDumbExp;
      //RealDumbExp.Init(newestEdge,TopAbs_EDGE);
      //int numWireEdges = 0;
      //for (int j=0;RealDumbExp.More();RealDumbExp.Next(),j++)
      //{
      //  numWireEdges++;
      //}
      //fprintf(stdout,"SHOUDL BE ONE FUSED EDGE: %d\n", numWireEdges);

      //BRepBuilderAPI_MakeWire wiremaker(edge0, edge1);
      //wiremaker.Build();
      //TopoDS_Wire new2EdgeWire = wiremaker.Wire();
      //BRepBuilderAPI_MakeWire newestWiremaker(edge1, newestEdge);
      BRepBuilderAPI_MakeWire newestWiremaker(edge0, edge1);
      newestWiremaker.Build();
      TopoDS_Wire new2EdgeWire = newestWiremaker.Wire();

      //// ===================== WRITE OUT THE THREE EDGES ===================
      //TopLoc_Location edge0Loc, edge1Loc, edge2Loc;
      //Standard_Real edge0First, edge0Last, edge1First, edge1Last, edge2First, edge2Last;
      //Handle(Geom_Curve) okayCurve0 = BRep_Tool::Curve (newestEdge, edge0Loc, edge0First, edge0Last);
      //Handle(Geom_Curve) okayCurve1 = BRep_Tool::Curve (edge1, edge1Loc, edge1First, edge1Last);
      //fprintf(stdout,"FIRST AND LAST 0: %.6f %.6f\n", edge0First, edge0Last);
      //fprintf(stdout,"FIRST AND LAST 1: %.6f %.6f\n", edge1First, edge1Last);

      //vtkNew(vtkPoints, okayPoints0);
      //vtkNew(vtkPoints, okayPoints1);
      //for (int j=0; j<20; j++)
      //{
      //  double val0 = (edge0Last - edge0First)*(j/20.) + edge0First;
      //  double val1 = (edge1Last - edge1First)*(j/20.) + edge1First;
      //  gp_Pnt okayPt0 = okayCurve0->Value(val0);
      //  gp_Pnt okayPt1 = okayCurve1->Value(val1);

      //  okayPoints0->InsertNextPoint(okayPt0.X(), okayPt0.Y(), okayPt0.Z());
      //  okayPoints1->InsertNextPoint(okayPt1.X(), okayPt1.Y(), okayPt1.Z());
      //}
      //vtkNew(vtkPolyData, okayPointsPd0);
      //okayPointsPd0->SetPoints(okayPoints0);
      //vtkNew(vtkPolyData, okayPointsPd1);
      //okayPointsPd1->SetPoints(okayPoints1);

      //std::string fn0 = "/Users/adamupdegrove/Desktop/tmp/USERCURVE"+std::to_string(surfer)+"_0.vtp";
      //std::string fn1 = "/Users/adamupdegrove/Desktop/tmp/USERCURVE"+std::to_string(surfer)+"_1.vtp";
      //vtkSVIOUtils::WriteVTPFile(fn0, okayPointsPd0);
      //vtkSVIOUtils::WriteVTPFile(fn1, okayPointsPd1);
      //// ===================== WRITE OUT THE THREE EDGES ===================

      // ==============================ANALYZER ================================
      BRepCheck_Analyzer preAnalyzer(*(loftedSurfs[surfer]->geom_), Standard_False);
      fprintf(stdout,"PRE ANALYZER RESULT: %d\n", preAnalyzer.IsValid());

      // =======================================================================

      //reshaper->Replace(tmpEdge,new2EdgeWire,Standard_False);
      //Standard_Integer shapeStatus = reshaper->Status(tmpEdge,new2EdgeWire,Standard_False);
      //TopoDS_Shape newShape;
      //Standard_Integer shapeStatus = reshaper->Status(*(loftedSurfs[surfer]->geom_),newShape,Standard_False);
      //reshaper->Replace(tmpEdge,edge0,Standard_False);
      //Standard_Integer shapeStatus = reshaper->Status(tmpEdge,edge0,Standard_False);
      //fprintf(stdout,"WHAT IS THE STATUS: %d\n", shapeStatus);

      TopExp_Explorer faceExp0;
      faceExp0.Init(*(loftedSurfs[surfer]->geom_), TopAbs_FACE);
      TopoDS_Face face0 = TopoDS::Face(faceExp0.Current());

      Handle(BRepTools_ReShape) reshaper =  new BRepTools_ReShape();
      reshaper->Replace(tmpEdge,new2EdgeWire,Standard_True);
      TopoDS_Shape newShape = reshaper->Apply(*(loftedSurfs[surfer]->geom_));
      //*(loftedSurfs[surfer]->geom_) = newShape;

      //// TRY TO MAKE OWN SHAPEEEEE=====================================
      //Handle(BRepTools_ReShape) fullReshaper = new BRepTools_ReShape();
      //fullReshaper->Replace(face0, newFace, Standard_True);
      //TopoDS_Shape newShape = fullReshaper->Apply(*(loftedSurfs[surfer]->geom_));
      //*(loftedSurfs[surfer]->geom_) = newShape;

      //TopExp_Explorer faceExp1;
      //faceExp1.Init(newShape, TopAbs_FACE);
      //TopoDS_Face face1 = TopoDS::Face(faceExp1.Current());

      BRepCheck_Analyzer postAnalyzer(*(loftedSurfs[surfer]->geom_), Standard_False);
      fprintf(stdout,"POST ANALYZER RESULT: %d\n", postAnalyzer.IsValid());

      BRep_Builder builder;
      TopoDS_Shell shell;
      builder.MakeShell(shell);

      BRepBuilderAPI_NurbsConvert nurbs(face0);
      Handle(Geom_Surface) geom_extrusion = BRepLib_FindSurface(nurbs).Surface();
      Handle(Geom_BSplineSurface) bsplinesurf = Handle(Geom_BSplineSurface)::DownCast(geom_extrusion);

      TopExp_Explorer NewWireExp;
      NewWireExp.Init(freeWires,TopAbs_EDGE);
      int numWires = 0;
      std::vector<TopoDS_Wire> allWires;
      for (int i=0;NewWireExp.More();NewWireExp.Next(),i++)
      {
        TopoDS_Edge toBeWireEdge = TopoDS::Edge(NewWireExp.Current());
        BRepBuilderAPI_MakeWire fWireMaker(toBeWireEdge);
        fWireMaker.Build();

        allWires.push_back(fWireMaker.Wire());
        numWires++;
      }
      fprintf(stdout,"NUM WIRES: %d\n", numWires);

      TopExp_Explorer DumbExp;
      std::vector<TopoDS_Edge> allEdges;
      DumbExp.Init(*(loftedSurfs[surfer]->geom_),TopAbs_EDGE);
      for (int j=0;DumbExp.More();DumbExp.Next(),j++)
      {
        TopoDS_Edge singleEdge = TopoDS::Edge(DumbExp.Current());
        allEdges.push_back(singleEdge);
      }

      if (surfer == -1)
      {
        //OCCTUtils_ShapeFromBSplineSurfaceWithSplitEdges(bsplinesurf, *(loftedSurfs[surfer]->geom_), allEdges);
        OCCTUtils_ShapeFromBSplineSurfaceWithEdges(bsplinesurf, *(loftedSurfs[surfer]->geom_), allEdges);
      }


      //// segmentation of TS
      //Standard_Real Ui1,Ui2,V0,V1;
      //Ui1 = 0;
      //Ui2 = 1;
      //Ui1 = OCCTUtils_PreciseUpar(Ui1, bsplinesurf);
      //Ui2 = OCCTUtils_PreciseUpar(Ui2, bsplinesurf);
      //V0  = bsplinesurf->VKnot(bsplinesurf->FirstVKnotIndex());
      //V1  = bsplinesurf->VKnot(bsplinesurf->LastVKnotIndex());
      //bsplinesurf->Segment(Ui1,Ui2,V0,V1);

      //TopoDS_Face face;
      //builder.MakeFace(face, bsplinesurf, Precision::Confusion());

      //TopoDS_Wire NEWWIIRE;
      //builder.MakeWire(NEWWIIRE);

      //Standard_Real f1, f2, l1, l2;
      //geom_extrusion->Bounds(f1,l1,f2,l2);
      //fprintf(stdout,"WHAT ARE THEY HERE: %.6f %.6f %.6f %.6f\n", f1, l1, f2, l2);

      //TopExp_Explorer DumbExp;
      //DumbExp.Init(*(loftedSurfs[surfer]->geom_),TopAbs_EDGE);
      //int numEdges = 0;
      //for (int j=0;DumbExp.More();DumbExp.Next(),j++)
      //{
      //  numEdges++;
      //}

      //std::vector<TopoDS_Edge> allEdges;
      //DumbExp.Init(*(loftedSurfs[surfer]->geom_),TopAbs_EDGE);
      //for (int j=0;DumbExp.More();DumbExp.Next(),j++)
      //{
      //  TopLoc_Location newLoc;
      //  Standard_Real newFirst, newLast;
      //  TopoDS_Edge theNEWEdge = TopoDS::Edge(DumbExp.Current());
      //  Handle(Geom_Curve) NEWCurve = BRep_Tool::Curve (theNEWEdge, newLoc, newFirst, newLast);

      //  TopoDS_Vertex v0, v1;
      //  TopExp::Vertices(theNEWEdge, v0, v1);

      //  TopoDS_Edge createdEdge;
      //  builder.MakeEdge(createdEdge, NEWCurve, Precision::Confusion());
      //  v0.Orientation(TopAbs_FORWARD);
      //  builder.Add(createdEdge, v0);
      //  v1.Orientation(TopAbs_REVERSED);
      //  builder.Add(createdEdge, v1);
      //  builder.Range(createdEdge, 0, 1);
      //  allEdges.push_back(createdEdge);
      //}

      //DumbExp.Init(*(loftedSurfs[surfer]->geom_),TopAbs_EDGE);
      //for (int j=0;DumbExp.More();DumbExp.Next(),j++)
      //{
      //  TopoDS_Edge theNEWEdge = TopoDS::Edge(DumbExp.Current());
      //  builder.Add(NEWWIIRE, allEdges[j]);
      //}

      //DumbExp.Init(*(loftedSurfs[surfer]->geom_),TopAbs_EDGE);
      //for (int j=0;DumbExp.More();DumbExp.Next(),j++)
      //{
      //  TopoDS_Edge theNEWEdge = TopoDS::Edge(DumbExp.Current());
      //  if (j == 2)
      //  {
      //    builder.UpdateEdge(allEdges[j], new Geom2d_Line(gp_Pnt2d(1,0), gp_Dir2d(0,1)),
      //        face, Precision::Confusion());
      //  }
      //  else
      //  {
      //    builder.UpdateEdge(allEdges[j], new Geom2d_Line(gp_Pnt2d(0,1), gp_Dir2d(1,0)),
      //        face, Precision::Confusion());
      //  }
      //  builder.Range(allEdges[j], face, 0, 1);

      //}
      //fprintf(stdout,"WAS NEW EDGE ADDED: %d\n", numEdges);

      //builder.Add(face, NEWWIIRE);
      //builder.Add(shell, face);

      //*(loftedSurfs[surfer]->geom_) = OCCTUtils_MakeShell(shell);
    }
  }

  // =====================================================================

  //cvOCCTSolidModel* sewSolid=loftedSurfs[0];

  //// CREATE USING SEWING
  //fprintf(stdout,"SEWING\n");
  //cvOCCTSolidModel *sewSolid = new cvOCCTSolidModel();
  //double sewTol = 1.0;
  //sewSolid->Sew(loftedSurfs, sewTol);
  //fprintf(stdout,"SEWED\n");

  // CREATE USING MAKE OF COMPOUND

  cvOCCTSolidModel* sewSolid=loftedSurfs[0];
  BRep_Builder compoundBuilder;
  TopoDS_Compound Compound;
  compoundBuilder.MakeCompound(Compound);
  for (int surfer=0; surfer<loftedSurfs.size(); surfer++)
  {
    compoundBuilder.Add(Compound,*(loftedSurfs[surfer]->geom_));
  }

  *(sewSolid->geom_) = Compound;

  //double sewTol = 1.0;
  //cvOCCTSolidModel* previousSewSolid=NULL;
  ////    SolidModel_SimplifyT smp = SM_Simplify_All;
  //for(int i=1;i<loftedSurfs.size();i++)
  //{
  //  previousSewSolid=sewSolid;
  //  sewSolid=new cvOCCTSolidModel();
  //  if (sewSolid->Sew(loftedSurfs[i],previousSewSolid, sewTol) != SV_OK)
  //  {
  //    MITK_ERROR << "Failed sewing patches together" << endl;
  //    return NULL;
  //  }
  //}

  // Check of shell DOESNT WORK
  //BRepCheck_Shell shellChecker(TopoDS::Shell(*(sewSolid->geom_)));
  //BRepCheck_ListOfStatus status = shellChecker.Status();
  //BRepCheck_ListIteratorOfListOfStatus statit;
  //statit.Initialize(status);
  //for (int j=0; statit.More(); statit.Next(), j++)
  //{
  //  BRepCheck_Status checker = statit.Value();
  //  fprintf(stdout,"WHAT IS SHELL STATUS: %d\n", checker);
  //}
  //
  // LETS TRY TO CREATE SOMETHING OF THIS

  //TopoDS_Shell shell;

  //BRep_Builder B;
  //fprintf(stdout,"MAKE SHELL\n");
  //B.MakeShell(shell);

  TopExp_Explorer thisFaceExp;
  thisFaceExp.Init(*(sewSolid->geom_),TopAbs_FACE);
  TopoDS_Face thisFace = TopoDS::Face(thisFaceExp.Current());
  int faceCount = 0;
  for (int j=0; thisFaceExp.More(); thisFaceExp.Next(), j++)
  {
    faceCount++;
  }
  fprintf(stdout,"NUM FACES: %d\n", faceCount);

  //BRepCheck_ListOfStatus thelist;
  //BRepCheck_DataMapOfShapeListOfStatus myMap;
  //fprintf(stdout,"BIND\n");
  //myMap.Bind(thisFace, thelist);
  //fprintf(stdout,"MAP LIST\n");
  //BRepCheck_ListOfStatus& lst = myMap(thisFace);

  //fprintf(stdout,"T FACE\n");
  //Handle(BRep_TFace)& TF = *((Handle(BRep_TFace)*) &thisFace.TShape());
  //fprintf(stdout,"CHECK NULL\n");
  //if (TF->Surface().IsNull()) {
  //fprintf(stdout,"IS NULL\n");
  //  BRepCheck::Add(lst,BRepCheck_NoSurface);
  //fprintf(stdout,"ADD NADA\n");
  //}
  //else {
  //  // Flag natural restriction???
  //}
  //if (lst.IsEmpty()) {
  //fprintf(stdout,"LIST EMPTY\n");
  //  lst.Append(BRepCheck_NoError);
  //fprintf(stdout,"APPEND NO ERROR\n");
  //}
  //fprintf(stdout,"AQUIIII\n");

  //BRepCheck_Face faceChecker(thisFace);
  //fprintf(stdout,"INTERSECT WIRES: %d\n", faceChecker.IntersectWires());
  //fprintf(stdout,"CLASSIFY WIRES:  %d\n", faceChecker.ClassifyWires());
  //fprintf(stdout,"ORIENTATION OF WIRES: %d\n", faceChecker.OrientationOfWires());
  //fprintf(stdout,"IS UNORIENTABLE: %d\n", faceChecker.IsUnorientable());
  //fprintf(stdout,"GEOMETRIC CONTROLS: %d\n", faceChecker.GeometricControls());

  //fprintf(stdout,"ADD FACE\n");
  //B.Add(shell, thisFace);

  //Standard_Real thissewtoler =  1.e-6;
  //Standard_Real thisclosetoler =  1.e-4;
  //ShapeFix_FreeBounds thisfindFree(*(sewSolid->geom_),thissewtoler,thisclosetoler,
  //          Standard_False,Standard_False);
  //TopoDS_Compound thisfreeWires = thisfindFree.GetClosedWires();
  //TopExp_Explorer thisEdgeExp;
  ////thisEdgeExp.Init(thisfreeWires,TopAbs_WIRE);
  //thisEdgeExp.Init(*(sewSolid->geom_),TopAbs_EDGE);

  //fprintf(stdout,"MAKE WIRE\n");
  //TopoDS_Wire wire;
  //B.MakeWire(wire);

  //for (int j=0; thisEdgeExp.More(); thisEdgeExp.Next(), j++)
  //{
  //  TopoDS_Edge thisEdge = TopoDS::Edge(thisEdgeExp.Current());

  //  //fprintf(stdout,"MAKE WIRE\n");
  //  //B.MakeWire(thisWire);
  //  fprintf(stdout,"ADD EDGE\n");
  //  B.Add(wire, thisEdge);
  //}

  //fprintf(stdout,"ADD WIRE\n");
  //B.Add(thisFace, wire);

  //fprintf(stdout,"CREATE SHELL\n");
  //TopoDS_Shape finalShape = OCCTUtils_MakeShell(shell);
  //*(sewSolid->geom_) = finalShape;

  // NOW CALCULATE HAUSDORFF
  fprintf(stdout,"CALCULATE DISTANCES\n");
  BRepExtrema_DistShapeShape distanceFinder;
  distanceFinder.LoadS1(*(sewSolid->geom_));

  // Average
  vtkNew(vtkDoubleArray, distanceArray);
  distanceArray->SetNumberOfTuples(wallPd->GetNumberOfPoints());
  distanceArray->SetName("Distance");
  double avgDist = 0.0;
  double maxDist = -1.0;
  double minDist = VTK_SV_LARGE_DOUBLE;

  double dist;
  double pt[3];
  TopoDS_Vertex vert;
  gp_Pnt pnt;
  for (int i=0; i<wallPd->GetNumberOfPoints(); i++)
  {
    wallPd->GetPoint(i, pt);
    pnt.SetCoord(pt[0], pt[1], pt[2]);
    vert = BRepBuilderAPI_MakeVertex(pnt);

    distanceFinder.LoadS2(vert);
    if (distanceFinder.Perform() != 1)
    {
      std::cerr << "Finding distance didnt complete." << endl;
      return wallPd;
    }

    dist = distanceFinder.Value();
    distanceArray->SetTuple1(i, dist);

    avgDist += dist;
    if (dist > maxDist)
    {
      maxDist = dist;
    }
    if (dist < minDist)
    {
      minDist = dist;
    }
  }
  avgDist /= wallPd->GetNumberOfPoints();

  wallPd->GetPointData()->AddArray(distanceArray);
  std::string fn = "/Users/adamupdegrove/Desktop/tmp/DISTANCEFILE.vtp";
  vtkSVIOUtils::WriteVTPFile(fn, wallPd);

  fprintf(stdout,"MAX DISTANCE: %.6f\n", maxDist);
  fprintf(stdout,"MIN DISTANCE: %.6f\n", minDist);
  fprintf(stdout,"AVG DISTANCE: %.6f\n", avgDist);


  fprintf(stdout,"GET VTK REP\n");
  //cvPolyData *wholePd = sewSolid->GetPolyData(0, 20.0);
  cvPolyData *wholePd = sewSolid->GetPolyData(1, 5.0);

  return wholePd->GetVtkPolyData();
  //return decomposedPolyData->GetVtkPolyData();

  ////setup face names
  //int numFaces;
  //int *ids;
  //int status=sewSolid->GetFaceIds( &numFaces, &ids);
  //if(status != SV_OK )
  //{
  //    MITK_ERROR << "GetFaceIds: error on object";
  //    return NULL;
  //}

  //std::vector<int> faceIDs;
  //std::vector<std::string> faceNames;
  //for(int i=0;i<numFaces;i++)
  //{
  //    faceIDs.push_back(ids[i]);
  //    char *value=NULL;
  //    char *parent=NULL;
  //    sewSolid->GetFaceAttribute("gdscName",ids[i],&value);
  //    std::string type(value);
  //    sewSolid->GetFaceAttribute("parent",ids[i],&parent);
  //    std::string groupName(parent);
  //    faceNames.push_back(type+"_"+groupName);
  //}

  //for(int i=0;i<faceNames.size()-1;i++)
  //{
  //    int idx=1;
  //    for(int j=i+1;j<faceNames.size();j++)
  //    {
  //        if(faceNames[i]==faceNames[j])
  //        {
  //            idx++;
  //            std::stringstream ss;
  //            ss << idx;
  //            std::string idxStr = ss.str();
  //            faceNames[j]=faceNames[j]+"_"+idxStr;
  //        }
  //    }
  //}

  //for(int i=0;i<numFaces;i++)
  //{
  //    char* fn=const_cast<char*>(faceNames[i].c_str());

  //    sewSolid->SetFaceAttribute("gdscName",ids[i],fn);
  //}

  //double maxDist = 1.0;
  //cvPolyData* cvwholevpd=sewSolid->GetPolyData(0,maxDist);
  //if(cvwholevpd==NULL || cvwholevpd->GetVtkPolyData()==NULL)
  //    return NULL;

  //std::vector<sv4guiModelElement::svFace*> faces;

  //for(int i=0;i<numFaces;i++)
  //{
  //    cvPolyData* cvfacevpd=sewSolid->GetFacePolyData(ids[i],1,maxDist);
  //    if(cvfacevpd==NULL || cvfacevpd->GetVtkPolyData()==NULL)
  //        return NULL;

  //    sv4guiModelElement::svFace* face =new sv4guiModelElement::svFace;
  //    face->id=ids[i];
  //    face->name=faceNames[i];
  //    face->vpd=cvfacevpd->GetVtkPolyData();

  //    if(face->name.substr(0,5)=="wall_")
  //        face->type="wall";
  //    else if(face->name.substr(0,4)=="cap_")
  //        face->type="cap";

  //    faces.push_back(face);
  //}

  //sv4guiModelElement* newModelElement=new sv4guiModelElement();
  ////newModelElement->SetSegNames(segNames);
  //newModelElement->SetFaces(faces);
  //newModelElement->SetWholeVtkPolyData(cvwholevpd->GetVtkPolyData());
  ////newModelElement->SetNumSampling(numSamplingPts);
  //newModelElement->SetInnerSolid(sewSolid);
  ////newModelElement->SetMaxDist(maxDist);

  //return newModelElement;
}

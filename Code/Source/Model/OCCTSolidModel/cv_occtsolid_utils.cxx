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

/** @file cv_occtsolid_utils.cxx
 *  @brief The implementations of functions in cv_polydatasolid_utils
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */
#include "SimVascular.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "cv_occtsolid_utils.h"

//OCCT Includes
#include "Precision.hxx"
#include "TopoDS.hxx"
#include "TopoDS_Face.hxx"
#include "TopoDS_Shell.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_Array1OfShape.hxx"
#include "TopTools_ListIteratorOfListOfShape.hxx"

#include "BSplCLib.hxx"
#include "BRepOffsetAPI_ThruSections.hxx"
#include "BRepClass3d_SolidClassifier.hxx"
#include "BRepTools_ReShape.hxx"
#include "BRepCheck_Analyzer.hxx"
#include "BRepCheck_Solid.hxx"
#include "BRepCheck_ListOfStatus.hxx"
#include "BRepCheck_ListIteratorOfListOfStatus.hxx"
#include "BRepExtrema_DistShapeShape.hxx"
#include "BRep_Tool.hxx"
#include "BRep_Builder.hxx"
#include "BRepFill_Filling.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "BRepBuilderAPI_FindPlane.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "BRepBuilderAPI_MakeShell.hxx"
#include "BRepBuilderAPI_MakeSolid.hxx"
#include "BRepBuilderAPI_MakeVertex.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_Sewing.hxx"
#include "BRepFilletAPI_MakeFillet.hxx"

//#include "ChFi3D_FilletShape.hxx"

#include "TDataStd_Integer.hxx"
#include "TDataStd_Name.hxx"
#include "Standard_Real.hxx"
#include "Standard_NullObject.hxx"
#include "StdFail_NotDone.hxx"
#include "Standard_Integer.hxx"
#include "TDataStd_Integer.hxx"
#include "TNaming_Builder.hxx"

#include "ShapeAnalysis_Surface.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "ShapeFix_FreeBounds.hxx"
#include "ShapeFix_Shape.hxx"

#include "GeomAPI_Interpolate.hxx"
#include "Geom_Plane.hxx"
#include "GeomFill_Line.hxx"
#include "GeomFill_AppSurf.hxx"
#include "GeomFill_SectionGenerator.hxx"
#include "GeomConvert.hxx"
#include "GeomConvert_ApproxCurve.hxx"
#include "GeomConvert_CompCurveToBSplineCurve.hxx"
#include "Geom_BSplineSurface.hxx"
#include "Geom_BoundedCurve.hxx"
#include "Geom_TrimmedCurve.hxx"
#include "Geom2d_Line.hxx"
#include "Geom_Conic.hxx"
#include "GeomLProp_SLProps.hxx"
#include "GCPnts_UniformAbscissa.hxx"

#include "vtkSmartPointer.h"
#include "vtkSVIOUtils.h"

#include <string>
#include <sstream>
#include <iostream>

//Function to turn an integer into a string
int charToInt(const char *in)
{
  int out = atoi(in);
  return out;
}

// ---------------------
// OCCTUtils_SetExtStringArray
// ---------------------
/**
 * @brief Helper function to create an occt string array from a character
 * array. Used to set labels on faces
 * @param &array the new string array in which to put the char array
 * @param &charstr the character to put in array
 * @return SV_OK if function completes properly
 */
int OCCTUtils_SetExtStringArrayFromChar(Handle(TDataStd_ExtStringArray) &array,
    char *charstr)
{
  int lower = 0;
  int upper = strlen(charstr);
  array->Init(lower,upper);
  for (int i=lower;i<upper;i++)
  {
    array->SetValue(i,charstr[i]);
  }

  return SV_OK;
}

// ---------------------
// OCCTUtils_GetExtStringArray
// ---------------------
/**
 * @brief Helper function get a character array from an occt string array.
 * @param &array the string array to be extracted
 * @param &charstr the char array to put in the information
 * @return SV_OK if function completes properly
 */
int OCCTUtils_GetExtStringArrayAsChar(Handle(TDataStd_ExtStringArray) &array,
    char *charstr)
{
    std::stringstream streamer;
    for (int i=array->Lower();i<array->Upper();i++)
    {
      TCollection_AsciiString asciiString(array->Value(i),'?');
      streamer << asciiString.ToCString();
    }
    std::string outstr = streamer.str();
    sprintf(charstr,"%s",outstr.c_str());

    return SV_OK;
}

// ---------------------
// Class LinearFunc
// ---------------------
/**
 * @brief Simple class to set up linear interpolation between points
 */
class LinearFunc
{
  public:
    LinearFunc(Standard_Real x1,Standard_Real y1,
	Standard_Real x2,Standard_Real y2);
    ~LinearFunc() {;}

    Standard_Real GetY(Standard_Real x);
  private:
    Standard_Real m_;
    Standard_Real b_;
};

//Constructor
LinearFunc::LinearFunc(Standard_Real x1,Standard_Real y1,
    		       Standard_Real x2,Standard_Real y2)
{
  m_ = (y2 - y1)/(x2 - x1);
  b_ = y2 - (m_*x2);
}

//Get output linear interpolated value
Standard_Real LinearFunc::GetY(Standard_Real x)
{
  Standard_Real y = m_*x + b_;
  return y;
}

// ---------------------
// OCCTUtils_CreateEdgeBlend
// ---------------------
/**
 * @brief Procedure to create edge blend between two faces faceA and faceB
 * @param &shape shape containing faces and resultant shape with blend
 * @param shapetool the XDEDoc manager that contains attribute info
 * @param shapelabel the label for the shape registered in XDEDoc
 * @param filletmake the occt API to make a fillet
 * @param faceA first integer face to blend
 * @param faceB second integer face to blend
 * @param radius Maximum radius to set anywhere on the fillet.
 * @param minRadius Minimum radius to set anywhere on the fillet. A linear
 * interpolation is created between the minimum and maximum radius specified
 * based on the angle created bewteen the two faces at a set number of points
 * around the fillet edge. The new fillet radius value will be somehwere
 * between the maximum and minimum radius values given.
 * @param blendname Name to be given the new face created for the shape
 * @return SV_OK if function completes properly
 */
int OCCTUtils_CreateEdgeBlend(TopoDS_Shape &shape,
		Handle(XCAFDoc_ShapeTool) &shapetool,TDF_Label &shapelabel,
		BRepFilletAPI_MakeFillet &filletmaker,
		int faceA,int faceB,double radius,
		char blendname[])
{
  //if (minRadius > radius)
  //  fprintf(stderr,"Minimum radius is larger than Maximum radius\n");

  char *name;
  char nameA[255];
  char nameB[255];
  TopTools_IndexedDataMapOfShapeListOfShape anEFsMap;
  TopExp::MapShapesAndAncestors (shape, TopAbs_EDGE,
      TopAbs_FACE, anEFsMap);
  int num = anEFsMap.Extent();
  int found = 0;

  TopoDS_Shape geomcopy = shape;
  TopoDS_Edge filletEdge;
  for (int i=1;i < num+1;i++)
  {
    TopTools_ListOfShape faces = anEFsMap.FindFromIndex(i);
    TopoDS_Shape face1 = faces.First();
    int faceId1,faceId2;
    OCCTUtils_GetFaceLabel(face1,shapetool,shapelabel,faceId1);
    TopoDS_Shape face2 = faces.Last();
    OCCTUtils_GetFaceLabel(face2,shapetool,shapelabel,faceId2);
    if ((faceId1 == faceA && faceId2 == faceB) ||
		    (faceId1 == faceB && faceId2 == faceA))
    {
      found++;
      char* checkid1;
      char* checkid2;
      int intcheck1,intcheck2;
      filletEdge = TopoDS::Edge(anEFsMap.FindKey(i));
      OCCTUtils_GetFaceAttribute(face1,shapetool,shapelabel,
	  "gdscName",&name);
      strncpy(nameA,name,sizeof(nameA));
      OCCTUtils_GetFaceAttribute(face2,shapetool,shapelabel,
	  "gdscName",&name);
      strncpy(nameB,name,sizeof(nameB));

      ////Trying to get a curvature measure
      //BRepAdaptor_Curve curveAdaptor;
      //curveAdaptor.Initialize(filletEdge);
      //Standard_Integer NbPts = 20;
      //GCPnts_UniformAbscissa uAbs;
      //uAbs.Initialize(curveAdaptor,NbPts);

      //if (!uAbs.IsDone())
      //{
      //  fprintf(stderr,"Could not create points on edge\n");
      //  return SV_ERROR;
      //}
      //TopoDS_Face face1ForCurve = TopoDS::Face(faces.First());
      //TopoDS_Face face2ForCurve = TopoDS::Face(faces.Last());
      //Handle(Geom_Surface) surfHand1 =
      //        BRep_Tool::Surface(face1ForCurve);
      //Handle(Geom_Surface) surfHand2 =
      //        BRep_Tool::Surface(face2ForCurve);
      //ShapeAnalysis_Surface analyzer1(surfHand1);
      //ShapeAnalysis_Surface analyzer2(surfHand2);
      //TColgp_Array1OfPnt2d radLawArray(1,uAbs.NbPoints());
      //Standard_Real minAng = M_PI + 1.0;
      //Standard_Real maxAng = -1.0;
      //for (int j=1;j<=uAbs.NbPoints();j++)
      //{
      //  gp_Pnt nextPnt = curveAdaptor.Value(uAbs.Parameter(j));
      //  gp_Pnt2d face1UV = analyzer1.ValueOfUV(nextPnt,1.0e-6);
      //  gp_Pnt2d face2UV = analyzer2.ValueOfUV(nextPnt,1.0e-6);
      //  GeomLProp_SLProps prop1(surfHand1,face1UV.X(),face1UV.Y(),1,1.0e-6);
      //  GeomLProp_SLProps prop2(surfHand2,face2UV.X(),face2UV.Y(),1,1.0e-6);
      //  fprintf(stderr,"Curvature 1 check %.4f\n",prop1.MeanCurvature());
      //  fprintf(stderr,"Curvature 2 check %.4f\n",prop2.MeanCurvature());
      //  gp_Vec f1tan1 = prop1.D1U();
      //  gp_Vec f1tan2 = prop1.D1V();
      //  gp_Vec f2tan1 = prop2.D1U();
      //  gp_Vec f2tan2 = prop2.D1V();
      //  gp_Vec norm1 = f1tan1.Crossed(f1tan2);
      //  gp_Vec norm2 = f2tan1.Crossed(f2tan2);
      //  Standard_Real ang = norm1.Angle(norm2);
      //  fprintf(stderr,"Angle between face normals at point is: %.4f\n",ang);
      //  gp_Pnt2d angSet(uAbs.Parameter(j),ang);
      //  radLawArray.SetValue(j,angSet);
      //  if (ang < minAng)
      //    minAng = ang;
      //  if (ang > maxAng)
      //    maxAng = ang;
      //}
      //fprintf(stderr,"Max angle: %.4f\n",maxAng);
      //fprintf(stderr,"Min angle: %.4f\n",minAng);
      //LinearFunc radLaw(minAng,minRadius,maxAng,radius);
      //fprintf(stderr,"Checking worked max radius: %.4f\n",radLaw.GetY(maxAng));
      //fprintf(stderr,"Checking worked min radius: %.4f\n",radLaw.GetY(minAng));
      //for (int j=1;j<=uAbs.NbPoints();j++)
      //{
      //  gp_Pnt2d currVal = radLawArray.Value(j);
      //  Standard_Real newRad = radLaw.GetY(currVal.Y());
      //  gp_Pnt2d newVal(currVal.X(),newRad);
      //  radLawArray.SetValue(j,newVal);
      //  fprintf(stderr,"%d point has new radius value of %.4f for parameter %.4f\n",j,newVal.Y(),newVal.X());
      //}
      //filletmaker.Add(radLawArray,filletEdge);
      filletmaker.Add(radius,filletEdge);
    }
  }

  fprintf(stderr,"Number of edges %d\n",found);
  if (found == 0)
  {
    fprintf(stderr,"No edges between faces\n");
    return SV_ERROR;
  }
  sprintf(blendname,"wall_blend_%s_%s",nameA,nameB);

  TopoDS_Shape tmpShape;
  try
  {
    filletmaker.Build();
  }
  catch (Standard_Failure)
  {
    fprintf(stderr,"Try different radius\n");
    return SV_ERROR;
  }
  try
  {
    tmpShape = filletmaker.Shape();
  }
  catch (StdFail_NotDone)
  {
    fprintf(stderr,"Try different radius\n");
    return SV_ERROR;
  }
  if (filletmaker.IsDone() != 1)
  {
    fprintf(stderr,"Not done\n");
    return SV_ERROR;
  }
  fprintf(stderr,"Number faulty contours: %d\n",filletmaker.NbFaultyContours());
  shape = tmpShape;

  return SV_OK;
}
// ---------------------
// OCCTUtils_CapShapeToSolid
// ---------------------
/**
 * @brief
 * @param
 * @return SV_OK if function completes properly
 */
int OCCTUtils_CapShapeToSolid(TopoDS_Shape &shape,TopoDS_Shape &geom,
    BRepBuilderAPI_Sewing &attacher,int &numFilled)
{
  if (shape.Closed())
  {
    geom = shape;
    fprintf(stdout,"Shape is closed, nothing to be done\n");
    return SV_OK;
  }

  //Attacher!
  attacher.Add(shape);
  Standard_Real sewtoler =  1.e-6;
  Standard_Real closetoler =  1.e-4;
  ShapeFix_FreeBounds findFree(shape,sewtoler,closetoler,
        	  Standard_False,Standard_False);
  TopoDS_Compound freeWires = findFree.GetClosedWires();
  TopExp_Explorer NewEdgeExp;
  NewEdgeExp.Init(freeWires,TopAbs_EDGE);
  for (int i=0;NewEdgeExp.More();NewEdgeExp.Next(),i++)
  {
    TopoDS_Edge tmpEdge = TopoDS::Edge(NewEdgeExp.Current());

    BRepBuilderAPI_MakeWire wiremaker(tmpEdge);
    wiremaker.Build();

    int degOfCap=2;
    int numPointsOnCurve = 20;
    BRepFill_Filling filler(degOfCap,numPointsOnCurve);
    filler.Add(tmpEdge,GeomAbs_C0,Standard_True);

    try {
      filler.Build();
    }
    catch (Standard_Failure)
    {
      fprintf(stderr,"Failure when filling holes\n");
      return SV_ERROR;
    }

    TopoDS_Face newFace = filler.Face();
    //if (!OCCTUtils_IsSameOrientedWEdge(newFace,shape,tmpEdge))
    //{
    //  fprintf(stderr,"Reversing\n");
    //  newFace.Reversed();
    //}
    attacher.Add(newFace);
    numFilled ++;
  }
  fprintf(stderr,"Number of holes found: %d\n",numFilled);
  if (numFilled == 0)
  {
    fprintf(stderr,"No holes found\n");
    return SV_OK;
  }
  attacher.Perform();

  TopoDS_Shell tmpShell;
  try {
    tmpShell = TopoDS::Shell(attacher.SewedShape());
  }
  catch (Standard_TypeMismatch) {
    fprintf(stderr,"No open boundaries found\n");
    return SV_ERROR;
  }

  BRepBuilderAPI_MakeSolid solidmaker(tmpShell);
  geom = solidmaker.Solid();
  geom.Closed(Standard_True);

  //BRep_Builder BB;
  //TopoDS_Solid solid;
  //BB.MakeSolid(solid);
  //BB.Add(solid,tmpShell);

  ////Set Orientation
  //BB.MakeSolid(solid);
  //BRepClass3d_SolidClassifier clas3d(solid);
  //clas3d.PerformInfinitePoint(Precision::Confusion());
  //fprintf(stderr,"Print state: %d\n",clas3d.State());
  //if (clas3d.State() == TopAbs_IN) {
  //  BB.MakeSolid(solid);
  //  TopoDS_Shape aLocalShape = tmpShell.Reversed();
  //  BB.Add(solid, TopoDS::Shell(aLocalShape));
  //}

  //solid.Closed(Standard_True);
  //geom = solid;
  return SV_OK;
}

// ---------------------
// OCCTUtils_MakeLoftedSurf
// ---------------------
/**
 * @brief Procedure to create a lofted surface from a set of wires
 * @param curves list of wires to be lofted. Wires should be created by
 * uniformly spaced points circumferentially and should be aligned with
 * first points matching
 * @param shape place to store the output lofted surface
 * @param numCurves number of wires contained in curves
 * @param continuity desired continuity of the ouptut surface. Typically, 2
 * works okay. Trouble with higher or lower
 * @param partype The parametrization method. 0 is chord, 1 is centripetal,
 * 2 is isoparametric. Depends on wires which works best.
 * @param w1, first weighting to be used only if smoothing
 * @param w2, second weighting to be used only if smoothing
 * @param w3, third weighting to be used only if smoothing
 * @param smoothing indicates whether smoothing should be used
 * @return SV_OK if function completes properly
 */
int OCCTUtils_MakeLoftedSurf(TopoDS_Wire *curves, TopoDS_Shape &shape,
		int numCurves,int continuity,
		int partype, double w1, double w2, double w3, int smoothing)
{
  //Methods using GeomFill_SectionGenerator
  GeomFill_SectionGenerator sectioner;
  Handle(Geom_BSplineSurface) surface;
  Handle(Geom_BSplineSurface) tmpSurface;
  Handle(Geom_BSplineCurve) BS, BS1;
  Handle(Geom_TrimmedCurve) curvTrim;

  Standard_Boolean checkDegenerate = Standard_False;
  for (int i = 0; i< numCurves;i++)
  {
    TopExp_Explorer getEdge(curves[i],TopAbs_EDGE);
    TopoDS_Edge tmpEdge = TopoDS::Edge(getEdge.Current());

    checkDegenerate = BRep_Tool::Degenerated(tmpEdge);
    if (checkDegenerate == Standard_True)
    {
      fprintf(stderr,"Degenerate wire detected\n");
      return SV_ERROR;
    }
    Handle(Geom_BSplineCurve) curvBS = OCCTUtils_EdgeToBSpline(tmpEdge);

    Standard_Real aTolV = Precision::Confusion();
    aTolV = 1.e-3;
    GeomConvert_CompCurveToBSplineCurve compBS(curvBS);
    compBS.Add(curvBS,aTolV,Standard_True,Standard_False,1);
    BS = compBS.BSplineCurve();
    sectioner.AddCurve(BS);
  }

  sectioner.Perform(Precision::PConfusion());
  Handle(GeomFill_Line) line = new GeomFill_Line(numCurves);

  Standard_Real pres3d = 1.e-6;
  Standard_Integer nbIt = 3;
  if(pres3d <= 1.e-3) nbIt = 0;

  Standard_Integer degmin = 2, degmax = 2;//Max(myDegMax, degmin);
  Standard_Boolean SpApprox = Standard_True;

  GeomFill_AppSurf anApprox(degmin, degmax, pres3d, pres3d, nbIt);
  anApprox.SetContinuity((GeomAbs_Shape) continuity);

  //anApprox.SetCriteriumWeight(w1, w2, w3);
  if(smoothing) {
    anApprox.SetCriteriumWeight(w1, w2, w3);
    try
    {
      anApprox.PerformSmoothing(line, sectioner);
    }
    catch (std::bad_alloc)
    {
      fprintf(stderr,"Not enough memory for this smoothing\n");
      return SV_ERROR;
    }
  }
  else
  {
    anApprox.SetParType((Approx_ParametrizationType) partype);
    anApprox.Perform(line, sectioner, SpApprox);
  }

  if(anApprox.IsDone()) {
    //fprintf(stderr,"UDegree %d\n",anApprox.UDegree());
    //fprintf(stderr,"VDegree %d\n",anApprox.VDegree());
    //TColStd_Array2OfReal surfweights = anApprox.SurfWeights();
    //fprintf(stderr,"RowLength %d\n",surfweights.RowLength());
    //fprintf(stderr,"ColLength %d\n",surfweights.ColLength());
    //fprintf(stderr,"SurfWeights\n");
    //for (int i=1;i<surfweights.ColLength();i++)
    //{
    //  for (int j=1;j<surfweights.RowLength();j++)
    //  {
    //    fprintf(stderr,"%.2f ",surfweights.Value(i,j));
    //  }
    //  fprintf(stderr,"\n");
    //}
    //TColStd_Array1OfReal uknots = anApprox.SurfUKnots();
    //TColStd_Array1OfReal vknots = anApprox.SurfVKnots();
    //fprintf(stderr,"Uknot length %d\n",uknots.Length());
    //for (int i=1;i<=uknots.Length();i++)
    //{
    //  fprintf(stderr,"%.8f ",uknots.Value(i));
    //}
    //fprintf(stderr,"\n");
    //fprintf(stderr,"Vknot length %d\n",vknots.Length());
    //for (int i=1;i<=vknots.Length();i++)
    //{
    //  fprintf(stderr,"%.2f ",vknots.Value(i));
    //}
    //fprintf(stderr,"\n");
    //TColStd_Array1OfInteger umults = anApprox.SurfUMults();
    //TColStd_Array1OfInteger vmults = anApprox.SurfVMults();
    //fprintf(stderr,"Umult length %d\n",umults.Length());
    //for (int i=1;i<=umults.Length();i++)
    //{
    //  fprintf(stderr,"%d ",umults.Value(i));
    //}
    //fprintf(stderr,"\n");
    //fprintf(stderr,"Vmult length %d\n",vmults.Length());
    //for (int i=1;i<=vmults.Length();i++)
    //{
    //  fprintf(stderr,"%d ",vmults.Value(i));
    //}
    //fprintf(stderr,"\n");
    surface =
      new Geom_BSplineSurface(anApprox.SurfPoles(), anApprox.SurfWeights(),
      anApprox.SurfUKnots(), anApprox.SurfVKnots(),
      anApprox.SurfUMults(), anApprox.SurfVMults(),
      anApprox.UDegree(), anApprox.VDegree());
    //surface->SetUPeriodic();
    //fprintf(stdout,"-----------------BSPLINE PARAMETERS----------------------\n");
    //fprintf(stdout,"U Degree:             %d\n",surface->UDegree());
    //fprintf(stdout,"Is U Closed?:         %d\n",surface->IsUClosed());
    //fprintf(stdout,"Is U Periodic?:       %d\n",surface->IsUPeriodic());
    //fprintf(stdout,"Is U Rational?:       %d\n",surface->IsURational());
    //fprintf(stdout,"Nb U Poles:           %d\n",surface->NbUPoles());
    //fprintf(stdout,"Nb U Knots:           %d\n",surface->NbUKnots());
    //fprintf(stdout,"First U Knot Index:   %d\n",surface->FirstUKnotIndex());
    //fprintf(stdout,"Last U Knot Index:    %d\n",surface->LastUKnotIndex());
    //fprintf(stdout,"_________________________________________________________\n");
    //fprintf(stdout,"V Degree:             %d\n",surface->VDegree());
    //fprintf(stdout,"Is V Closed?:         %d\n",surface->IsVClosed());
    //fprintf(stdout,"Is V Periodic?:       %d\n",surface->IsVPeriodic());
    //fprintf(stdout,"Is U Rational?:       %d\n",surface->IsVRational());
    //fprintf(stdout,"Nb V Poles:           %d\n",surface->NbVPoles());
    //fprintf(stdout,"Nb V Knots:           %d\n",surface->NbVKnots());
    //fprintf(stdout,"First V Knot Index:   %d\n",surface->FirstVKnotIndex());
    //fprintf(stdout,"Last V Knot Index:    %d\n",surface->LastVKnotIndex());
    //fprintf(stdout,"_________________________________________________________\n");
  }

  if (OCCTUtils_ShapeFromBSplineSurface(surface,shape,curves[0],curves[numCurves-1]) != SV_OK)
  {
    fprintf(stderr,"Error in conversion from bspline surface to shape\n");
    return SV_ERROR;
  }

  return SV_OK;
}

// ---------------------
// OCCTUtils_ShapeFromBSplineSurfaceWithSplitEdges
// ---------------------
int OCCTUtils_ShapeFromBSplineSurfaceWithSplitEdges(const Handle(Geom_BSplineSurface) surface,
    		TopoDS_Shape &shape, std::vector<TopoDS_Edge> edges)
{
  if(surface.IsNull()) {
    fprintf(stderr,"Lofting did not complete\n");
    return SV_ERROR;
  }
  fprintf(stdout,"IS U PERIODIC? %d\n", surface->IsUPeriodic());
  fprintf(stdout,"IS V PERIODIC? %d\n", surface->IsVPeriodic());

  // WRITE THE EDGES
  for (int j=0; j<edges.size(); j++)
  {
    TopLoc_Location newLoc;
    Standard_Real newFirst, newLast;
    Handle(Geom_Curve) newCurve = BRep_Tool::Curve (edges[j], newLoc, newFirst, newLast);
    fprintf(stdout,"IS PERIODIC? %d\n", newCurve->IsPeriodic());

    vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
    for (int k=0; k<20; k++)
    {
      double newVal = (newLast-newFirst)*(k/20.) + newFirst;
      gp_Pnt newPnt = newCurve->Value(newVal);
      newPoints->InsertNextPoint(newPnt.X(), newPnt.Y(), newPnt.Z());
    }
    vtkSmartPointer<vtkPolyData> newPointsPd = vtkSmartPointer<vtkPolyData>::New();
    newPointsPd->SetPoints(newPoints);
    std::string newFn = "/Users/adamupdegrove/Desktop/tmp/REALEDGE_"+std::to_string(j)+".vtp";
    vtkSVIOUtils::WriteVTPFile(newFn, newPointsPd);
  }

  int numEdges = edges.size();

  // create the new surface
  TopoDS_Shell shell;
  TopoDS_Face face;
  TopoDS_Wire W;
  TopoDS_Edge edge, couture;
  std::vector<TopoDS_Edge> theseEdges(numEdges);
  std::vector<Handle(Geom_Curve)> theseCurves(numEdges);

  for (int j=0; j<numEdges; j++)
  {
    TopLoc_Location newLoc;
    Standard_Real newFirst, newLast;
    Handle(Geom_Curve) newCurve = BRep_Tool::Curve (edges[j], newLoc, newFirst, newLast);
    fprintf(stdout,"NEWFIRST: %.6f, NEWLAST: %.6f\n", newFirst, newLast);
    theseCurves[j] = newCurve;
  }

  BRep_Builder builder;
  builder.MakeShell(shell);

  std::vector<TopoDS_Vertex> vf(numEdges);
  std::vector<TopoDS_Vertex> vl(numEdges);
  std::vector<gp_Pnt> vpf(numEdges);
  std::vector<gp_Pnt> vpl(numEdges);

  for (int j=0; j<numEdges; j++)
  {
    TopExp::Vertices(edges[j], vf[j], vl[j]);
    vpf[j] = BRep_Tool::Pnt(vf[j]);
    vpl[j] = BRep_Tool::Pnt(vl[j]);
  }

  fprintf(stdout,"CHECK VERTICES: %.6f %.6f %.6f\n", vpf[0].X(), vpf[0].Y(), vpf[0].Z());
  fprintf(stdout,"CHECK VERTICES: %.6f %.6f %.6f\n", vpl[0].X(), vpl[0].Y(), vpl[0].Z());
  fprintf(stdout,"CHECK VERTICES: %.6f %.6f %.6f\n", vpf[1].X(), vpf[1].Y(), vpf[1].Z());
  fprintf(stdout,"CHECK VERTICES: %.6f %.6f %.6f\n", vpl[1].X(), vpl[1].Y(), vpl[1].Z());

  // make the face
  builder.MakeFace(face, surface, Precision::Confusion());

  // make the wire
  builder.MakeWire(W);

  // make the missing edges
  Standard_Real f1, f2, l1, l2;
  surface->Bounds(f1,l1,f2,l2);
  fprintf(stdout,"WHAT ARE THESE: %.6f %.6f %.6f %.6f\n", f1, l1, f2, l2);

  // --- edge 1
  builder.MakeEdge(theseEdges[0], theseCurves[0], Precision::Confusion());
  vf[0].Orientation(TopAbs_FORWARD);
  builder.Add(theseEdges[0], vf[0]);
  vl[0].Orientation(TopAbs_REVERSED);
  builder.Add(theseEdges[0], vl[0]);
  builder.Range(theseEdges[0], 0, 0.5);
  // processing of looping sections
  // store edges of the 1st section

  //// --- edge 2
  builder.MakeEdge(theseEdges[1], theseCurves[1], Precision::Confusion());
  vf[1].Orientation(TopAbs_FORWARD);
  builder.Add(theseEdges[1], vf[1]);
  vl[1].Orientation(TopAbs_REVERSED);
  builder.Add(theseEdges[1], vl[1]);
  builder.Range(theseEdges[1], 0.5, 1);
  // processing of looping sections
  // store edges of the 1st section

  // --- edge 3
  builder.MakeEdge(theseEdges[2], theseCurves[2], Precision::Confusion());
  vf[2].Orientation(TopAbs_FORWARD);
  builder.Add(theseEdges[2], vf[2]);
  vl[2].Orientation(TopAbs_REVERSED);
  builder.Add(theseEdges[2], vl[2]);
  builder.Range(theseEdges[2], f2, l2);
  couture = theseEdges[2];

  // --- edge 5
  theseEdges[4] = couture;
  theseEdges[4].Reverse();

  // --- edge 4
  builder.MakeEdge(theseEdges[3], theseCurves[3], Precision::Confusion());
  vf[3].Orientation(TopAbs_FORWARD);
  builder.Add(theseEdges[3], vf[3]);
  vl[3].Orientation(TopAbs_REVERSED);
  builder.Add(theseEdges[3], vl[3]);
  builder.Range(theseEdges[3], f1, l1);
  theseEdges[3].Reverse();

  builder.Add(W,theseEdges[0]);
  builder.Add(W,theseEdges[1]);
  builder.Add(W,theseEdges[2]);
  builder.Add(W,theseEdges[3]);
  builder.Add(W,theseEdges[4]);

  // set PCurve
  builder.UpdateEdge(theseEdges[0],new Geom2d_Line(gp_Pnt2d(0,0),gp_Dir2d(1,0)),face,
    Precision::Confusion());
  builder.Range(theseEdges[0],face,0,1);

  builder.UpdateEdge(theseEdges[1],new Geom2d_Line(gp_Pnt2d(0.5,1),gp_Dir2d(1,0)),face,
    Precision::Confusion());
  builder.Range(theseEdges[1],face,0.5,1);

  builder.UpdateEdge(theseEdges[3],new Geom2d_Line(gp_Pnt2d(0,l2),gp_Dir2d(1,0)),face,
    Precision::Confusion());
  builder.Range(theseEdges[3],face,f1,l1);

  builder.UpdateEdge(theseEdges[4],
    new Geom2d_Line(gp_Pnt2d(l1,0),gp_Dir2d(0,1)),
    new Geom2d_Line(gp_Pnt2d(f1,0),gp_Dir2d(0,1)),face,
    Precision::Confusion());
  builder.Range(theseEdges[4],face,f2,l2);

  //builder.UpdateEdge(edges[0],new Geom2d_Line(gp_Pnt2d(0,f2),gp_Dir2d(1,0)),face,
  //  Precision::Confusion());
  //builder.Range(edges[0],face,f1,l1);
  //builder.UpdateEdge(edges[2],new Geom2d_Line(gp_Pnt2d(0,l2),gp_Dir2d(1,0)),face,
  //  Precision::Confusion());
  //builder.Range(edges[2],face,f1,l1);

  //builder.UpdateEdge(edges[1],
  //  new Geom2d_Line(gp_Pnt2d(l1,0),gp_Dir2d(0,1)),
  //  new Geom2d_Line(gp_Pnt2d(f1,0),gp_Dir2d(0,1)),face,
  //  Precision::Confusion());
  //builder.Range(edges[1],face,f2,l2);

  builder.Add(face,W);
  builder.Add(shell, face);

  shape = shell;
  TopoDS_Face first,last;
  shape = OCCTUtils_MakeShell(shell);

  OCCTUtils_AnalyzeShape(shape);

  return SV_OK;
}

int OCCTUtils_AnalyzeShape(TopoDS_Shape shape)
{
  std::cout<<"Analyzing Shape"<<endl;
  Handle(TopTools_HSequenceOfShape) sl,slv,sle,slw,slf,sls,slo;
  sl = new TopTools_HSequenceOfShape();
  Handle(TColStd_HArray1OfInteger) NbProblems = new
                              TColStd_HArray1OfInteger(1, 36);
  for(int i=1; i<=36; i++)
  {
    NbProblems->SetValue(i,0);
  }

  BRepCheck_Analyzer analyzer(shape, Standard_True);
  TopTools_DataMapOfShapeListOfShape theMap;
  OCCTUtils_GetProblemShapes(analyzer, shape, sl, NbProblems, theMap);

  int isProblem = 0;

  Standard_Integer aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidPointOnCurve);
  if(NbProblems->Value(aProblemID) > 0)
  {
    std::cout<<"  Invalid Point on Curve ................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidPointOnCurveOnSurface);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid Point on CurveOnSurface .......... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidPointOnSurface);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid Point on Surface ................. "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_No3DCurve);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  No 3D Curve .............................. "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_Multiple3DCurve);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Multiple 3D Curve ........................ "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_Invalid3DCurve);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid 3D Curve ......................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_NoCurveOnSurface);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  No Curve on Surface ...................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidCurveOnSurface);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid Curve on Surface ................. "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidCurveOnClosedSurface);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid Curve on closed Surface .......... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidSameRangeFlag);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid SameRange Flag ................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidSameParameterFlag);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid SameParameter Flag ............... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidDegeneratedFlag);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid Degenerated Flag ................. "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_FreeEdge);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Free Edge ................................ "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidMultiConnexity);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid MultiConnexity ................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidRange);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid Range ............................ "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_EmptyWire);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Empty Wire ............................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_RedundantEdge);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Redundant Edge ........................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_SelfIntersectingWire);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Self Intersecting Wire ................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_NoSurface);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  No Surface ............................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidWire);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid Wire ............................. "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_RedundantWire);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Redundant Wire ........................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_IntersectingWires);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Intersecting Wires ....................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidImbricationOfWires);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid Imbrication of Wires ............. "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_EmptyShell);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Empty Shell .............................. "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_RedundantFace);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Redundant Face ........................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_UnorientableShape);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Unorientable Shape ....................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  // Not considering this a problem in case we have open shape
  aProblemID = static_cast<Standard_Integer>(BRepCheck_NotClosed);
  if(NbProblems->Value(aProblemID)>0)
    std::cout<<"  Not Closed ............................... "<<NbProblems->Value(aProblemID)<<"\n";

  aProblemID = static_cast<Standard_Integer>(BRepCheck_NotConnected);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Not Connected ............................ "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_SubshapeNotInShape);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Subshape not in Shape .................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_BadOrientation);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Bad Orientation .......................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_BadOrientationOfSubshape);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Bad Orientation of Subshape .............. "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidToleranceValue);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid tolerance value................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidPolygonOnTriangulation);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid polygon on triangulation.......... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_InvalidImbricationOfShells);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Invalid Imbrication of Shells............. "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

 aProblemID = static_cast<Standard_Integer>(BRepCheck_EnclosedRegion);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  Enclosed Region........................... "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  aProblemID = static_cast<Standard_Integer>(BRepCheck_CheckFail);
  if(NbProblems->Value(aProblemID)>0)
  {
    std::cout<<"  checkshape failure........................ "<<NbProblems->Value(aProblemID)<<"\n";
    isProblem = 1;
  }

  if (isProblem)
  {
    return SV_ERROR;
  }

  return SV_OK;
}

// ---------------------
// OCCTUtils_ShapeFromBSplineSurfaceWithEdges
// ---------------------
int OCCTUtils_ShapeFromBSplineSurfaceWithEdges(const Handle(Geom_BSplineSurface) surface,
    		TopoDS_Shape &shape, std::vector<TopoDS_Edge> edges)
{
  if(surface.IsNull()) {
    fprintf(stderr,"Lofting did not complete\n");
    return SV_ERROR;
  }

  // create the new surface
  TopoDS_Shell shell;
  TopoDS_Face face;
  TopoDS_Wire W;
  TopoDS_Edge edge, couture;
  std::vector<TopoDS_Edge> theseEdges(4);
  std::vector<Handle(Geom_Curve)> theseCurves(4);

  for (int j=0; j<4; j++)
  {
    TopLoc_Location newLoc;
    Standard_Real newFirst, newLast;
    Handle(Geom_Curve) newCurve = BRep_Tool::Curve (edges[j], newLoc, newFirst, newLast);
    theseCurves[j] = newCurve;
  }

  BRep_Builder builder;
  builder.MakeShell(shell);

  TopoDS_Vertex v0f,v0l,v1f,v1l,v2f,v2l,v3f,v3l;

  TopExp::Vertices(edges[0], v0f, v0l);
  TopExp::Vertices(edges[1], v1f, v1l);
  TopExp::Vertices(edges[2], v2f, v2l);
  TopExp::Vertices(edges[3], v3f, v3l);

  // make the face
  builder.MakeFace(face, surface, Precision::Confusion());

  // make the wire
  builder.MakeWire(W);

  // make the missing edges
  Standard_Real f1, f2, l1, l2;
  surface->Bounds(f1,l1,f2,l2);
  fprintf(stdout,"WHAT ARE THESE: %.6f %.6f %.6f %.6f\n", f1, l1, f2, l2);

  // --- edge 1
  //builder.MakeEdge(theseEdges[0], surface->VIso(f2), Precision::Confusion());
  builder.MakeEdge(theseEdges[0], theseCurves[0], Precision::Confusion());
  v0f.Orientation(TopAbs_FORWARD);
  builder.Add(theseEdges[0], v0f);
  v0l.Orientation(TopAbs_REVERSED);
  builder.Add(theseEdges[0], v0l);
  builder.Range(theseEdges[0], f1, l1);
  // processing of looping sections
  // store edges of the 1st section

  // --- edge 2
  //builder.MakeEdge(theseEdges[2], surface->UIso(f1), Precision::Confusion());
  builder.MakeEdge(theseEdges[1], theseCurves[1], Precision::Confusion());
  v1f.Orientation(TopAbs_FORWARD);
  builder.Add(theseEdges[1], v1f);
  v1l.Orientation(TopAbs_REVERSED);
  builder.Add(theseEdges[1], v1l);
  builder.Range(theseEdges[1], f2, l2);
  couture = theseEdges[1];

  // --- edge 4
  theseEdges[3] = couture;
  theseEdges[3].Reverse();

  // --- edge 3
  //builder.MakeEdge(theseEdges[1], surface->VIso(l2), Precision::Confusion());
  builder.MakeEdge(theseEdges[2], theseCurves[2], Precision::Confusion());
  v2f.Orientation(TopAbs_FORWARD);
  builder.Add(theseEdges[2], v2f);
  v2l.Orientation(TopAbs_REVERSED);
  builder.Add(theseEdges[2], v2l);
  builder.Range(theseEdges[2], f1, l1);
  theseEdges[2].Reverse();

  builder.Add(W,theseEdges[0]);
  builder.Add(W,theseEdges[1]);
  builder.Add(W,theseEdges[2]);
  builder.Add(W,theseEdges[3]);

  //// WRITE THE EDGES
  //for (int j=0; j<4; j++)
  //{
  //  TopLoc_Location newLoc;
  //  Standard_Real newFirst, newLast;
  //  Handle(Geom_Curve) newCurve = BRep_Tool::Curve (theseEdges[j], newLoc, newFirst, newLast);

  //  vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
  //  for (int k=0; k<20; k++)
  //  {
  //    double newVal = (newLast-newFirst)*(k/20.) + newFirst;
  //    gp_Pnt newPnt = newCurve->Value(newVal);
  //    newPoints->InsertNextPoint(newPnt.X(), newPnt.Y(), newPnt.Z());
  //  }
  //  vtkSmartPointer<vtkPolyData> newPointsPd = vtkSmartPointer<vtkPolyData>::New();
  //  newPointsPd->SetPoints(newPoints);
  //  std::string newFn = "/Users/adamupdegrove/Desktop/tmp/SODUMBPOINTS_"+std::to_string(j)+".vtp";
  //  vtkSVIOUtils::WriteVTPFile(newFn, newPointsPd);
  //}

  // set PCurve
  builder.UpdateEdge(theseEdges[0],new Geom2d_Line(gp_Pnt2d(0,f2),gp_Dir2d(1,0)),face,
    Precision::Confusion());
  builder.Range(theseEdges[0],face,f1,l1);


  builder.UpdateEdge(theseEdges[2],new Geom2d_Line(gp_Pnt2d(0,l2),gp_Dir2d(1,0)),face,
    Precision::Confusion());
  builder.Range(theseEdges[2],face,f1,l1);

  builder.UpdateEdge(theseEdges[3],
    new Geom2d_Line(gp_Pnt2d(l1,0),gp_Dir2d(0,1)),
    new Geom2d_Line(gp_Pnt2d(f1,0),gp_Dir2d(0,1)),face,
    Precision::Confusion());
  builder.Range(theseEdges[3],face,f2,l2);

  //builder.UpdateEdge(edges[0],new Geom2d_Line(gp_Pnt2d(0,f2),gp_Dir2d(1,0)),face,
  //  Precision::Confusion());
  //builder.Range(edges[0],face,f1,l1);
  //builder.UpdateEdge(edges[2],new Geom2d_Line(gp_Pnt2d(0,l2),gp_Dir2d(1,0)),face,
  //  Precision::Confusion());
  //builder.Range(edges[2],face,f1,l1);

  //builder.UpdateEdge(edges[1],
  //  new Geom2d_Line(gp_Pnt2d(l1,0),gp_Dir2d(0,1)),
  //  new Geom2d_Line(gp_Pnt2d(f1,0),gp_Dir2d(0,1)),face,
  //  Precision::Confusion());
  //builder.Range(edges[1],face,f2,l2);

  builder.Add(face,W);
  builder.Add(shell, face);

  shape = shell;
  TopoDS_Face first,last;
  shape = OCCTUtils_MakeShell(shell);

  OCCTUtils_AnalyzeShape(shape);

  return SV_OK;
}

// ---------------------
// OCCTUtils_ShapeFromBSplineSurface
// ---------------------
int OCCTUtils_ShapeFromBSplineSurface(const Handle(Geom_BSplineSurface) surface,
    		TopoDS_Shape &shape,const TopoDS_Wire &first_wire,
		const TopoDS_Wire &last_wire)
{
  if(surface.IsNull()) {
    fprintf(stderr,"Lofting did not complete\n");
    return SV_ERROR;
  }

  // create the new surface
  TopoDS_Shell shell;
  TopoDS_Face face;
  TopoDS_Wire W;
  TopoDS_Edge edge, edge1, edge2, edge3, edge4, couture;

  BRep_Builder builder;
  builder.MakeShell(shell);

  TopoDS_Wire newW1, newW2;
  BRep_Builder BW1, BW2;
  BW1.MakeWire(newW1);
  BW2.MakeWire(newW2);

  TopLoc_Location loc;
  TopoDS_Vertex v1f,v1l,v2f,v2l;

  Standard_Integer nbPnts = 21;

  TopoDS_Shape firstEdge;

  // segmentation of TS
  Standard_Real Ui1,Ui2,V0,V1;
  Ui1 = 0;
  Ui2 = 1;
  Ui1 = OCCTUtils_PreciseUpar(Ui1, surface);
  Ui2 = OCCTUtils_PreciseUpar(Ui2, surface);
  V0  = surface->VKnot(surface->FirstVKnotIndex());
  V1  = surface->VKnot(surface->LastVKnotIndex());
  surface->Segment(Ui1,Ui2,V0,V1);

  // return vertices
  TopExp_Explorer edge_one(first_wire,TopAbs_EDGE);
  edge =  TopoDS::Edge(edge_one.Current());
  TopExp::Vertices(edge,v1f,v1l);
  if (edge.Orientation() == TopAbs_REVERSED)
    TopExp::Vertices(edge,v1l,v1f);
  firstEdge = edge;

  TopExp_Explorer edge_last(last_wire,TopAbs_EDGE);
  edge =  TopoDS::Edge(edge_last.Current());
  TopExp::Vertices(edge,v2f,v2l);
  if (edge.Orientation() == TopAbs_REVERSED)
    TopExp::Vertices(edge,v2l,v2f);

  // make the face
  builder.MakeFace(face, surface, Precision::Confusion());

  // make the wire
  builder.MakeWire(W);

  // make the missing edges
  Standard_Real f1, f2, l1, l2;
  surface->Bounds(f1,l1,f2,l2);
  fprintf(stdout,"WHAT ARE THESE: %.6f %.6f %.6f %.6f\n", f1, l1, f2, l2);

  // --- edge 1
  builder.MakeEdge(edge1, surface->VIso(f2), Precision::Confusion());
  v1f.Orientation(TopAbs_FORWARD);
  builder.Add(edge1, v1f);
  v1l.Orientation(TopAbs_REVERSED);
  builder.Add(edge1, v1l);
  builder.Range(edge1, f1, l1);
  // processing of looping sections
  // store edges of the 1st section

  // --- edge 2
  builder.MakeEdge(edge2, surface->VIso(l2), Precision::Confusion());
  v2f.Orientation(TopAbs_FORWARD);
  builder.Add(edge2, v2f);
  v2l.Orientation(TopAbs_REVERSED);
  builder.Add(edge2, v2l);
  builder.Range(edge2, f1, l1);
  edge2.Reverse();

  // --- edge 3
  builder.MakeEdge(edge3, surface->UIso(f1), Precision::Confusion());
  v1f.Orientation(TopAbs_FORWARD);
  builder.Add(edge3, v1f);
  v2f.Orientation(TopAbs_REVERSED);
  builder.Add(edge3, v2f);
  builder.Range(edge3, f2, l2);
  couture = edge3;
  edge3.Reverse();

  // --- edge 4
  edge4 = couture;

  builder.Add(W,edge1);
  builder.Add(W,edge4);
  builder.Add(W,edge2);
  builder.Add(W,edge3);

  // set PCurve
  builder.UpdateEdge(edge1,new Geom2d_Line(gp_Pnt2d(0,f2),gp_Dir2d(1,0)),face,
    Precision::Confusion());
  builder.Range(edge1,face,f1,l1);
  builder.UpdateEdge(edge2,new Geom2d_Line(gp_Pnt2d(0,l2),gp_Dir2d(1,0)),face,
    Precision::Confusion());
  builder.Range(edge2,face,f1,l1);

  builder.UpdateEdge(edge3,
    new Geom2d_Line(gp_Pnt2d(l1,0),gp_Dir2d(0,1)),
    new Geom2d_Line(gp_Pnt2d(f1,0),gp_Dir2d(0,1)),face,
    Precision::Confusion());
  builder.Range(edge3,face,f2,l2);

  builder.Add(face,W);
  builder.Add(shell, face);

  //// complete newW1 newW2
  //TopoDS_Edge edge12 = edge1;
  //TopoDS_Edge edge22 = edge2;
  //edge12.Reverse();
  //edge22.Reverse();
  //BW1.Add(newW1, edge12);
  //BW2.Add(newW2, edge22);

  //// history
  //TopTools_DataMapOfShapeShape generated;
  //generated.Bind(firstEdge, face);

  shape = shell;
  TopoDS_Face first,last;
  shape = OCCTUtils_MakeShell(shell);

  return SV_OK;
}

// ---------------------
// OCCTUtils_IsSameOriented
// ---------------------
/**
 * @brief Taken from BRepOffsetAPI_ThruSections
 */
Standard_Boolean OCCTUtils_IsSameOriented(const TopoDS_Shape& aFace,
  const TopoDS_Shape& aShell)
{
  TopExp_Explorer Explo(aFace, TopAbs_EDGE);
  TopoDS_Shape anEdge = Explo.Current();
  TopAbs_Orientation Or1 = anEdge.Orientation();

  TopTools_IndexedDataMapOfShapeListOfShape EFmap;
  TopExp::MapShapesAndAncestors( aShell, TopAbs_EDGE, TopAbs_FACE, EFmap );

  const TopoDS_Shape& AdjacentFace = EFmap.FindFromKey(anEdge).First();
  TopoDS_Shape theEdge;
  for (Explo.Init(AdjacentFace, TopAbs_EDGE); Explo.More(); Explo.Next())
  {
    theEdge = Explo.Current();
    if (theEdge.IsSame(anEdge))
      break;
  }

  TopAbs_Orientation Or2 = theEdge.Orientation();
  if (Or1 == Or2)
    return Standard_False;
  return Standard_True;
}

// ---------------------
// OCCTUtils_IsSameOrientedWEdge
// ---------------------
/**
 * @brief Taken from BRepOffsetAPI_ThruSections
 */
Standard_Boolean OCCTUtils_IsSameOrientedWEdge(const TopoDS_Shape& aFace,
  const TopoDS_Shape& aShell,const TopoDS_Shape &anEdge)
{
  TopAbs_Orientation Or1 = anEdge.Orientation();

  TopTools_IndexedDataMapOfShapeListOfShape EFmap;
  TopExp::MapShapesAndAncestors( aShell, TopAbs_EDGE, TopAbs_FACE, EFmap );

  const TopoDS_Shape& AdjacentFace = EFmap.FindFromKey(anEdge).First();
  TopoDS_Shape theEdge;
  TopExp_Explorer Explo(AdjacentFace, TopAbs_EDGE);
  for (Explo.Init(AdjacentFace, TopAbs_EDGE); Explo.More(); Explo.Next())
  {
    theEdge = Explo.Current();
    if (theEdge.IsSame(anEdge))
      break;
  }

  TopAbs_Orientation Or2 = theEdge.Orientation();
  if (Or1 == Or2)
    return Standard_False;
  return Standard_True;
}

// ---------------------
// OCCTUtils_MakeSolid
// ---------------------
/**
 * @brief Taken from BRepOffsetAPI_ThruSections
 */
TopoDS_Solid OCCTUtils_MakeSolid(TopoDS_Shell& shell, const TopoDS_Wire& wire1,
  const TopoDS_Wire& wire2, const Standard_Real presPln,
  TopoDS_Face& face1, TopoDS_Face& face2)
{
  if (shell.IsNull())
    StdFail_NotDone::Raise("Thrusections is not build");
  Standard_Boolean B = shell.Closed();
  BRep_Builder BB;

  if (!B)
  {
    // It is necessary to close the extremities
    B =  OCCTUtils_PerformPlan(wire1, presPln, face1);
    if (B) {
      B =  OCCTUtils_PerformPlan(wire2, presPln, face2);
      if (B) {
        if (!face1.IsNull() && !OCCTUtils_IsSameOriented( face1, shell ))
          face1.Reverse();
        if (!face2.IsNull() && !OCCTUtils_IsSameOriented( face2, shell ))
          face2.Reverse();

        if (!face1.IsNull())
          BB.Add(shell, face1);
        if (!face2.IsNull())
          BB.Add(shell, face2);

        shell.Closed(Standard_True);
      }
    }
  }

  TopoDS_Solid solid;
  BB.MakeSolid(solid);
  BB.Add(solid, shell);

  // verify the orientation the solid
  BRepClass3d_SolidClassifier clas3d(solid);
  clas3d.PerformInfinitePoint(Precision::Confusion());
  if (clas3d.State() == TopAbs_IN) {
    BB.MakeSolid(solid);
    TopoDS_Shape aLocalShape = shell.Reversed();
    BB.Add(solid, TopoDS::Shell(aLocalShape));
    //    B.Add(solid, TopoDS::Shell(newShell.Reversed()));
  }

  solid.Closed(Standard_True);
  return solid;
}

// ---------------------
// OCCTUtils_PerformPlan
// ---------------------
/**
 * @brief Taken from BRepOffsetAPI_ThruSections
 */
Standard_Boolean OCCTUtils_PerformPlan(const TopoDS_Wire& W,
  const Standard_Real presPln,
  TopoDS_Face& theFace)
{
  Standard_Boolean isDegen = Standard_True;
  TopoDS_Iterator iter(W);
  for (; iter.More(); iter.Next())
  {
    const TopoDS_Edge& anEdge = TopoDS::Edge(iter.Value());
    if (!BRep_Tool::Degenerated(anEdge))
      isDegen = Standard_False;
  }
  if (isDegen)
    return Standard_True;

  Standard_Boolean Ok = Standard_False;
  if (!W.IsNull()) {
    BRepBuilderAPI_FindPlane Searcher( W, presPln );
    if (Searcher.Found())
    {
      theFace = BRepBuilderAPI_MakeFace(Searcher.Plane(), W);
      Ok = Standard_True;
    }
    else // try to find another surface
    {
      BRepBuilderAPI_MakeFace MF( W );
      if (MF.IsDone())
      {
        theFace = MF.Face();
        Ok = Standard_True;
      }
    }
  }

  return Ok;
}

// ---------------------
// OCCTUtils_MakeShell
// ---------------------
/**
 * @brief Taken from BRepOffsetAPI_ThruSections
 */
TopoDS_Solid OCCTUtils_MakeShell(TopoDS_Shell& shell)
{
  if (shell.IsNull())
    StdFail_NotDone::Raise("Thrusections is not build");
  Standard_Boolean B = shell.Closed();
  BRep_Builder BB;

  //if (!B)
  //{
  //  // It is necessary to close the extremities
  //  B =  OCCTUtils_PerformPlan(wire1, presPln, face1);
  //  if (B) {
  //    B =  OCCTUtils_PerformPlan(wire2, presPln, face2);
  //    if (B) {
  //      if (!face1.IsNull() && !OCCTUtils_IsSameOriented( face1, shell ))
  //        face1.Reverse();
  //      if (!face2.IsNull() && !OCCTUtils_IsSameOriented( face2, shell ))
  //        face2.Reverse();

  //      if (!face1.IsNull())
  //        BB.Add(shell, face1);
  //      if (!face2.IsNull())
  //        BB.Add(shell, face2);

  //      shell.Closed(Standard_True);
  //    }
  //  }
  //}

  TopoDS_Solid orientedshell;
  BB.MakeSolid(orientedshell);
  BB.Add(orientedshell, shell);

  // verify the orientation the solid
  BRepClass3d_SolidClassifier clas3d(orientedshell);
  clas3d.PerformInfinitePoint(Precision::Confusion());
  if (clas3d.State() == TopAbs_IN) {
    BB.MakeSolid(orientedshell);
    TopoDS_Shape aLocalShape = shell.Reversed();
    BB.Add(orientedshell, TopoDS::Shell(aLocalShape));
    //    B.Add(orientedshell, TopoDS::Shell(newShell.Reversed()));
  }

  //orientedshell.Closed(Standard_True);
  return orientedshell;
}

// ---------------------
// OCCTUtils_PreciseUpar
// ---------------------
/**
 * @brief Taken from BRepOffsetAPI_ThruSections
 */
Standard_Real OCCTUtils_PreciseUpar(const Standard_Real anUpar,
  const Handle(Geom_BSplineSurface)& aSurface)
{
  Standard_Real Tol = Precision::PConfusion();
  Standard_Integer i1, i2;

  aSurface->LocateU(anUpar, Tol, i1, i2);
  Standard_Real U1 = aSurface->UKnot(i1);
  Standard_Real U2 = aSurface->UKnot(i2);

  Standard_Real NewU = anUpar;

  NewU = (anUpar - U1 < U2 - anUpar)? U1 : U2;
  return NewU;
}

// ---------------------
// OCCTUtils_EdgeToBSpline
// ---------------------
/**
 * @brief Taken from BRepOffsetAPI_ThruSections
 */
Handle(Geom_BSplineCurve) OCCTUtils_EdgeToBSpline(const TopoDS_Edge& theEdge)
{
  Handle(Geom_BSplineCurve) aBSCurve;
  if (BRep_Tool::Degenerated(theEdge)) {
    // degenerated edge : construction of a point curve
    TColStd_Array1OfReal aKnots (1,2);
    aKnots(1) = 0.;
    aKnots(2) = 1.;

    TColStd_Array1OfInteger aMults (1,2);
    aMults(1) = 2;
    aMults(2) = 2;

    TColgp_Array1OfPnt aPoles(1,2);
    TopoDS_Vertex vf, vl;
    TopExp::Vertices(theEdge,vl,vf);
    aPoles(1) = BRep_Tool::Pnt(vf);
    aPoles(2) = BRep_Tool::Pnt(vl);

    aBSCurve = new Geom_BSplineCurve (aPoles, aKnots, aMults, 1);
  }
  else
  {
    // get the curve of the edge
    TopLoc_Location aLoc;
    Standard_Real aFirst, aLast;
    Handle(Geom_Curve) aCurve = BRep_Tool::Curve (theEdge, aLoc, aFirst, aLast);
    if (aCurve.IsNull())
      Standard_NullObject::Raise("Null 3D curve in edge");

    // convert its part used by edge to bspline; note that if edge curve is bspline,
    // conversion made via trimmed curve is still needed -- it will copy it, segment
    // as appropriate, and remove periodicity if it is periodic (deadly for approximator)
    Handle(Geom_TrimmedCurve) aTrimCurve = new Geom_TrimmedCurve (aCurve, aFirst, aLast);

    // special treatment of conic curve
    if (aTrimCurve->BasisCurve()->IsKind(STANDARD_TYPE(Geom_Conic)))
    {
      const Handle(Geom_Curve)& aCurveTrimmed = aTrimCurve; // to avoid ambiguity
      GeomConvert_ApproxCurve anAppr (aCurveTrimmed, Precision::Confusion(), GeomAbs_C1, 16, 14);
      if (anAppr.HasResult())
        aBSCurve = anAppr.Curve();
    }

    // general case
    if (aBSCurve.IsNull())
      aBSCurve = GeomConvert::CurveToBSplineCurve (aTrimCurve);

    // apply transformation if needed
    if (! aLoc.IsIdentity())
      aBSCurve->Transform (aLoc.Transformation());

    // reparameterize to [0,1]
    TColStd_Array1OfReal aKnots (1, aBSCurve->NbKnots());
    aBSCurve->Knots (aKnots);
    BSplCLib::Reparametrize (0., 1., aKnots);
    aBSCurve->SetKnots (aKnots);
  }

  // reverse curve if edge is reversed
  if (theEdge.Orientation() == TopAbs_REVERSED)
    aBSCurve->Reverse();

  return aBSCurve;
}

// ---------------------
// OCCTUtils_GetFaceIds
// ---------------------
/**
 * @brief Procedure to get face numbers that correspond to the scalars
 * assigned to the geometry
 * @param *geom input TopoDS_Shape on which to get the face ids
 * @param *v_num_faces int that contains the number of total face regions
 * @param **v_faces vector containing the array of numerical values
 * corresponding to each face region
 * @return SV_OK if function completes properly
 */
int OCCTUtils_GetFaceIds( const TopoDS_Shape &geom,
		Handle(XCAFDoc_ShapeTool) &shapetool,TDF_Label &shapelabel,
	       	int *v_num_faces, int **v_faces)
{
  int num = 0;

  const TopoDS_Shape& aShape = geom;

  TopExp_Explorer anExp (aShape, TopAbs_FACE);
  for (; anExp.More(); anExp.Next()) {
   const TopoDS_Face& aFace = TopoDS::Face (anExp.Current());
   num++;
  }

  *v_num_faces = num;

  if (num == 0) return SV_ERROR;

  (*v_faces) = new int [num];

  TopExp_Explorer anExp2 (aShape, TopAbs_FACE);

  int j = 0;

  for (; anExp2.More(); anExp2.Next()) {
    TopoDS_Face aFace = TopoDS::Face (anExp2.Current());
    int faceId= -1;
    OCCTUtils_GetFaceLabel(aFace,shapetool,shapelabel,faceId);
    //(*v_faces)[j] = aFace.HashCode(9999999999);
    (*v_faces)[j] = faceId;
    j++;
  }

  return SV_OK;
}

// -----------------
// OCCTUtils_GetFaceLabel
// -----------------
/**
 * @brief Procedure to get a singular face number
 * @param *geom input TopoDS_Shape on which to get the face ids
 * @param id is the scalar of the face region
 * @param returns input value and an error if the face does not have id
 * @return SV_OK if function completes properly
 */
int OCCTUtils_GetFaceLabel(const TopoDS_Shape &geom,
		const Handle(XCAFDoc_ShapeTool) &shapetool,TDF_Label &shapelabel,
	       	int &id)
{
  TDF_Label tmpLabel;
  shapetool->FindSubShape(shapelabel,geom,tmpLabel);
  if (tmpLabel.IsNull())
  {
    //fprintf(stderr,"Face does not have label\n");
    return SV_ERROR;
  }

  TDF_Label idLabel = tmpLabel.FindChild(0);
  if (idLabel.IsNull())
  {
    fprintf(stderr,"Face does not have id\n");
    return SV_ERROR;
  }
  //Retrive attribute
  Handle(TDataStd_Integer) INT = new TDataStd_Integer();
  idLabel.FindAttribute(TDataStd_Integer::GetID(),INT);
  id = INT->Get();

  return SV_OK;
}

// -------------------
// OCCTUtils_RenumberFaces
// -------------------
/**
 * @brief Procedure to renumber faces (not used)
 * @return SV_OK if function completes properly
 */
int OCCTUtils_RenumberFaces(TopoDS_Shape &shape,
		Handle(XCAFDoc_ShapeTool) &shapetool,TDF_Label &shapelabel)
{
  int *faces;
  int numFaces;
  int facerange;
  OCCTUtils_GetFaceIds(shape,shapetool,shapelabel,&numFaces,&faces);
  OCCTUtils_GetFaceRange(shape,shapetool,shapelabel,facerange);

  //Initialize to zero
  int *newmap = new int[numFaces];

  //Increase map one at spot for each face id
  for (int i=0;i<numFaces;i++)
    newmap[i] = -2;

  int checkid=-1;
  int currentid=1;
  TopExp_Explorer anExp(shape,TopAbs_FACE);
  while (checkid != facerange)
  {
    anExp.Init(shape,TopAbs_FACE);
    int found=0;
    for (int j=0;anExp.More();anExp.Next(),j++)
    {
      TopoDS_Face tmpFace = TopoDS::Face(anExp.Current());
      int faceid=-1;
      OCCTUtils_GetFaceLabel(tmpFace,shapetool,shapelabel,faceid);
      if (faceid == checkid)
      {
	if (newmap[j] == -2)
	{
	  newmap[j] = currentid++;
	  found=1;
	  break;
	}
      }
    }
    if (!found)
      checkid++;
  }

  anExp.Init(shape,TopAbs_FACE);
  for (int i=0;anExp.More();anExp.Next(),i++)
  {
    TopoDS_Face tmpFace = TopoDS::Face(anExp.Current());
    if (OCCTUtils_ReLabelFace(tmpFace,shapetool,shapelabel,newmap[i]) != SV_OK)
    {
      fprintf(stderr,"Could not label face\n");
      return SV_ERROR;
    }
  }

  delete [] newmap;
  delete [] faces;
  return SV_OK;
}

// ----------------
// GetFaceRange
// ----------------
/**
 * @brief Procedure to get the range of id values on a shape
 * @param face_range the returned range of faces
 * @return SV_OK if function completes properly
 */
int OCCTUtils_GetFaceRange(const TopoDS_Shape &shape,
		Handle(XCAFDoc_ShapeTool) &shapetool,TDF_Label &shapelabel,
    		int &face_range)
{
  face_range = 0;
  TopExp_Explorer anExp(shape,TopAbs_FACE);
  for (int i=0;anExp.More();anExp.Next())
  {
    const TopoDS_Face &tmpFace = TopoDS::Face(anExp.Current());
    int faceid =-1;
    OCCTUtils_GetFaceLabel(tmpFace,shapetool,shapelabel,faceid);
    if (faceid > face_range)
      face_range = faceid;
  }

  return SV_OK;
}


// -------------------
// OCCTUtils_GetOrientation
// -------------------
/**
 * @brief Procedure to get a shape orientation
 * @param *geom input TopoDS_Shape on which to get orientation
 * @return SV_OK if function completes properly
 */
int OCCTUtils_GetOrientation(const TopoDS_Shape &shape,int &orientation)
{
  orientation = (int) shape.Orientation();
  return SV_OK;
}

// -------------------
// OCCTUtils_SetOrientation
// -------------------
/**
 * @brief Procedure to set a face orientation
 * @param shape input TopoDS_Shape on which to set a faces orientation
 * @param face the face on which to set the orientation
 * @return SV_OK if function completes properly
 */
int OCCTUtils_SetOrientation(TopoDS_Shape &shape,TopoDS_Shape &face,int &orientation)
{
  Handle(BRepTools_ReShape) reshaper =  new BRepTools_ReShape();
  reshaper->ModeConsiderOrientation() = Standard_True;

  TopoDS_Shape compFace = face.Complemented();
  reshaper->Replace(face,compFace,Standard_True);
  TopoDS_Shape tmpShape = reshaper->Apply(shape,TopAbs_FACE);
  shape = tmpShape;

  return SV_OK;
}

// -------------------
// OCCTUtils_ReLabelFace
// -------------------
/**
 * @brief Procedure to relabel the face of a shape
 * @param shape input TopoDS_Shape face that needs to be relabled
 * @param id desired id for face
 * @return SV_OK if function completes properly
 */
int OCCTUtils_ReLabelFace( TopoDS_Shape &shape,
		Handle(XCAFDoc_ShapeTool) &shapetool,TDF_Label &shapelabel,
		int &id)
{
  if (shape.IsNull())
  {
    fprintf(stderr,"Face is NULL, cannot add\n");
    return SV_ERROR;
  }

  TDF_Label tmpLabel;
  shapetool->FindSubShape(shapelabel,shape,tmpLabel);
  if (tmpLabel.IsNull())
  {
    fprintf(stderr,"Face has not been given a label\n");
    return SV_ERROR;
  }
  TDF_Label idLabel = tmpLabel.FindChild(0);
  if (idLabel.IsNull())
  {
    fprintf(stderr,"Face has not been given an id\n");
    return SV_ERROR;
  }
  Handle(TDataStd_Integer) INT = new TDataStd_Integer();
  idLabel.FindAttribute(TDataStd_Integer::GetID(),INT);
  INT->Set(id);

  return SV_OK;
}

// -------------------
// OCCTUtils_GetNumberFaces
// -------------------
/**
 * @brief Procedure to return the number of faces in shape
 * @param shape input TopoDS_Shape to find number of faces
 * @param num_faces returns the number of faces found
 * @return SV_OK if function completes properly
 */
int OCCTUtils_GetNumberOfFaces(const TopoDS_Shape &shape,int &num_faces)
{
  num_faces = 0;
  TopExp_Explorer anExp(shape,TopAbs_FACE);
  for (int i=0;anExp.More();anExp.Next())
    num_faces++;

  return SV_OK;
}

// -------------------
// OCCTUtils_GetFaceAttribute
// -------------------
/**
 * @brief Procedure to get an attribute of the shape
 * @param shape input TopoDS_Shape to get attribute
 * @param shapetool the XDEDoc manager that contains attribute info
 * @param shapelabel the label for the shape registered in XDEDoc
 * @note attributes includ id, gdscName, and parent
 * @return SV_OK if function completes properly
 */
int OCCTUtils_GetFaceAttribute(const TopoDS_Shape &face,
		Handle(XCAFDoc_ShapeTool) &shapetool,TDF_Label &shapelabel,
    				char *attr, char **value)
{
  static char returnString[255];
  TDF_Label tmpLabel;
  shapetool->FindSubShape(shapelabel,face,tmpLabel);
  if (tmpLabel.IsNull())
  {
    fprintf(stderr,"Face is not labelled and thus has no attribute\n");
    return SV_ERROR;
  }
  if (!strncmp(attr,"gdscName",4))
  {
    TDF_Label nameLabel = tmpLabel.FindChild(1,Standard_False);
    if (nameLabel.IsNull())
    {
      fprintf(stderr,"gdscName label doesn't exist, cannot retrive name\n");
      return SV_ERROR;
    }
    Handle(TDataStd_ExtStringArray) NSTRING = new
      TDataStd_ExtStringArray();
    int isLabel = nameLabel.FindAttribute(TDataStd_ExtStringArray::GetID(),NSTRING);
    if (isLabel == 0)
    {
      fprintf(stderr,"gdscName attribute does not exist on face\n");
      return SV_ERROR;
    }
    returnString[0]='\0';
    OCCTUtils_GetExtStringArrayAsChar(NSTRING,returnString);
    *value = returnString;
  }
  else if (!strncmp(attr,"parent",6))
  {
    TDF_Label parentLabel = tmpLabel.FindChild(2,Standard_False);
    if (parentLabel.IsNull())
    {
      fprintf(stderr,"gdscName label doesn't exist, cannot retrive name\n");
      return SV_ERROR;
    }
    Handle(TDataStd_ExtStringArray) PSTRING = new
      TDataStd_ExtStringArray();
    int isLabel = parentLabel.FindAttribute(TDataStd_ExtStringArray::GetID(),PSTRING);
    if (isLabel == 0)
    {
      fprintf(stderr,"parent attribute does not exist on face\n");
      return SV_ERROR;
    }
    returnString[0]='\0';
    OCCTUtils_GetExtStringArrayAsChar(PSTRING,returnString);
    *value = returnString;
  }
  else if (!strncmp(attr,"id",2))
  {
    TDF_Label idLabel = tmpLabel.FindChild(0);
    Handle(TDataStd_Integer) INT = new TDataStd_Integer();
    int isLabel = idLabel.FindAttribute(TDataStd_Integer::GetID(),INT);
    if (isLabel == 0)
    {
      fprintf(stderr,"id attribute does not exist on face\n");
      return SV_ERROR;
    }
    fprintf(stderr,"Inside id and want to check the actual id!! %d\n",INT->Get());
    returnString[0]='\0';
    sprintf(returnString,"%d",INT->Get());
    *value = returnString;
  }
  else
  {
    fprintf(stderr,"Attribute %s is not attribute of shape. Options are gdscName, parent, id\n",attr);
    return SV_ERROR;
  }

  return SV_OK;
}

// -------------------
// OCCTUtils_SetFaceAttribute
// -------------------
/**
 * @brief Procedure to set an attribute of the shape
 * @param shape input TopoDS_Shape to set attribute
 * @param shapetool the XDEDoc manager that contains attribute info
 * @param shapelabel the label for the shape registered in XDEDoc
 * @note attributes include id, name, and parent
 * @return SV_OK if function completes properly
 */
int OCCTUtils_SetFaceAttribute(const TopoDS_Shape &face,
		Handle(XCAFDoc_ShapeTool) &shapetool,TDF_Label &shapelabel,
    				char *attr, char *value)
{
  TDF_Label tmpLabel;
  shapetool->FindSubShape(shapelabel,face,tmpLabel);
  if (tmpLabel.IsNull())
  {
    fprintf(stderr,"Face is not labelled and thus has no attribute\n");
    return SV_ERROR;
  }
  if (!strncmp(attr,"gdscName",4))
  {
    TDF_Label nameLabel = tmpLabel.FindChild(1,Standard_False);
    Handle(TDataStd_ExtStringArray) NSTRING = new
      TDataStd_ExtStringArray();
    int isLabel = nameLabel.FindAttribute(TDataStd_ExtStringArray::GetID(),NSTRING);

    OCCTUtils_SetExtStringArrayFromChar(NSTRING,value);
  }
  else if (!strncmp(attr,"parent",6))
  {
    TDF_Label parentLabel = tmpLabel.FindChild(2,Standard_False);
    Handle(TDataStd_ExtStringArray) PSTRING = new
      TDataStd_ExtStringArray();
    int isLabel = parentLabel.FindAttribute(TDataStd_ExtStringArray::GetID(),PSTRING);

    OCCTUtils_SetExtStringArrayFromChar(PSTRING,value);
  }
  else if (!strncmp(attr,"id",2))
  {
    TDF_Label idLabel = tmpLabel.FindChild(0,Standard_False);
    Handle(TDataStd_Integer) INT = new TDataStd_Integer();
    int isLabel = idLabel.FindAttribute(TDataStd_Integer::GetID(),INT);

    INT->Set(charToInt(value));
  }
  else
  {
    fprintf(stderr,"Attribute %s is not attribute of shape. Options are gdscName, parent, id\n",attr);
    return SV_ERROR;
  }

  return SV_OK;
}

// -------------------
// OCCTUtils_PassFaceAttributes
// -------------------
/**
 * @brief Procedure to get an attribute of the shape
 * @param shape input TopoDS_Shape to get attribute
 * @param shapetool the XDEDoc manager that contains attribute info
 * @param shapelabel the label for the shape registered in XDEDoc
 * @note attributes includ id, name, and parent
 * @return SV_OK if function completes properly
 */
int OCCTUtils_PassFaceAttributes(TopoDS_Shape &faceSrc,TopoDS_Shape &faceDst,
		Handle(XCAFDoc_ShapeTool) &shapetool,TDF_Label &labelSrc,
		TDF_Label &labelDst)
{
  //Pass the face names first
  char *name;
  if (OCCTUtils_GetFaceAttribute(
	faceSrc,shapetool,labelSrc,"gdscName",&name) != SV_OK)
  {
    fprintf(stderr,"Failure in getting gdscName for shape\n");
    return SV_ERROR;
  }
  if (OCCTUtils_SetFaceAttribute(
	faceDst,shapetool,labelDst,"gdscName",name) != SV_OK)
  {
    fprintf(stderr,"Failure in setting gdscName for shape\n");
    return SV_ERROR;
  }

  //Now parent name
  char *parent;
  if (OCCTUtils_GetFaceAttribute(
	faceSrc,shapetool,labelSrc,"parent",&parent) != SV_OK)
  {
    fprintf(stderr,"Failure in getting parent for shape\n");
    return SV_ERROR;
  }
  if (OCCTUtils_SetFaceAttribute(
	faceDst,shapetool,labelDst,"parent",parent) != SV_OK)
  {
    fprintf(stderr,"Failure in setting parent for shape\n");
    return SV_ERROR;
  }
  return SV_OK;
}

// -------------------
// OCCTUtils_CheckIsSolid
// -------------------
/**
 * @brief Procedure to check and see if it is solid
 * @param shape input TopoDS_Shape to check
 * @param issue contains integer for issue. BRepCheck_Status doc of occt
 * @return SV_OK if solid, SV_ERROR if it is not solid
 */
int OCCTUtils_CheckIsSolid(const TopoDS_Shape &shape,int &issue)
{
  try {
    BRepCheck_Solid solidchecker(TopoDS::Solid(shape));
    BRepCheck_ListOfStatus status = solidchecker.Status();
    BRepCheck_ListIteratorOfListOfStatus statit;
    statit.Initialize(status);
    for (int i=0;statit.More();statit.Next())
    {
      BRepCheck_Status checker = statit.Value();
      if (checker != 0)
      {
        issue = checker;
        fprintf(stderr,"Shape is not solid!\n");
        return SV_ERROR;
      }
    }
    issue = 0;
  }
  catch (Standard_TypeMismatch) {
    fprintf(stderr,"Shape caused error in solid check\n");
    return SV_ERROR;
  }

  return SV_OK;
}

// ---------------------
// OCCTUtils_SewShapes
// ---------------------
/**
 * @brief
 * @param
 * @return SV_OK if function completes properly
 */
int OCCTUtils_SewShapes(std::vector<TopoDS_Shape> shapeList, double tolerance, TopoDS_Shape &newShape)
{
  for (int i=0; i<shapeList.size(); i++)
  {
    if (shapeList[i].Closed())
    {
      fprintf(stdout,"Shape %d is closed, shouldn't be\n", i);
      return SV_OK;
    }
  }

  //Attacher!
  BRepBuilderAPI_Sewing attacher(tolerance);
  for (int i=0; i<shapeList.size(); i++)
  {
    attacher.Add(shapeList[i]);
  }
  attacher.Perform();

  TopoDS_Shell tmpShell;
  try {
    tmpShell = TopoDS::Shell(attacher.SewedShape());
  }
  catch (Standard_TypeMismatch) {
    fprintf(stderr,"Couldn't sew surfaces together\n");
    return SV_ERROR;
  }

  BRepBuilderAPI_MakeSolid solidmaker(tmpShell);
  newShape = solidmaker.Solid();

  return SV_OK;
}

void OCCTUtils_GetProblemShapes(const BRepCheck_Analyzer& Ana,
                      const TopoDS_Shape& Shape,
                      Handle(TopTools_HSequenceOfShape)& sl,
                      Handle(TColStd_HArray1OfInteger)& NbProblems,
                      TopTools_DataMapOfShapeListOfShape &theMap
                      )
{
  for (TopoDS_Iterator iter(Shape); iter.More(); iter.Next())
  {
    OCCTUtils_GetProblemShapes(Ana,iter.Value(),sl, NbProblems, theMap);
  }
  TopAbs_ShapeEnum styp = Shape.ShapeType();
  BRepCheck_ListIteratorOfListOfStatus itl;
  if (!Ana.Result(Shape).IsNull() && !theMap.IsBound(Shape))
  {
    itl.Initialize(Ana.Result(Shape)->Status());

    if (itl.Value() != BRepCheck_NoError)
    {
      sl->Append(Shape);
      OCCTUtils_FillProblems(itl.Value(),NbProblems);
    }
  }
  if (!theMap.IsBound(Shape))
  {
    TopTools_ListOfShape thelist;
    theMap.Bind(Shape, thelist);
  }

  switch (styp) {
  case TopAbs_EDGE:
    OCCTUtils_GetProblemSub(Ana, Shape, sl, NbProblems, theMap, TopAbs_VERTEX);
    break;
  case TopAbs_FACE:
    OCCTUtils_GetProblemSub(Ana, Shape, sl, NbProblems, theMap, TopAbs_WIRE);
    OCCTUtils_GetProblemSub(Ana, Shape, sl, NbProblems, theMap, TopAbs_EDGE);
    OCCTUtils_GetProblemSub(Ana, Shape, sl, NbProblems, theMap, TopAbs_VERTEX);
    break;
  case TopAbs_SHELL:
    break;
  case TopAbs_SOLID:
    OCCTUtils_GetProblemSub(Ana, Shape, sl, NbProblems, theMap, TopAbs_SHELL);
    break;
  default:
    break;
  }
}

void OCCTUtils_GetProblemSub(const BRepCheck_Analyzer& Ana,
                   const TopoDS_Shape& Shape,
                   Handle(TopTools_HSequenceOfShape)& sl,
                   Handle(TColStd_HArray1OfInteger)& NbProblems,
                   TopTools_DataMapOfShapeListOfShape &theMap,
                   const TopAbs_ShapeEnum Subtype)
{
  BRepCheck_ListIteratorOfListOfStatus itl;
  TopExp_Explorer exp;
  for (exp.Init(Shape,Subtype); exp.More(); exp.Next()) {
    const Handle(BRepCheck_Result)& res = Ana.Result(exp.Current());

    const TopoDS_Shape& sub = exp.Current();
    for (res->InitContextIterator(); res->MoreShapeInContext(); res->NextShapeInContext())
    {
      if (res->ContextualShape().IsSame(Shape) && !OCCTUtils_Contains(theMap(sub),Shape))
      {
        theMap(sub).Append(Shape);

        itl.Initialize(res->StatusOnShape());

        if (itl.Value() != BRepCheck_NoError)
        {
          Standard_Integer ii = 0;

          for(ii=1; ii<=sl->Length(); ii++)
          {
            if(sl->Value(ii).IsSame(sub))
            {
              break;
            }
          }

          if(ii>sl->Length())
          {
            sl->Append(sub);
            OCCTUtils_FillProblems(itl.Value(),NbProblems);
          }
          for(ii=1; ii<=sl->Length(); ii++)
          {
            if(sl->Value(ii).IsSame(Shape))
            {
              break;
            }
          }
          if(ii>sl->Length())
          {
            sl->Append(Shape);
            OCCTUtils_FillProblems(itl.Value(),NbProblems);
          }
        }
        break;
      }
    }
  }
}

void OCCTUtils_FillProblems(const BRepCheck_Status stat,
                            Handle(TColStd_HArray1OfInteger)& NbProblems)
{

  const Standard_Integer anID = static_cast<Standard_Integer> (stat);

  if((NbProblems->Upper() < anID) || (NbProblems->Lower() > anID))
    return;

  NbProblems->SetValue(anID, NbProblems->Value(anID)+1);

}

Standard_Boolean OCCTUtils_Contains(const TopTools_ListOfShape& L,
				 const TopoDS_Shape& S)
{
  TopTools_ListIteratorOfListOfShape it;
  for (it.Initialize(L); it.More(); it.Next()) {
    if (it.Value().IsSame(S)) {
      return Standard_True;
    }
  }
  return Standard_False;
}

TopoDS_Shape OCCTUtils_GetFirstType(const TopoDS_Shape &shape, TopAbs_ShapeEnum type)
{
  if (!shape.IsNull())
  {
      TopExp_Explorer it;
      for (it.Init(shape, type); it.More(); it.Next())
          return it.Current();
  }
  throw std::runtime_error("Couldn't find first object of type");
}

// ---------------------
// OCCTUtils_SplitEdgeOnShapeNearPoints
// ---------------------
/**
 * @brief
 * @param
 * @return SV_OK if function completes properly
 */
int OCCTUtils_SplitEdgeOnShapeNearPoints(TopoDS_Shape &shape, const int edgeNumber, vtkPoints *points)
{
  int numSplitPoints = points->GetNumberOfPoints();

  if (numSplitPoints == 0)
  {
    std::cout << "No points given" << endl;
    return SV_ERROR;
  }

  TopoDS_Edge edge;
  if (OCCTUtils_GetOpenEdge(shape, edgeNumber, edge) != SV_OK)
  {
    std::cout << "Failed getting open edge on shape" << endl;
    return SV_ERROR;
  }

  try
  {
    // Get curve from edge
    Standard_Real pFirst, pLast;
    Handle(Geom_Curve) curve = BRep_Tool::Curve (edge, pFirst, pLast);

    // Loop through points and get all param split locations
    double pt[3];
    std::vector<double> splitParams(numSplitPoints);
    for (int i=0; i<numSplitPoints; i++)
    {
      // Get vertex
      points->GetPoint(i, pt);
      gp_Pnt occtPnt(pt[0], pt[1], pt[2]);
      TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(occtPnt);

      // Get closest point on curve
      BRepExtrema_DistShapeShape closestPointFinder(edge, vertex);
      closestPointFinder.Perform();

      //fprintf(stdout,"ACTUAL CLOSE POINT 0: %.6f %.6f %.6f\n", closestPointFinder.PointOnShape1(1).X(), closestPointFinder.PointOnShape1(1).Y(), closestPointFinder.PointOnShape1(1).Z());

      Standard_Real newParam;
      closestPointFinder.ParOnEdgeS1(1, newParam);
      splitParams[i] = newParam;
      //fprintf(stdout,"ACTUAL CLOSE PARAMETER 0: %.6f\n", newParam);
    }

    // Order the params so that we split in order
    std::sort(splitParams.begin(), splitParams.end());

    int numFullPoints = numSplitPoints+2;
    std::vector<double> fullSplitParams(numFullPoints);
    fullSplitParams[0] = pFirst;
    for (int i=0; i<numSplitPoints; i++)
    {
      fullSplitParams[i+1] = splitParams[i];
    }
    fullSplitParams[numFullPoints-1] = pLast;

    // Start and end vertices
    TopoDS_Vertex vStart, vEnd;
    vStart = TopExp::FirstVertex(edge);
    vEnd = TopExp::LastVertex(edge);

    // Construct new split edge with builder
    BRep_Builder builder;

    // Get  vertices
    std::vector<TopoDS_Vertex> vParams(numFullPoints);
    vParams[0] = vStart;
    for (int i=0; i<splitParams.size(); i++)
    {
      gp_Pnt tmpPnt = curve->Value(splitParams[i]);
      builder.MakeVertex(vParams[i+1], tmpPnt, Precision::Confusion());
    }
    vParams[numFullPoints-1] = vEnd;

    // Stitch together all edges
    TopoDS_Wire newWire;
    builder.MakeWire(newWire);

    // create new edges
    int numSplitEdges = numSplitPoints+1;
    std::vector<TopoDS_Edge> newEdges(numSplitEdges);

    // build up new edges
    for (int i=0; i<numSplitEdges; i++)
    {
      builder.MakeEdge(newEdges[i], curve, Precision::Confusion());
      vParams[i].Orientation(TopAbs_FORWARD);
      builder.Add(newEdges[i], vParams[i]);
      vParams[i+1].Orientation(TopAbs_REVERSED);
      builder.Add(newEdges[i], vParams[i+1]);
      builder.Range(newEdges[i], fullSplitParams[i], fullSplitParams[i+1]);
      newEdges[i].Orientation(edge.Orientation());
    }

    // add edges
    for (int i=0; i<numSplitEdges; i++)
    {
      builder.Add(newWire, newEdges[i]);
    }

    fprintf(stdout,"CHECKING WIRE\n");
    if (OCCTUtils_AnalyzeShape(newWire) != SV_OK)
    {
      std::cerr << "Error forming new wire" << endl;
      return SV_ERROR;
    }

    // Assuming there is only one face. If there are multiple, this needs to be changed
    TopoDS_Face face = TopoDS::Face(OCCTUtils_GetFirstType(shape, TopAbs_FACE));

    // UV points of edge
    gp_Pnt2d pntFirst, pntLast;
    BRep_Tool::UVPoints(edge, face, pntFirst, pntLast);

    // Add pcurves for all edges
    // Assuming that the edge wraps around in U, need to modify if edge goes the other way
    for (int i=0; i<numSplitEdges; i++)
    {
      builder.UpdateEdge(newEdges[i], new Geom2d_Line(gp_Pnt2d(pntFirst.X(), pntFirst.Y()), gp_Dir2d(1,0)), face, Precision::Confusion());
      builder.Range(newEdges[i], face, fullSplitParams[i], fullSplitParams[i+1]);
    }

    // Now replace edge with split edge
    Handle(ShapeBuild_ReShape) context = new ShapeBuild_ReShape();
    context->Replace(edge, newWire);
    TopoDS_Shape newShape = context->Apply(shape);

    if (OCCTUtils_AnalyzeShape(newShape) != SV_OK)
    {
      std::cerr << "Issue with shape after edge split" << endl;
      return SV_ERROR;
    }

    //TopExp_Explorer edgeExp;
    //edgeExp.Init(*(loftedSurfs[surfer]->geom_), TopAbs_EDGE);

    //for (int r=0; edgeExp.More(); edgeExp.Next(), r++)
    //{
    //  TopoDS_Edge thisEdge = TopoDS::Edge(edgeExp.Current());

    //  Standard_Real paraFirst, paraLast;
    //  Handle(Geom_Curve) thisCurve = BRep_Tool::Curve (thisEdge, paraFirst, paraLast);
    //  Standard_Real faceFirst, faceLast;
    //  Handle(Geom2d_Curve) aPCurve = BRep_Tool::CurveOnSurface (thisEdge, face, faceFirst, faceLast);

    //  fprintf(stdout,"EDGE %d HAS ORIEN %d\n", r, thisEdge.Orientation());
    //  fprintf(stdout,"EDGE %d HAS RANGE %.6f %.6f\n", r, paraFirst, paraLast);
    //  fprintf(stdout,"EDGE %d HAS RANGE ON FACE %.6f %.6f\n", r, faceFirst, faceLast);

    //  gp_Pnt2d pntFirst, pntLast;
    //  BRep_Tool::UVPoints(thisEdge, face, pntFirst, pntLast);

    //  fprintf(stdout,"EDGE %d HAS STARTING UV VALUE: %.6f %.6f\n", r, pntFirst.X(), pntFirst.Y());
    //  fprintf(stdout,"EDGE %d HAS ENDERING UV VALUE: %.6f %.6f\n", r, pntLast.X(), pntLast.Y());


    //  vtkNew(vtkPoints, writePoints);
    //  for (int s=0; s<20; s++)
    //  {
    //    Standard_Real pointVal = paraFirst + (s/20.)*(paraLast-paraFirst);
    //    gp_Pnt writePoint = thisCurve->Value(pointVal);

    //    writePoints->InsertNextPoint(writePoint.X(), writePoint.Y(), writePoint.Z());
    //  }

    //  vtkNew(vtkPolyData, writePointsPd);
    //  writePointsPd->SetPoints(writePoints);

    //  std::string fn = "/Users/adamupdegrove/Desktop/tmp/BEFOREEDGES_"+std::to_string(surfer)+"_"+std::to_string(r)+".vtp";
    //  vtkSVIOUtils::WriteVTPFile(fn, writePointsPd);
    //}

    //edgeExp.Init(newShape, TopAbs_EDGE);

    //fprintf(stdout,"\n");
    //fprintf(stdout,"NEW EDGES: \n");
    //int numEdges = 0;
    //for (int r=0; edgeExp.More(); edgeExp.Next(), r++)
    //{
    //  TopoDS_Edge thisEdge = TopoDS::Edge(edgeExp.Current());

    //  Standard_Real paraFirst, paraLast;
    //  Handle(Geom_Curve) thisCurve = BRep_Tool::Curve (thisEdge, paraFirst, paraLast);
    //  Standard_Real faceFirst, faceLast;
    //  Handle(Geom2d_Curve) aPCurve = BRep_Tool::CurveOnSurface (thisEdge, face, faceFirst, faceLast);

    //  fprintf(stdout,"EDGE %d HAS ORIEN %d\n", r, thisEdge.Orientation());
    //  fprintf(stdout,"EDGE %d HAS RANGE %.6f %.6f\n", r, paraFirst, paraLast);
    //  fprintf(stdout,"EDGE %d HAS RANGE ON FACE %.6f %.6f\n", r, faceFirst, faceLast);

    //  gp_Pnt2d pntFirst, pntLast;
    //  BRep_Tool::UVPoints(thisEdge, face, pntFirst, pntLast);

    //  fprintf(stdout,"EDGE %d HAS STARTING UV VALUE: %.6f %.6f\n", r, pntFirst.X(), pntFirst.Y());
    //  fprintf(stdout,"EDGE %d HAS ENDERING UV VALUE: %.6f %.6f\n", r, pntLast.X(), pntLast.Y());

    //  numEdges++;

    //  vtkNew(vtkPoints, writePoints);
    //  for (int s=0; s<20; s++)
    //  {
    //    Standard_Real pointVal = paraFirst + (s/20.)*(paraLast-paraFirst);
    //    gp_Pnt writePoint = thisCurve->Value(pointVal);

    //    writePoints->InsertNextPoint(writePoint.X(), writePoint.Y(), writePoint.Z());
    //  }

    //  vtkNew(vtkPolyData, writePointsPd);
    //  writePointsPd->SetPoints(writePoints);

    //  std::string fn = "/Users/adamupdegrove/Desktop/tmp/AFTEREEDGES_"+std::to_string(surfer)+"_"+std::to_string(r)+".vtp";
    //  vtkSVIOUtils::WriteVTPFile(fn, writePointsPd);

    //  fprintf(stdout,"ANALYZING EDGE\n");
    //  OCCTUtils_AnalyzeShape(thisEdge);

    //  //TopExp_Explorer vertexExp;
    //  //vertexExp.Init(thisEdge, TopAbs_VERTEX);
    //  //fprintf(stdout,"\n");
    //  //for (int l=0; vertexExp.More(); vertexExp.Next(), l++)
    //  //{
    //  //  TopoDS_Vertex thisVertex = TopoDS::Vertex(vertexExp.Current());
    //  //  gp_Pnt thisVertexPoint = BRep_Tool::Pnt(thisVertex);
    //  //  fprintf(stdout,"THIS VERTEX POSITION: %.6f %.6f %.6f\n", thisVertexPoint.X(), thisVertexPoint.Y(), thisVertexPoint.Z());

    //  //  fprintf(stdout, "POINT %d ORIEN: %d\n", l, thisVertex.Orientation());
    //  //}

    //}
    //fprintf(stdout,"TOT NUM OF EDGES: %d\n", numEdges);

    shape = newShape;

    //Standard_Real length = 8.0;
    //Standard_Real radius = 2.0;
    //TopoDS_Solid cyl = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 0.1)), radius, length).Solid();

    //TopoDS_Face face = TopoDS::Face(OCCTUtils_GetFirstType(cyl, TopAbs_FACE));
    //TopoDS_Edge edge = TopoDS::Edge(OCCTUtils_GetFirstType(cyl, TopAbs_EDGE));

    //TopExp_Explorer edgeExp;
    //edgeExp.Init(cyl, TopAbs_EDGE);

    //for (int r=0; edgeExp.More(); edgeExp.Next(), r++)
    //{
    //  TopoDS_Edge thisEdge = TopoDS::Edge(edgeExp.Current());

    //  Standard_Real paraFirst, paraLast;
    //  Handle(Geom_Curve) thisCurve = BRep_Tool::Curve (thisEdge, paraFirst, paraLast);
    //  Standard_Real faceFirst, faceLast;
    //  Handle(Geom2d_Curve) aPCurve = BRep_Tool::CurveOnSurface (thisEdge, face, faceFirst, faceLast);

    //  fprintf(stdout,"EDGE %d HAS ORIEN %d\n", r, thisEdge.Orientation());
    //  fprintf(stdout,"EDGE %d HAS RANGE %.6f %.6f\n", r, paraFirst, paraLast);
    //  fprintf(stdout,"EDGE %d HAS RANGE ON FACE %.6f %.6f\n", r, faceFirst, faceLast);

    //  gp_Pnt2d pntFirst, pntLast;
    //  BRep_Tool::UVPoints(thisEdge, face, pntFirst, pntLast);

    //  fprintf(stdout,"EDGE %d HAS STARTING UV VALUE: %.6f %.6f\n", r, pntFirst.X(), pntFirst.Y());
    //  fprintf(stdout,"EDGE %d HAS ENDERING UV VALUE: %.6f %.6f\n", r, pntLast.X(), pntLast.Y());
    //}

    //BRep_Builder builder;

    //Standard_Real pFirst, pLast;
    //Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, pFirst, pLast);
    //gp_Pnt midPoint = curve->Value(pFirst + (pLast-pFirst)/2);

    //TopoDS_Vertex vStart, vEnd, vMiddle;
    //vStart = TopExp::FirstVertex(edge);
    //vEnd = TopExp::LastVertex(edge);
    //builder.MakeVertex(vMiddle, midPoint, Precision::Confusion());

    //TopoDS_Edge newEdge1 = BRepBuilderAPI_MakeEdge(curve, vStart, TopoDS::Vertex(vMiddle.Reversed()));
    //TopoDS_Edge newEdge2 = BRepBuilderAPI_MakeEdge(curve, vMiddle, TopoDS::Vertex(vEnd.Reversed()));

    //TopoDS_Wire wire;
    //builder.MakeWire(wire);
    //newEdge1.Orientation(TopAbs_REVERSED);
    //builder.Add(wire, newEdge1);
    //newEdge2.Orientation(TopAbs_REVERSED);
    //builder.Add(wire, newEdge2);

    //builder.UpdateEdge(newEdge1, new Geom2d_Line(gp_Pnt2d(0.0, length), gp_Dir2d(1,0)), face, Precision::Confusion());
    //builder.Range(newEdge1, face, 0, pFirst + (pLast-pFirst)/2);
    //builder.UpdateEdge(newEdge2, new Geom2d_Line(gp_Pnt2d(0.0, length), gp_Dir2d(1,0)), face, Precision::Confusion());
    //builder.Range(newEdge2, face, pFirst + (pLast-pFirst)/2, pLast);

    //Handle(ShapeBuild_ReShape) context = new ShapeBuild_ReShape();
    //context->Replace(edge, wire);
    //TopoDS_Shape output = context->Apply(cyl);
    //*(loftedSurfs[surfer]->geom_) = output;
    //OCCTUtils_AnalyzeShape(*(loftedSurfs[surfer]->geom_));

    //edgeExp.Init(output, TopAbs_EDGE);

    //fprintf(stdout,"\n");
    //fprintf(stdout,"NEW EDGES: \n");
    //int numEdges = 0;
    //for (int r=0; edgeExp.More(); edgeExp.Next(), r++)
    //{
    //  TopoDS_Edge thisEdge = TopoDS::Edge(edgeExp.Current());

    //  Standard_Real paraFirst, paraLast;
    //  Handle(Geom_Curve) thisCurve = BRep_Tool::Curve (thisEdge, paraFirst, paraLast);
    //  Standard_Real faceFirst, faceLast;
    //  Handle(Geom2d_Curve) aPCurve = BRep_Tool::CurveOnSurface (thisEdge, face, faceFirst, faceLast);

    //  fprintf(stdout,"EDGE %d HAS ORIEN %d\n", r, thisEdge.Orientation());
    //  fprintf(stdout,"EDGE %d HAS RANGE %.6f %.6f\n", r, paraFirst, paraLast);
    //  fprintf(stdout,"EDGE %d HAS RANGE ON FACE %.6f %.6f\n", r, faceFirst, faceLast);

    //  gp_Pnt2d pntFirst, pntLast;
    //  BRep_Tool::UVPoints(thisEdge, face, pntFirst, pntLast);

    //  fprintf(stdout,"EDGE %d HAS STARTING UV VALUE: %.6f %.6f\n", r, pntFirst.X(), pntFirst.Y());
    //  fprintf(stdout,"EDGE %d HAS ENDERING UV VALUE: %.6f %.6f\n", r, pntLast.X(), pntLast.Y());

    //  numEdges++;
    //}
    //fprintf(stdout,"TOT NUM OF EDGES: %d\n", numEdges);

    //std::cout << std::endl << "Program finished normally" << std::endl;
  }
  catch (Standard_Failure)
  {
    Handle_Standard_Failure e = Standard_Failure::Caught();
    std::cout << "OCC Error: " << e->GetMessageString() << std::endl;
    return SV_ERROR;
  }
  catch (const std::exception &error)
  {
    std::cout << "Error: " << error.what() << std::endl;
    return SV_ERROR;
  }

  return SV_OK;
}

// ---------------------
// OCCTUtils_GetOpenEdge
// ---------------------
/**
 * @brief
 * @param
 * @return SV_OK if function completes properly
 */
int OCCTUtils_GetOpenEdge(TopoDS_Shape &shape, const int edgeNumber, TopoDS_Edge &edge)
{
  Standard_Real sewtoler =  1.e-6;
  Standard_Real closetoler =  1.e-4;
  ShapeFix_FreeBounds findFree(shape,sewtoler,closetoler,
            Standard_False,Standard_False);
  TopoDS_Compound freeWires = findFree.GetClosedWires();
  TopExp_Explorer NewEdgeExp;
  NewEdgeExp.Init(freeWires,TopAbs_EDGE);

  int foundEdge = 0;
  int numEdges = 0;
  for (int i=0;NewEdgeExp.More();NewEdgeExp.Next(),i++)
  {
    if (i == edgeNumber)
    {
      edge = TopoDS::Edge(NewEdgeExp.Current());
      foundEdge = 1;
      break;
    }
    numEdges++;
  }

  if (!foundEdge)
  {
    std::cerr << "Error getting edge " << edgeNumber << " off of shape. There are only " << numEdges << " open edges on the shape.";
  }

  return SV_OK;
}

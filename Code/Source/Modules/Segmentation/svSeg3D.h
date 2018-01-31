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

#ifndef SVSEG3D_H
#define SVSEG3D_H

#include <svSegmentationExports.h>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <map>

struct SVSEGMENTATION_EXPORT svSeed
{
    int id;
    std::string type;//begin, end

    double x;
    double y;
    double z;

    double radius;

    bool selected;

    std::string status;

    svSeed()
        : id(-1)
        , type("")
        , x(0)
        , y(0)
        , z(0)
        , radius(0.2)
        , selected(false)
        , status("")
    {
    }

    svSeed(double xx, double yy, double zz, std::string ttype="")
        : id(-1)
        , type(ttype)
        , x(xx)
        , y(yy)
        , z(zz)
        , radius(0.2)
        , selected(false)
        , status("")
    {
    }

    svSeed(double xx, double yy, double zz, double r, std::string ttype="")
        : id(-1)
        , type(ttype)
        , x(xx)
        , y(yy)
        , z(zz)
        , radius(r)
        , selected(false)
        , status("")
    {
    }

    svSeed(const svSeed &other)
        : id(other.id)
        , type(other.type)
        , x(other.x)
        , y(other.y)
        , z(other.z)
        , radius(other.radius)
        , selected(false)
        , status("")
    {
    }

};

struct SVSEGMENTATION_EXPORT svSeg3DParam
{
    std::string method;

    double lowerThreshold;
    double upperThreshold;

    std::map<int, svSeed> seedMap;

    svSeg3DParam()
        : method("")
        , lowerThreshold(0)
        , upperThreshold(0)
    {
    }

    svSeg3DParam(const svSeg3DParam &other)
        : method(other.method)
        , lowerThreshold(other.lowerThreshold)
        , upperThreshold(other.upperThreshold)
        , seedMap(other.seedMap)
    {
    }

    std::map<int,svSeed>& GetSeedMap()
    {
        return seedMap;
    }

    int AddSeed(svSeed seed)
    {    int newID=seed.id;

         if(newID<0)
         {
             int idmax=0;
             for(auto s:seedMap)
             {
                 if(s.first>idmax)
                     idmax=s.first;
             }

             newID=idmax+1;

             seed.id=newID;
         }

         seedMap[newID]=seed;

         return newID;
    }

    void RemoveSeed(int id)
    {
       seedMap.erase(id);
    }

};



class SVSEGMENTATION_EXPORT svSeg3D
{
public:

    svSeg3D();

    svSeg3D(const svSeg3D &other, bool copyVpd=true);

    virtual ~svSeg3D();

    virtual svSeg3D* Clone();

    svSeg3DParam& GetParam();

    svSeg3DParam& GetInnerParam();

    void SetParam(svSeg3DParam param, bool copyToInner=true);

    vtkSmartPointer<vtkPolyData> GetVtkPolyData(){return m_Vpd;}

    void SetVtkPolyData(vtkSmartPointer<vtkPolyData> vpd) {m_Vpd=vpd;}

protected:

    svSeg3DParam m_Param;

    svSeg3DParam m_InnerParam;

    vtkSmartPointer<vtkPolyData> m_Vpd;

};


#endif // SVSEG3D_H

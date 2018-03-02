# Copyright (c) Stanford University, The Regents of the University of
#               California, and others.
#
# All Rights Reserved.
#
# See Copyright-SimVascular.txt for additional details.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

set(H_FILES
    svSegmentationUtils.h
    svContour.h
    svContourCircle.h
    svContourEllipse.h
    svContourPolygon.h
    svContourTensionPolygon.h
    svContourSplinePolygon.h
    svContourOperation.h
    svContourModel.h
    svContourModelVtkMapper2D.h
    svContourModelThresholdInteractor.h
    svContourGroup.h
    svContourGroupVtkMapper2D.h
    svContourGroupVtkMapper3D.h
    svContourGroupDataInteractor.h
    svContourGroupIO.h
    svSegmentationLegacyIO.h
    svSurface.h
    svSurfaceVtkMapper3D.h
    svSeg3D.h
    svMitkSeg3D.h
    svMitkSeg3DOperation.h
    svMitkSeg3DIO.h
    svMitkSeg3DVtkMapper3D.h
    svMitkSeg3DDataInteractor.h
    svSegmentationObjectFactory.h
    svSeg3DUtils.h
)

set(CPP_FILES
    svSegmentationUtils.cxx
    svContour.cxx
    svContourCircle.cxx
    svContourEllipse.cxx
    svContourPolygon.cxx
    svContourTensionPolygon.cxx
    svContourSplinePolygon.cxx
    svContourOperation.cxx
    svContourModel.cxx
    svContourModelVtkMapper2D.cxx
    svContourModelThresholdInteractor.cxx
    svContourGroup.cxx
    svContourGroupVtkMapper2D.cxx
    svContourGroupVtkMapper3D.cxx
    svContourGroupDataInteractor.cxx
    svContourGroupIO.cxx
    svSegmentationLegacyIO.cxx
    svSurface.cxx
    svSurfaceVtkMapper3D.cxx
    svSeg3D.cxx
    svMitkSeg3D.cxx
    svMitkSeg3DOperation.cxx
    svMitkSeg3DIO.cxx
    svMitkSeg3DVtkMapper3D.cxx
    svMitkSeg3DDataInteractor.cxx
    svSegmentationObjectFactory.cxx
    svSeg3DUtils.cxx
)

set(RESOURCE_FILES
    Interactions/svContourGroupInteraction.xml
    Interactions/svContourModelThresholdInteraction.xml
    Interactions/svMitkSeg3DInteraction.xml
    Interactions/svSegmentationConfig.xml
)

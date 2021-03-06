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

TCL_MAJOR_VERSION=8
TCL_MINOR_VERSION=6
TCL_PATCH_VERSION=8
TCL_VERSION=$(TCL_MAJOR_VERSION).$(TCL_MINOR_VERSION).$(TCL_PATCH_VERSION)

TK_MAJOR_VERSION=8
TK_MINOR_VERSION=6
TK_PATCH_VERSION=8
TK_VERSION=$(TK_MAJOR_VERSION).$(TK_MINOR_VERSION).$(TK_PATCH_VERSION)

TCLLIB_VERSION=1.17
TKLIB_VERSION=0.6

BUILDFLAGS     += -D__NON_STD_TCL_INSTALL
TCL_BASE       = $(OPEN_SOFTWARE_BINARIES_TOPLEVEL)/tcltk-$(TCL_VERSION)
TK_BASE        = $(OPEN_SOFTWARE_BINARIES_TOPLEVEL)/tcltk-$(TK_VERSION)
TCLTK_INCDIR   = -I$(TCL_BASE)/include -I$(TK_BASE)/include
TCLTK_LIBDIR   = -L$(TCL_BASE)/lib -L$(TK_BASE)/lib
TCLTK_DLLS     = $(TCL_BASE)/bin/tcl$(TCL_MAJOR_VERSION).$(TCL_MINOR_VERSION).$(SOEXT) $(TCL_BASE)/bin/tk$(TK_MAJOR_VERSION).$(TK_MINOR_VERSION).$(SOEXT)
TCLTK_LIBS     = $(TCLTK_LIBDIR) -ltcl$(TCL_MAJOR_VERSION).$(TCL_MINOR_VERSION) -ltk$(TK_MAJOR_VERSION).$(TK_MINOR_VERSION)
TKCXIMAGE_BASE = $(OPEN_SOFTWARE_BINARIES_TOPLEVEL)/tkcximage-0.98.9/tcltk-$(TK_VERSION)
TKCXIMAGE_DLL  = $(TKCXIMAGE_BASE)/bin/Tkcximage.$(SOEXT)
TCLTK_SO_PATH  = $(TCL_BASE)/lib
TCL_LIBRARY    = $(TCL_BASE)/lib/tcl$(TCL_MAJOR_VERSION).$(TCL_MINOR_VERSION)
TK_LIBRARY     = $(TCL_BASE)/lib/tk$(TK_MAJOR_VERSION).$(TK_MINOR_VERSION)
TCLSH          = $(TCL_BASE)/bin/tclsh$(TCL_MAJOR_VERSION).$(TCL_MINOR_VERSION)

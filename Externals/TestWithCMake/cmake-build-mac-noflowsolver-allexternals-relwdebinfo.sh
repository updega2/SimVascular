rm -Rf ../../Code/SeparateSandBoxBuild
mkdir -p ../../Code/SeparateSandBoxBuild
pushd ../../Code/SeparateSandBoxBuild

#compilers
export CC="clang"
export CXX="clang++"

#cmake
export REPLACEME_SV_CMAKE_CMD="cmake"
export REPLACEME_SV_CMAKE_GENERATOR="Unix Makefiles"
export REPLACEME_SV_CMAKE_BUILD_TYPE="RelWithDebInfo"
export REPLACEME_SV_MAKE_CMD="make -j8"
export REPLACEME_SV_TOP_SRC_DIR_SV="../"

#externals
export REPLACEME_SV_EXTERNALS_TOPLEVEL_DIR="/usr/local/sv/ext"

"$REPLACEME_SV_CMAKE_CMD" \
\
    --debug-output \
\
   -G "$REPLACEME_SV_CMAKE_GENERATOR" \
\
   -DCMAKE_BUILD_TYPE="$REPLACEME_SV_CMAKE_BUILD_TYPE" \
   -DBUILD_SHARED_LIBS=ON \
   -DBUILD_TESTING=OFF \
\
   -DSV_SUPERBUILD=OFF \
\
   -DSV_USE_FREETYPE=ON \
   -DSV_USE_GDCM=ON \
   -DSV_USE_ITK=ON \
   -DSV_USE_MPICH2=OFF \
   -DSV_USE_OpenCASCADE=ON \
   -DSV_USE_MMG=ON \
\
   -DSV_USE_SYSTEM_FREETYPE=ON \
   -DSV_USE_SYSTEM_GDCM=ON \
   -DSV_USE_SYSTEM_ITK=ON \
   -DSV_USE_SYSTEM_OpenCASCADE=ON \
   -DSV_USE_SYSTEM_TCL=ON \
   -DSV_USE_SYSTEM_VTK=ON \
   -DSV_USE_SYSTEM_MMG=ON \
\
   -DSV_USE_GDCM_SHARED=ON \
   -DSV_USE_FREETYPE_SHARED=ON \
   -DSV_USE_ITK_SHARED=ON \
   -DSV_USE_OpenCASCADE_SHARED=ON \
   -DSV_USE_TCL_SHARED=ON \
   -DSV_USE_VTK_SHARED=ON \
\
   -DSV_EXTERNALS_USE_TOPLEVEL_DIR=ON \
   -DSV_EXTERNALS_TOPLEVEL_DIR="$REPLACEME_SV_EXTERNALS_TOPLEVEL_DIR" \
\
 "$REPLACEME_SV_TOP_SRC_DIR_SV" >& stdout-cmake-config.txt

$REPLACEME_SV_MAKE_CMD >& stdout-compile.txt

popd

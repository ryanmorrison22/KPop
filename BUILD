#!/usr/bin/env bash

set -euo pipefail

PROFILE="${1:-dev}"

# Default compilers (override via env if needed)
CC="${CC:-cc}"
CXX="${CXX:-c++}"
FC="${FC:-gfortran}"

# Detect platform
OS="$(uname -s)"
ARCH="$(uname -m)"

# Portable CPU count
if command -v nproc >/dev/null 2>&1; then
  NPROC=$(nproc)
else
  NPROC=$(sysctl -n hw.ncpu)
fi

# BLAS target
if [[ -z "${BLAS_TARGET:-}" ]]; then
  if [[ "$ARCH" == "arm64" || "$ARCH" == "aarch64" ]]; then
    BLAS_TARGET="ARMV8"
  else
    BLAS_TARGET="HASWELL"
  fi
fi

echo "== Build configuration =="
echo "PROFILE=$PROFILE"
echo "CC=$CC"
echo "CXX=$CXX"
echo "FC=$FC"
echo "NPROC=$NPROC"
echo "BLAS_TARGET=$BLAS_TARGET"
echo "ARCH=$ARCH"
echo "OS=$OS"
echo "=========================="

# Clean OCaml build
rm -rf _build

# Keep build cache dir but ensure it exists
mkdir -p build
mkdir -p lib

rm -f lib/libopenblas.a lib/libfaiss.a lib/libinterfaiss.a

########################################
# OpenBLAS
########################################
if [[ -f OpenBLAS/libopenblas.a ]]; then
  echo "Using cached OpenBLAS"
  cp OpenBLAS/libopenblas.a lib/
else
  echo "Building OpenBLAS..."
  pushd OpenBLAS

  if [[ "$OS" == "Darwin" ]]; then
    EXTRA_FFLAGS="-fno-lto"
  else
    EXTRA_FFLAGS=""
  fi

  make -j "$NPROC" libs netlib \
    CC="$CC" \
    FC="$FC" \
    HOSTCC="$CC" \
    TARGET="$BLAS_TARGET" \
    NO_AVX="${NO_AVX:-0}" \
    USE_OPENMP="${USE_OPENMP:-1}" \
    FFLAGS="$EXTRA_FFLAGS"

  cp libopenblas.a ../lib/
  popd
fi

########################################
# FAISS
########################################

# Disable AVX on ARM
if [[ "$ARCH" == "arm64" || "$ARCH" == "aarch64" ]]; then
  FAISS_OPT_LEVEL="generic"
  FAISS_TARGET="faiss"
else
  FAISS_OPT_LEVEL="avx2"
  FAISS_TARGET="faiss_avx2"
fi

if [[ -f build/faiss/faiss/lib${FAISS_TARGET}.a ]]; then
  echo "Using cached FAISS"
  cp build/faiss/faiss/lib${FAISS_TARGET}.a lib/libfaiss.a
else
  echo "Building FAISS..."
  mkdir -p build/faiss

  pushd faiss

  if [[ "$OS" == "Darwin" ]]; then
    LIBOMP_PREFIX=$(brew --prefix libomp)
    OMP_FLAGS=(
      -D OpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include"
      -D OpenMP_CXX_LIB_NAMES="omp"
      -D OpenMP_omp_LIBRARY="${LIBOMP_PREFIX}/lib/libomp.a"
    )
  else
    OMP_FLAGS=()
  fi

  cmake \
    -D CMAKE_VERBOSE_MAKEFILE=true \
    -D CMAKE_CXX_COMPILER="$CXX" \
    -D BLAS_LIBRARIES="$(pwd)/../lib/libopenblas.a" \
    -D LAPACK_LIBRARIES="$(pwd)/../lib/libopenblas.a" \
    -D FAISS_ENABLE_GPU=false \
    -D FAISS_ENABLE_PYTHON=false \
    -D BUILD_TESTING=false \
    -D CMAKE_BUILD_TYPE=Release \
    -D FAISS_OPT_LEVEL="$FAISS_OPT_LEVEL" \
    "${OMP_FLAGS[@]}" \
    -B ../build/faiss .

  popd

  pushd build/faiss
  make -j "$NPROC" "$FAISS_TARGET"
  cp faiss/lib${FAISS_TARGET}.a ../../lib/libfaiss.a
  popd
fi

########################################
# interfaiss
########################################
echo "Building interfaiss..."
pushd lib

if [[ "$OS" == "Darwin" ]]; then
  LIBOMP_PREFIX=$(brew --prefix libomp)
  "$CXX" -std=c++17 -I ../faiss/ -O3 -fPIC \
    -Xpreprocessor -fopenmp -I"${LIBOMP_PREFIX}/include" \
    -c interfaiss.cpp -o libinterfaiss.o
else
  "$CXX" -std=c++17 -I ../faiss/ -O3 -fPIC -fopenmp -c interfaiss.cpp -o libinterfaiss.o
fi
ar rcs libinterfaiss.a libinterfaiss.o
rm -f libinterfaiss.o

popd

########################################
# Version info (portable date)
########################################

LAST_COMMIT_DATE=$(git log -1 --format="%cd" --date=format:"%d-%b-%Y")
VERSION=$(git rev-list --count HEAD)

pushd BiOCamLib
cat > lib/Info.ml <<EOF
include (
  struct
    let info = {
      Tools.Argv.name = "BiOCamLib";
      version = "$VERSION";
      date = "$LAST_COMMIT_DATE"
    }
  end
)
EOF
popd

cat > lib/Info.ml <<EOF
include (
  struct
    let info = {
      BiOCamLib.Tools.Argv.name = "KPop";
      version = "$VERSION";
      date = "$LAST_COMMIT_DATE"
    }
  end
)
EOF

########################################
# OCaml build
########################################

dune build --profile="$PROFILE" bin/KPopCount.exe
dune build --profile="$PROFILE" bin/KPopCountDB.exe
dune build --profile="$PROFILE" bin/KPopTwist_.exe
dune build --profile="$PROFILE" bin/KPopTwistDB.exe

########################################
# Collect binaries
########################################

mv _build/default/bin/KPopCount.exe build/KPopCount
mv _build/default/bin/KPopCountDB.exe build/KPopCountDB
mv _build/default/bin/KPopTwist_.exe build/KPopTwist_
mv _build/default/bin/KPopTwistDB.exe build/KPopTwistDB

chmod 755 build/*

########################################
# Strip (release only)
########################################

if [[ "$PROFILE" == "release" || "$PROFILE" == "release-static" ]]; then
  strip build/* || true
  rm -rf _build
fi

cp src/KPop* build

echo "✅ Build completed successfully"

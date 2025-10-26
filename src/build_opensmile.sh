#!/bin/bash
# Build OpenSMILE static library for R package integration
# This script builds OpenSMILE as a static library that can be linked by the R package

set -e  # Exit on error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OPENSMILE_DIR="$SCRIPT_DIR/opensmile"
BUILD_DIR="$OPENSMILE_DIR/build_r"

echo "===================================================="
echo "Building OpenSMILE static library for R package"
echo "===================================================="

# Check if OpenSMILE directory exists
if [ ! -d "$OPENSMILE_DIR" ]; then
    echo "ERROR: OpenSMILE directory not found at $OPENSMILE_DIR"
    echo "Make sure the OpenSMILE submodule is initialized:"
    echo "  git submodule update --init --recursive"
    exit 1
fi

# Create build directory
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

echo "Configuring OpenSMILE with CMake..."

# Configure with minimal dependencies
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DSTATIC_LINK=ON \
    -DWITH_PORTAUDIO=OFF \
    -DWITH_FFMPEG=OFF \
    -DWITH_OPENSLES=OFF \
    -DWITH_OPENCV=OFF \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    ..

echo "Building OpenSMILE..."
cmake --build . --config Release -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

echo "===================================================="
echo "OpenSMILE build complete!"
echo "Static libraries created in: $BUILD_DIR"
echo "===================================================="

# List created libraries
echo "Built libraries:"
ls -lh "$BUILD_DIR"/*.a 2>/dev/null || ls -lh "$BUILD_DIR/Release"/*.lib 2>/dev/null || echo "Library files not found in expected locations"

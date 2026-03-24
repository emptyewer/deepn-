#!/bin/bash

# Build script for DESeq2 project (Statistics Library + GUI)

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Building DESeq2 Project (Statistics Library + GUI)...${NC}"

# Check if we're in the right directory
if [ ! -f "CMakeLists.txt" ]; then
    echo -e "${RED}Error: CMakeLists.txt not found. Please run this script from the deseq2 directory.${NC}"
    exit 1
fi

# Create build directory
BUILD_DIR="build"
if [ ! -d "$BUILD_DIR" ]; then
    echo -e "${YELLOW}Creating build directory...${NC}"
    mkdir -p "$BUILD_DIR"
else
    echo -e "${GREEN}Build directory already exists.${NC}"
    rm -rf "$BUILD_DIR"
    mkdir -p "$BUILD_DIR"
fi

# Navigate to build directory
cd "$BUILD_DIR"

# Configure with CMake
echo -e "${YELLOW}Configuring with CMake...${NC}"
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
echo -e "${YELLOW}Building...${NC}"
make -j$(nproc)

echo -e "${GREEN}Build completed successfully!${NC}"
echo -e "${GREEN}Executable location: $BUILD_DIR/bin/DESeq2.app${NC}"

# Check if Qt5 is available
if command -v qmake &>/dev/null; then
    QT_VERSION=$(qmake -query QT_VERSION)
    echo -e "${GREEN}Qt version: $QT_VERSION${NC}"
else
    echo -e "${YELLOW}Warning: qmake not found. Make sure Qt5 is properly installed.${NC}"
fi

# Check if Eigen3 is available
if pkg-config --exists eigen3; then
    EIGEN_VERSION=$(pkg-config --modversion eigen3)
    echo -e "${GREEN}Eigen3 version: $EIGEN_VERSION${NC}"
else
    echo -e "${YELLOW}Warning: Eigen3 not found via pkg-config.${NC}"
fi

# Check if Boost is available
if pkg-config --exists boost; then
    BOOST_VERSION=$(pkg-config --modversion boost)
    echo -e "${GREEN}Boost version: $BOOST_VERSION${NC}"
else
    echo -e "${YELLOW}Warning: Boost not found via pkg-config.${NC}"
fi

echo -e "${BLUE}Project structure:${NC}"
echo -e "  ${GREEN}✓${NC} Statistics Library (deseq2_statistics)"
echo -e "  ${GREEN}✓${NC} GUI Application (deseq2_gui)"
echo -e "  ${GREEN}✓${NC} Test Executables"
echo -e "  ${GREEN}✓${NC} Example Executables"

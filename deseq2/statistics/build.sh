#!/bin/bash

# DESeq2 C++ Library Build Script
# This script builds the DESeq2 library with various options and configurations

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values
BUILD_TYPE="Release"
BUILD_DIR="build"
CLEAN_BUILD=false
RUN_TESTS=false
RUN_EXAMPLES=false
INSTALL=false
INSTALL_PREFIX="/usr/local"
VERBOSE=false
JOBS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to show usage
show_usage() {
    cat << EOF
DESeq2 C++ Library Build Script

Usage: $0 [OPTIONS]

Options:
    -h, --help              Show this help message
    -c, --clean             Clean build directory before building
    -d, --debug             Build in Debug mode (default: Release)
    -t, --test              Run tests after building
    -e, --examples          Run examples after building
    -i, --install           Install library after building
    -p, --prefix DIR        Installation prefix (default: /usr/local)
    -j, --jobs N            Number of parallel jobs (default: auto-detect)
    -v, --verbose           Verbose output
    -b, --build-dir DIR     Build directory (default: build)

Examples:
    $0                      # Build in Release mode
    $0 -c -d               # Clean build in Debug mode
    $0 -t -e               # Build and run tests and examples
    $0 -i -p /opt/deseq2   # Build and install to /opt/deseq2

EOF
}

# Function to detect OS
detect_os() {
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        OS="linux"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        OS="macos"
    elif [[ "$OSTYPE" == "cygwin" ]] || [[ "$OSTYPE" == "msys" ]]; then
        OS="windows"
    else
        OS="unknown"
    fi
    print_status "Detected OS: $OS"
}

# Function to check dependencies
check_dependencies() {
    print_status "Checking dependencies..."
    
    # Check CMake
    if ! command -v cmake &> /dev/null; then
        print_error "CMake not found. Please install CMake 3.16 or higher."
        exit 1
    fi
    
    CMAKE_VERSION=$(cmake --version | head -n1 | cut -d' ' -f3)
    print_status "CMake version: $CMAKE_VERSION"
    
    # Check C++ compiler
    if command -v g++ &> /dev/null; then
        COMPILER="g++"
        COMPILER_VERSION=$(g++ --version | head -n1)
    elif command -v clang++ &> /dev/null; then
        COMPILER="clang++"
        COMPILER_VERSION=$(clang++ --version | head -n1)
    else
        print_error "No C++ compiler found. Please install g++ or clang++."
        exit 1
    fi
    
    print_status "Using compiler: $COMPILER_VERSION"
    
    # Check for Eigen3
    if [[ "$OS" == "macos" ]]; then
        if brew list eigen &> /dev/null; then
            print_success "Eigen3 found via Homebrew"
        else
            print_warning "Eigen3 not found. Installing via Homebrew..."
            brew install eigen
        fi
    elif [[ "$OS" == "linux" ]]; then
        if pkg-config --exists eigen3; then
            print_success "Eigen3 found via pkg-config"
        else
            print_warning "Eigen3 not found. Please install libeigen3-dev"
            print_status "On Ubuntu/Debian: sudo apt-get install libeigen3-dev"
            print_status "On CentOS/RHEL: sudo yum install eigen3-devel"
        fi
    fi
    
    # Check for Boost
    if [[ "$OS" == "macos" ]]; then
        if brew list boost &> /dev/null; then
            print_success "Boost found via Homebrew"
        else
            print_warning "Boost not found. Installing via Homebrew..."
            brew install boost
        fi
    elif [[ "$OS" == "linux" ]]; then
        if pkg-config --exists boost; then
            print_success "Boost found via pkg-config"
        else
            print_warning "Boost not found. Please install libboost-all-dev"
            print_status "On Ubuntu/Debian: sudo apt-get install libboost-all-dev"
            print_status "On CentOS/RHEL: sudo yum install boost-devel"
        fi
    fi
}

# Function to clean build directory
clean_build() {
    if [[ "$CLEAN_BUILD" == true ]]; then
        print_status "Cleaning build directory..."
        if [[ -d "$BUILD_DIR" ]]; then
            rm -rf "$BUILD_DIR"
            print_success "Build directory cleaned"
        fi
    fi
}

# Function to create build directory
create_build_dir() {
    print_status "Creating build directory: $BUILD_DIR"
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
}

# Function to configure with CMake
configure_cmake() {
    print_status "Configuring with CMake..."
    
    CMAKE_ARGS=(
        "-DCMAKE_BUILD_TYPE=$BUILD_TYPE"
        "-DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX"
    )
    
    if [[ "$VERBOSE" == true ]]; then
        CMAKE_ARGS+=("--verbose")
    fi
    
    # Add platform-specific options
    if [[ "$OS" == "macos" ]]; then
        # On macOS, help CMake find Homebrew packages
        CMAKE_ARGS+=(
            "-DEigen3_DIR=$(brew --prefix eigen)/lib/cmake/eigen3"
            "-DBOOST_ROOT=$(brew --prefix boost)"
        )
    fi
    
    cmake "${CMAKE_ARGS[@]}" ..
    
    if [[ $? -eq 0 ]]; then
        print_success "CMake configuration successful"
    else
        print_error "CMake configuration failed"
        exit 1
    fi
}

# Function to build the project
build_project() {
    print_status "Building DESeq2 library..."
    print_status "Build type: $BUILD_TYPE"
    print_status "Jobs: $JOBS"
    
    if [[ "$VERBOSE" == true ]]; then
        make -j"$JOBS" VERBOSE=1
    else
        make -j"$JOBS"
    fi
    
    if [[ $? -eq 0 ]]; then
        print_success "Build completed successfully"
    else
        print_error "Build failed"
        exit 1
    fi
}

# Function to run tests
run_tests() {
    if [[ "$RUN_TESTS" == true ]]; then
        print_status "Running tests..."
        if [[ -f "test_deseq2" ]]; then
            ./test_deseq2
            if [[ $? -eq 0 ]]; then
                print_success "Tests passed"
            else
                print_error "Tests failed"
                exit 1
            fi
        else
            print_error "Test executable not found"
            exit 1
        fi
    fi
}

# Function to run examples
run_examples() {
    if [[ "$RUN_EXAMPLES" == true ]]; then
        print_status "Running examples..."
        if [[ -f "example_usage" ]]; then
            ./example_usage
            if [[ $? -eq 0 ]]; then
                print_success "Examples completed successfully"
            else
                print_error "Examples failed"
                exit 1
            fi
        else
            print_error "Example executable not found"
            exit 1
        fi
    fi
}

# Function to install the library
install_library() {
    if [[ "$INSTALL" == true ]]; then
        print_status "Installing library to $INSTALL_PREFIX..."
        
        # Check if we have write permissions
        if [[ ! -w "$INSTALL_PREFIX" ]] && [[ "$INSTALL_PREFIX" == "/usr/local" ]]; then
            print_warning "No write permission to $INSTALL_PREFIX. Using sudo..."
            sudo make install
        else
            make install
        fi
        
        if [[ $? -eq 0 ]]; then
            print_success "Library installed successfully"
            print_status "Headers installed to: $INSTALL_PREFIX/include"
            print_status "Library installed to: $INSTALL_PREFIX/lib"
        else
            print_error "Installation failed"
            exit 1
        fi
    fi
}

# Function to show build summary
show_summary() {
    print_status "Build Summary:"
    echo "  Build type: $BUILD_TYPE"
    echo "  Build directory: $BUILD_DIR"
    echo "  Compiler: $COMPILER"
    echo "  Jobs: $JOBS"
    echo "  Tests run: $RUN_TESTS"
    echo "  Examples run: $RUN_EXAMPLES"
    echo "  Installed: $INSTALL"
    if [[ "$INSTALL" == true ]]; then
        echo "  Install prefix: $INSTALL_PREFIX"
    fi
}

# Function to create pkg-config file
create_pkgconfig() {
    if [[ "$INSTALL" == true ]]; then
        print_status "Creating pkg-config file..."
        
        PC_FILE="$INSTALL_PREFIX/lib/pkgconfig/deseq2.pc"
        mkdir -p "$(dirname "$PC_FILE")"
        
        cat > "$PC_FILE" << EOF
prefix=$INSTALL_PREFIX
exec_prefix=\${prefix}
libdir=\${prefix}/lib
includedir=\${prefix}/include

Name: DESeq2
Description: C++ implementation of DESeq2 for differential expression analysis
Version: 1.0.0
Libs: -L\${libdir} -ldeseq2
Cflags: -I\${includedir}
Requires: eigen3
EOF
        
        print_success "pkg-config file created: $PC_FILE"
    fi
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_usage
            exit 0
            ;;
        -c|--clean)
            CLEAN_BUILD=true
            shift
            ;;
        -d|--debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        -t|--test)
            RUN_TESTS=true
            shift
            ;;
        -e|--examples)
            RUN_EXAMPLES=true
            shift
            ;;
        -i|--install)
            INSTALL=true
            shift
            ;;
        -p|--prefix)
            INSTALL_PREFIX="$2"
            shift 2
            ;;
        -j|--jobs)
            JOBS="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -b|--build-dir)
            BUILD_DIR="$2"
            shift 2
            ;;
        *)
            print_error "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Main build process
main() {
    print_status "Starting DESeq2 C++ Library build process..."
    
    # Store original directory
    ORIGINAL_DIR=$(pwd)
    
    # Detect OS
    detect_os
    
    # Check dependencies
    check_dependencies
    
    # Clean build if requested
    clean_build
    
    # Create and enter build directory
    create_build_dir
    
    # Configure with CMake
    configure_cmake
    
    # Build the project
    build_project
    
    # Run tests if requested
    run_tests
    
    # Run examples if requested
    run_examples
    
    # Install if requested
    install_library
    
    # Create pkg-config file if installing
    create_pkgconfig
    
    # Return to original directory
    cd "$ORIGINAL_DIR"
    
    # Show summary
    show_summary
    
    print_success "Build process completed successfully!"
    
    # Show usage instructions
    if [[ "$INSTALL" == true ]]; then
        echo
        print_status "To use the library in your project:"
        echo "  #include <deseq_dataset.h>"
        echo "  #include <deseq_stats.h>"
        echo "  #include <utils.h>"
        echo
        echo "  // Link with: -ldeseq2 -leigen3 -lboost_filesystem -lboost_system"
        echo "  // Or use pkg-config: pkg-config --cflags --libs deseq2"
    fi
}

# Run main function
main "$@" 
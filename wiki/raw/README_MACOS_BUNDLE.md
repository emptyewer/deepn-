# DESeq2 macOS Bundle - Static Qt5 Build

## Overview
Successfully created a self-contained DESeq2 GUI application as a macOS bundle (.app) using static Qt5 libraries. The application has no Qt5 runtime dependencies and can be distributed without requiring Qt5 installation on target systems.

## Build Configuration

### Static Qt5 Setup
- **Qt5 Static Path**: `/Users/venky/Softwares/qt5-static`
- **Qt Version**: 5.15.13
- **Build Type**: Hybrid Static/Dynamic (Qt5 static, system libraries dynamic)

### Application Details
- **Name**: DESeq2
- **Bundle Identifier**: `com.deseq2.app`
- **Version**: 1.0.0
- **Minimum macOS Version**: 10.14
- **Bundle Location**: `build/bin/DESeq2.app`

## Key Features

### ✅ Static Qt5 Integration
- Qt5 Core, Widgets, Gui, UiTools libraries embedded
- Qt5 plugins (image formats, platforms, styles) included
- No Qt5 runtime installation required

### ✅ macOS Bundle Structure
```
DESeq2.app/
├── Contents/
│   ├── Info.plist          # Bundle metadata
│   ├── MacOS/
│   │   └── DESeq2          # Main executable (15MB)
│   └── Resources/
│       └── DESeq2.icns     # Application icon
```

### ✅ Self-Contained Components
- DESeq2 Statistics Library (static)
- Eigen3 (header-only)
- All Qt5 dependencies embedded

### ✅ System Integration
- Proper macOS bundle structure
- Application icon support
- Document type associations (CSV, TXT, TSV)
- URL scheme support (`deseq2://`)

## Dependencies Analysis

### Static Libraries (Embedded)
- ✓ Qt5 Core, Widgets, Gui, UiTools
- ✓ Qt5 Plugins (image formats, platforms, styles)
- ✓ DESeq2 Statistics Library
- ✓ Eigen3 (header-only)

### Dynamic Dependencies (System)
- macOS System Frameworks (CoreFoundation, Cocoa, etc.)
- Boost Libraries (filesystem, system, atomic)
- Standard C++ Libraries (libc++, libSystem)

## Usage

### Running the Application
```bash
# Method 1: Direct execution
./build/bin/DESeq2.app/Contents/MacOS/DESeq2

# Method 2: Using macOS open command
open build/bin/DESeq2.app

# Method 3: Double-click in Finder
# Navigate to build/bin/ and double-click DESeq2.app
```

### Distribution
The `DESeq2.app` bundle can be:
- Copied to `/Applications/` for system-wide installation
- Distributed as a standalone application
- Packaged in a DMG for easy installation

## Build Process

### Prerequisites
1. Static Qt5 installation at `/Users/venky/Softwares/qt5-static`
2. CMake 3.16+
3. C++17 compatible compiler
4. Eigen3 and Boost libraries

### Build Commands
```bash
# Clean build
rm -rf build

# Build with static Qt5
./build.sh

# Create application icon (if needed)
cd ui && ./create_icon.sh
```

## Benefits

### ✅ Portability
- No Qt5 runtime dependencies
- Works on systems without Qt5 installed
- Consistent Qt5 version across systems

### ✅ Deployment
- Single .app bundle for distribution
- Smaller deployment package
- Professional macOS application appearance

### ✅ Maintenance
- No version conflicts with system Qt5
- Self-contained updates
- Simplified dependency management

## Limitations

### ⚠️ System Requirements
- Still requires macOS system frameworks
- Boost libraries remain dynamic (can be made static if needed)
- Larger executable size due to embedded Qt5

### ⚠️ Platform Specific
- macOS bundle format only
- Linux/Windows builds would need separate configurations

## Future Enhancements

### Potential Improvements
1. **Fully Static Build**: Include Boost libraries statically
2. **Universal Binary**: Support both Intel and Apple Silicon
3. **Code Signing**: Add developer signature for distribution
4. **DMG Packaging**: Create installer package
5. **Cross-Platform**: Extend to Linux and Windows

### Code Signing (Optional)
```bash
# For distribution, consider code signing
codesign --force --deep --sign "Developer ID Application: Your Name" build/bin/DESeq2.app
```

## Troubleshooting

### Common Issues
1. **Qt5 Static Path**: Ensure `/Users/venky/Softwares/qt5-static` exists
2. **Icon Issues**: Run `./create_icon.sh` to regenerate icon
3. **Permission Issues**: Ensure executable permissions on DESeq2 binary

### Verification Commands
```bash
# Check bundle structure
ls -la build/bin/DESeq2.app/Contents/

# Verify dependencies
otool -L build/bin/DESeq2.app/Contents/MacOS/DESeq2

# Test execution
./build/bin/DESeq2.app/Contents/MacOS/DESeq2
```

## Conclusion

The DESeq2 macOS bundle successfully combines:
- Static Qt5 libraries for portability
- Professional macOS application structure
- Self-contained deployment capability
- Full DESeq2 differential expression analysis functionality

This creates a robust, distributable application that maintains the benefits of static linking while providing a native macOS experience. 
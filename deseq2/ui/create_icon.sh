#!/bin/bash

# Create a basic DESeq2 icon using ImageMagick
# This creates a simple icon with "DESeq2" text

echo "Creating DESeq2 icon..."

# Check if ImageMagick is available
if ! command -v convert &>/dev/null; then
    echo "ImageMagick not found. Creating a placeholder icon..."
    # Create a simple text file as placeholder
    echo "DESeq2" >icon_placeholder.txt
    echo "Icon placeholder created. Please replace with proper icon later."
    exit 0
fi

# Create icon directory structure
mkdir -p DESeq2.iconset

# Create different sizes of the icon
sizes=(16 32 64 128 256 512 1024)

for size in "${sizes[@]}"; do
    echo "Creating ${size}x${size} icon..."

    # Create a simple icon with DESeq2 text
    convert -size ${size}x${size} \
        -background '#2E86AB' \
        -fill white \
        -gravity center \
        -pointsize $((size / 4)) \
        -font Arial-Bold \
        label:"DESeq2" \
        "DESeq2.iconset/icon_${size}x${size}.png"

    # Also create @2x versions for retina displays
    if [ $size -le 512 ]; then
        convert -size $((size * 2))x$((size * 2)) \
            -background '#2E86AB' \
            -fill white \
            -gravity center \
            -pointsize $((size / 2)) \
            -font Arial-Bold \
            label:"DESeq2" \
            "DESeq2.iconset/icon_${size}x${size}@2x.png"
    fi
done

# Create the .icns file
echo "Creating DESeq2.icns..."
iconutil -c icns DESeq2.iconset

# Clean up
rm -rf DESeq2.iconset

echo "Icon created: DESeq2.icns"

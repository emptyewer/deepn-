/bin/cp pythonqt/lib/*.dylib deepn/DEEPN++.app/Contents/MacOS

/bin/mkdir -p pythonqt/lib/GeneCount++.app/Contents/Scripts
/bin/cp -r ../gene_count/python/* pythonqt/lib/GeneCount++.app/Contents/Scripts

/bin/mkdir -p pythonqt/lib/JunctionMake++.app/Contents/Scripts
/bin/cp -r ../junction_make/python/* pythonqt/lib/JunctionMake++.app/Contents/Scripts

/bin/cp -r pythonqt/lib/GeneCount++.app deepn/DEEPN++.app/Contents/Resources
/bin/cp -r pythonqt/lib/JunctionMake++.app deepn/DEEPN++.app/Contents/Resources

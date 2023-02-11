# /bin/cp pythonqt/lib/*.dylib deepn/DEEPN++.app/Contents/MacOS

PYTHON_PATH=/usr/local/opt/python@3.11/Frameworks/Python.framework/Versions/3.11
SITEPACKAGE=/Users/vkrishnamani/.virtualenvs/deepn/lib/python3.11/site-packages

files=(GeneCount++ JunctionMake++)
for file in "${files[@]}"; do  # loop through the array
    # samtools
    /bin/mkdir -p $(pwd)/pythonqt/lib/${file}.app/Contents/Tools
    /bin/cp ../samtools/samtools $(pwd)/pythonqt/lib/${file}.app/Contents/Tools
    # scripts
    /bin/mkdir -p $(pwd)/pythonqt/lib/${file}.app/Contents/Scripts/data
    /bin/cp -r ../scripts/${file}/* $(pwd)/pythonqt/lib/${file}.app/Contents/Scripts
    # python
    rsync -avh --ignore-errors $PYTHON_PATH $(pwd)/pythonqt/lib/${file}.app/Contents/MacOS
    rsync -avh --ignore-errors $(pwd)/pythonqt/lib/*.dylib $(pwd)/pythonqt/lib/${file}.app/Contents/MacOS
    rsync -avh --ignore-errors $SITEPACKAGE/ $(pwd)/pythonqt/lib/${file}.app/Contents/MacOS/3.11/lib/python3.11
    # change lib loader paths
    install_name_tool -change libPythonQt-Qt5-Python3.11.3.dylib @executable_path/libPythonQt-Qt5-Python3.11.3.dylib $(pwd)/pythonqt/lib/${file}.app/Contents/MacOS/${file}
    install_name_tool -change libPythonQt_QtAll-Qt5-Python3.11.3.dylib @executable_path/libPythonQt_QtAll-Qt5-Python3.11.3.dylib $(pwd)/pythonqt/lib/${file}.app/Contents/MacOS/${file}
    install_name_tool -change $PYTHON_PATH/Python @executable_path/3.11/Python $(pwd)/pythonqt/lib/${file}.app/Contents/MacOS/${file}
done

rsync -avh --ignore-errors $(pwd)/pythonqt/lib/GeneCount++.app  $(pwd)/deepn/DEEPN++.app/Contents/Resources
rsync -avh --ignore-errors $(pwd)/pythonqt/lib/JunctionMake++.app  $(pwd)/deepn/DEEPN++.app/Contents/Resources
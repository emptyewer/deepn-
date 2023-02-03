CONFIG += sdk_no_version_check
TEMPLATE = subdirs
SUBDIRS = deepn gene_count junction_make pythonqt
QMAKE_MACOSX_DEPLOYMENT_TARGET = 13.0

gene_count.depends = pythonqt
junction_make.depends = pythonqt
deepn.depends = gene_count junction_make

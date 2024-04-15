CONFIG += sdk_no_version_check
TEMPLATE = subdirs
SUBDIRS = deepn \
    junction_dice \
    gene_count \
    pythonqt \
    query \
    read_depth \

QMAKE_MACOSX_DEPLOYMENT_TARGET = 13.0

gene_count.depends = pythonqt
junction_make.depends = pythonqt
deepn.depends = junction_dice gene_count query read_depth

CONFIG += sdk_no_version_check
TEMPLATE = subdirs
SUBDIRS = deepn \
    junction_dice \
    gene_count \
    query \
    read_depth \

QMAKE_MACOSX_DEPLOYMENT_TARGET = 14.0

deepn.depends = junction_dice gene_count query read_depth

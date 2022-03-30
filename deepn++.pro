CONFIG+=sdk_no_version_check
TEMPLATE = subdirs
SUBDIRS = deepn gene_count junction_make pythonqt
gene_count.depends = pythonqt
junction_make.depends = pythonqt
deepn.depends = gene_count junction_make

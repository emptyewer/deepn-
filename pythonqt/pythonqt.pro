TEMPLATE = subdirs

SUBDIRS = generator src extensions
extensions.depends += src

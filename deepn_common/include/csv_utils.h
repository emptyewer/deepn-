#ifndef DEEPN_CSV_UTILS_H
#define DEEPN_CSV_UTILS_H

#include "data_structures.h"
#include "gene_annotation_db.h"

#include <QChar>
#include <QString>
#include <QStringList>

namespace deepn {

// Detect whether a CSV header line uses tabs or commas
QChar detectSeparator(const QString& headerLine);

// Split a delimited line respecting quoted fields
QStringList splitDelimitedLine(const QString& line, QChar separator);

// Find the first matching column name (case-insensitive) in headers
int findColumnIndex(const QStringList& headers, const QStringList& names);

// Safely extract and clean a field from a split line
QString safeField(const QStringList& fields, int idx);

// Resolve a gene identifier to its preferred display name via the annotation DB
QString selectionKeyFor(const QString& gene, const GeneAnnotationDB& annotationDB);

// Parse a DESeq2 results CSV (handles both comma and tab delimiters, quoted fields)
QVector<DESeq2Result> parseDESeq2CSV(const QString& path);

}  // namespace deepn

#endif  // DEEPN_CSV_UTILS_H

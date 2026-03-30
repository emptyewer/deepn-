#ifndef DEEPN_SQLITE_JUNCTION_LOADER_H
#define DEEPN_SQLITE_JUNCTION_LOADER_H

#include "data_structures.h"

#include <QString>
#include <QStringList>
#include <QVector>

namespace deepn {

class SqliteJunctionLoader {
public:
    // Quick schema check used during workdir auto-discovery.
    static bool looksLikeJunctionDatabase(const QString& dbPath);

    // Open a depth database
    bool open(const QString& dbPath);
    void close();

    // Check if database has v2 schema (rstart/rend columns)
    bool hasPositionColumns() const;

    // Get schema version (0 if no metadata table, 1 if no rstart/rend, 2 if full)
    int schemaVersion() const;

    // Load all junctions for a specific gene
    // If gene annotation is provided, PPM is calculated
    GeneJunctionProfile loadGeneJunctions(const QString& gene,
                                           const GeneAnnotation& annotation = {}) const;

    // Collapse junctions by position
    static QVector<CollapsedJunction> collapseByPosition(const QVector<JunctionSite>& sites);

    // Get list of all genes in this database
    QStringList availableGenes() const;

    // Get total distinct reads in database
    int totalDistinctReads() const;

    // Get total reads for a specific gene
    int geneReadCount(const QString& gene) const;

    QString lastError() const;
    QString databasePath() const;

private:
    QString m_dbPath;
    QString m_connName;
    mutable QString m_error;
    mutable int m_totalReads = -1;  // cached
};

}  // namespace deepn

#endif  // DEEPN_SQLITE_JUNCTION_LOADER_H

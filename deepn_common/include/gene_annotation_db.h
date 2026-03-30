#ifndef DEEPN_GENE_ANNOTATION_DB_H
#define DEEPN_GENE_ANNOTATION_DB_H

#include "data_structures.h"

#include <QMap>
#include <QString>
#include <QStringList>
#include <QVector>

namespace deepn {

class GeneAnnotationDB {
public:
    // Load/create gene annotation cache from a FASTA file
    // On first call, parses FASTA and creates a SQLite cache at fastaPath + ".annotations.sqlite"
    // On subsequent calls, loads from cache (validates via SHA-256 checksum of FASTA)
    bool loadFromFasta(const QString& fastaPath);

    // Load from pre-built gene annotation SQLite (schema: table "gene" with
    // columns id, name, nm_number, start, stop, sequence)
    bool loadFromSqlite(const QString& sqlitePath);

    // Lookup by refseq accession (NM_*)
    GeneAnnotation findByRefseq(const QString& refseq) const;

    // Lookup by gene name
    GeneAnnotation findByGeneName(const QString& geneName) const;

    // Fetch sequence on demand (not loaded into memory during loadFromSqlite)
    QString sequenceForRefseq(const QString& refseq) const;

    // Fill in mRNALength and sequence for a gene loaded without them (from SQLite bulk load)
    void populateGeneDetails(GeneAnnotation& annotation) const;

    // Search by partial name (for autocomplete)
    QVector<GeneAnnotation> search(const QString& query, int maxResults = 20) const;

    // Get all gene names (for combo boxes)
    QStringList allGeneNames() const;
    QStringList allRefseqs() const;

    // Number of genes loaded
    int count() const;

    QString lastError() const;

private:
    bool parseFasta(const QString& fastaPath, const QString& cachePath);
    bool loadCache(const QString& cachePath);
    QString computeChecksum(const QString& filePath) const;

    QMap<QString, GeneAnnotation> m_byRefseq;
    QMap<QString, GeneAnnotation> m_byGeneName;
    QString m_sqlitePath;  // For on-demand sequence queries
    QString m_error;
};

}  // namespace deepn

#endif  // DEEPN_GENE_ANNOTATION_DB_H

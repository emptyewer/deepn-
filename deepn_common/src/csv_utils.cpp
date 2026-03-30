#include "csv_utils.h"

#include <QFile>
#include <QTextStream>

namespace deepn {

QChar detectSeparator(const QString& headerLine)
{
    return headerLine.contains('\t') ? QChar('\t') : QChar(',');
}

QStringList splitDelimitedLine(const QString& line, QChar separator)
{
    QStringList fields;
    QString current;
    bool inQuotes = false;

    for (int i = 0; i < line.size(); ++i) {
        const QChar ch = line.at(i);
        if (ch == '"') {
            if (inQuotes && i + 1 < line.size() && line.at(i + 1) == '"') {
                current += '"';
                ++i;
            } else {
                inQuotes = !inQuotes;
            }
        } else if (ch == separator && !inQuotes) {
            fields.append(current.trimmed());
            current.clear();
        } else {
            current += ch;
        }
    }

    fields.append(current.trimmed());
    for (QString& field : fields) {
        if (field.size() >= 2 && field.startsWith('"') && field.endsWith('"')) {
            field = field.mid(1, field.size() - 2);
        }
        field = field.trimmed();
    }
    return fields;
}

int findColumnIndex(const QStringList& headers, const QStringList& names)
{
    for (const QString& name : names) {
        for (int i = 0; i < headers.size(); ++i) {
            QString candidate = headers.at(i).trimmed();
            candidate.remove('"');
            if (candidate.compare(name, Qt::CaseInsensitive) == 0) {
                return i;
            }
        }
    }
    return -1;
}

QString safeField(const QStringList& fields, int idx)
{
    if (idx >= 0 && idx < fields.size()) {
        QString value = fields.at(idx).trimmed();
        value.remove('"');
        return value.trimmed();
    }
    return {};
}

QString selectionKeyFor(const QString& gene, const GeneAnnotationDB& annotationDB)
{
    const GeneAnnotation annotation = annotationDB.findByRefseq(gene);
    if (annotation.isValid() && !annotation.geneName.isEmpty()) {
        return annotation.geneName;
    }
    return gene.trimmed();
}

QVector<DESeq2Result> parseDESeq2CSV(const QString& path)
{
    QVector<DESeq2Result> results;

    QFile file(path);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return results;
    }

    QTextStream in(&file);
    const QString headerLine = in.readLine();
    if (headerLine.isEmpty()) {
        return results;
    }

    const QChar sep = detectSeparator(headerLine);
    const QStringList headers = splitDelimitedLine(headerLine, sep);
    int geneCol = findColumnIndex(headers, {"gene", "Gene", "gene_name", "geneName"});
    const int baseMeanCol = findColumnIndex(headers, {"baseMean", "basemean", "base_mean"});
    const int l2fcCol = findColumnIndex(headers, {"log2FoldChange", "log2foldchange", "l2fc", "log2FC"});
    const int lfcSECol = findColumnIndex(headers, {"lfcSE", "lfcse", "lfc_se"});
    const int statCol = findColumnIndex(headers, {"stat", "Stat"});
    const int pvalCol = findColumnIndex(headers, {"pvalue", "pval", "p_value"});
    const int padjCol = findColumnIndex(headers, {"padj", "p_adj", "adjusted_pvalue"});
    const int enrichCol = findColumnIndex(headers, {"enrichment", "Enrichment", "direction"});

    if (geneCol < 0) {
        geneCol = 0;
    }

    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        if (line.isEmpty()) continue;

        const QStringList fields = splitDelimitedLine(line, sep);
        if (fields.size() <= geneCol) continue;

        DESeq2Result r;
        r.gene = safeField(fields, geneCol);
        if (r.gene.isEmpty()) continue;
        if (baseMeanCol >= 0) r.baseMean = safeField(fields, baseMeanCol).toDouble();
        if (l2fcCol >= 0) r.log2FoldChange = safeField(fields, l2fcCol).toDouble();
        if (lfcSECol >= 0) r.lfcSE = safeField(fields, lfcSECol).toDouble();
        if (statCol >= 0) r.stat = safeField(fields, statCol).toDouble();
        if (pvalCol >= 0) r.pvalue = safeField(fields, pvalCol).toDouble();
        if (padjCol >= 0) r.padj = safeField(fields, padjCol).toDouble();
        if (enrichCol >= 0) r.enrichment = safeField(fields, enrichCol);

        results.append(r);
    }

    return results;
}

}  // namespace deepn

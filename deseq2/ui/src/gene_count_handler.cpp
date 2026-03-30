#include "gene_count_handler.h"
#include <QFile>
#include <QTextStream>
#include <QFileInfo>
#include <QDir>
#include <QHeaderView>
#include <QDebug>
#include <QRegularExpression>
#include <QSqlDatabase>
#include <QSqlQuery>
#include <QSqlError>
#include <algorithm>

namespace deseq2
{

    GeneCountHandler::GeneCountHandler()
    {
    }

    GeneCountData GeneCountHandler::parseGeneCountFile(const QString &filePath)
    {
        GeneCountData data;
        data.filePath = QFileInfo(filePath).absoluteFilePath();
        data.fileName = QFileInfo(filePath).fileName();
        data.sampleName = extractSampleName(data.fileName);
        data.isValid = false;

        QFile file(filePath);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            data.errorMessage = QString("Cannot open file: %1").arg(filePath);
            return data;
        }

        QTextStream in(&file);
        QStringList lines;
        while (!in.atEnd())
        {
            lines.append(in.readLine());
        }
        file.close();

        if (lines.isEmpty())
        {
            data.errorMessage = "File is empty";
            return data;
        }

        // Extract header information
        extractHeaderInfo(lines, data.totalReads, data.totalHits);

        // Find the start of data lines
        int dataStartLine = findDataStartLine(lines);
        if (dataStartLine == -1)
        {
            data.errorMessage = "Cannot find data section in file";
            return data;
        }

        // Parse data lines
        for (int i = dataStartLine; i < lines.size(); ++i)
        {
            QString line = lines[i].trimmed();
            if (line.isEmpty())
            {
                continue;
            }

            QStringList fields = parseCSVLine(line);
            if (fields.size() < 3)
            {
                continue; // Skip malformed lines
            }

            QString chromosome = fields[0].trimmed();
            QString geneName = fields[1].trimmed();
            QString ppmStr = fields[2].trimmed();

            // Skip header-like lines
            if (chromosome.toLower() == "chromosome" || geneName.toLower() == "genename")
            {
                continue;
            }

            // Parse PPM value
            bool ok;
            double ppm = ppmStr.toDouble(&ok);
            if (!ok)
            {
                continue; // Skip lines with invalid PPM values
            }

            // Store gene count data
            if (!geneName.isEmpty())
            {
                data.geneCounts[geneName] = ppm;
            }
        }

        data.isValid = true;
        return data;
    }

    GeneCountData GeneCountHandler::parseGeneCountXlsxFile(const QString &filePath)
    {
        // DEPRECATED: XLSX support removed (QXlsx dependency eliminated).
        // Use parseGeneCountSqliteFile() for GeneCount++ output instead.
        GeneCountData data;
        data.filePath = QFileInfo(filePath).absoluteFilePath();
        data.fileName = QFileInfo(filePath).fileName();
        data.sampleName = extractSampleName(data.fileName);
        data.isValid = false;
        data.totalReads = 0;
        data.totalHits = 0;
        data.errorMessage = QString("XLSX format is no longer supported. "
                                    "Please use the SQLite (.db) output from GeneCount++ instead: %1")
                                .arg(filePath);
        qWarning() << data.errorMessage;
        return data;
    }

    GeneCountData GeneCountHandler::parseGeneCountSqliteFile(const QString &filePath)
    {
        GeneCountData data;
        data.filePath = QFileInfo(filePath).absoluteFilePath();
        data.fileName = QFileInfo(filePath).fileName();
        data.sampleName = extractSampleName(data.fileName);
        data.isValid = false;
        data.totalReads = 0;
        data.totalHits = 0;

        // Use a unique connection name per file to allow concurrent opens
        QString connectionName = QString("genecount_%1").arg(filePath);

        {
            QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connectionName);
            db.setDatabaseName(filePath);

            if (!db.open())
            {
                data.errorMessage = QString("Cannot open SQLite database: %1 (%2)")
                                        .arg(filePath, db.lastError().text());
                QSqlDatabase::removeDatabase(connectionName);
                return data;
            }

            // Read summary table for metadata
            QSqlQuery summaryQuery(db);
            if (summaryQuery.exec("SELECT key, value FROM summary"))
            {
                while (summaryQuery.next())
                {
                    QString key = summaryQuery.value(0).toString();
                    QString value = summaryQuery.value(1).toString();

                    if (key == "total_reads")
                        data.totalReads = value.toInt();
                    else if (key == "total_hits")
                        data.totalHits = value.toInt();
                }
            }
            else
            {
                qDebug() << "No summary table or query failed:" << summaryQuery.lastError().text();
                // Not fatal -- we can still read gene counts
            }

            // Read gene_counts table
            QSqlQuery geneQuery(db);
            if (!geneQuery.exec("SELECT gene, count, ppm FROM gene_counts"))
            {
                data.errorMessage = QString("Failed to query gene_counts table: %1")
                                        .arg(geneQuery.lastError().text());
                db.close();
                QSqlDatabase::removeDatabase(connectionName);
                return data;
            }

            int totalCounts = 0;
            while (geneQuery.next())
            {
                QString geneName = geneQuery.value(0).toString().trimmed();
                int rawCount = geneQuery.value(1).toInt();
                double ppm = geneQuery.value(2).toDouble();

                if (geneName.isEmpty())
                    continue;

                data.geneRawCounts[geneName] = rawCount;
                data.geneCounts[geneName] = ppm;
                totalCounts += rawCount;
            }

            // If totalReads was not in summary, estimate from sum of counts
            if (data.totalReads == 0)
                data.totalReads = totalCounts;

            db.close();
        }

        // Remove the connection outside the scope where QSqlDatabase was in use
        QSqlDatabase::removeDatabase(connectionName);

        if (data.geneCounts.isEmpty())
        {
            data.errorMessage = "No gene data found in SQLite database";
            return data;
        }

        data.isValid = true;
        return data;
    }

    bool GeneCountHandler::addGeneCountFile(const QString &filePath)
    {
        // Detect file type and parse accordingly
        GeneCountData data;
        if (filePath.endsWith(".db", Qt::CaseInsensitive) ||
            filePath.endsWith(".sqlite", Qt::CaseInsensitive))
        {
            data = parseGeneCountSqliteFile(filePath);
        }
        else if (filePath.endsWith(".xlsx", Qt::CaseInsensitive))
        {
            data = parseGeneCountXlsxFile(filePath); // Returns invalid -- xlsx no longer supported
        }
        else
        {
            data = parseGeneCountFile(filePath);
        }
        if (!data.isValid)
        {
            qDebug() << "Failed to parse gene count file:" << data.errorMessage;
            return false;
        }

        // Check if file already exists
        for (const auto &existingData : m_geneCountFiles)
        {
            if (existingData.filePath == data.filePath)
            {
                qDebug() << "File already exists:" << data.filePath;
                return false;
            }
        }

        m_geneCountFiles.append(data);
        return true;
    }

    bool GeneCountHandler::removeGeneCountFile(const QString &fileName)
    {
        for (int i = 0; i < m_geneCountFiles.size(); ++i)
        {
            if (m_geneCountFiles[i].filePath == fileName || m_geneCountFiles[i].fileName == fileName)
            {
                m_geneCountFiles.removeAt(i);
                return true;
            }
        }
        return false;
    }

    void GeneCountHandler::clearAllFiles()
    {
        m_geneCountFiles.clear();
    }

    const QVector<GeneCountData> &GeneCountHandler::getGeneCountFiles() const
    {
        return m_geneCountFiles;
    }

    bool GeneCountHandler::assignGroupToFile(const QString &fileName, const QString &groupName)
    {
        for (int i = 0; i < m_geneCountFiles.size(); ++i)
        {
            auto &data = m_geneCountFiles[i];
            if (data.filePath == fileName || data.fileName == fileName)
            {
                if (data.groupName.isEmpty())
                {
                    // No group yet — assign directly
                    data.groupName = groupName;
                    return true;
                }
                if (data.groupName == groupName)
                {
                    // Already assigned to this group
                    return true;
                }
                // File already has a different group — check if a duplicate
                // with this group already exists
                bool alreadyDuplicated = false;
                for (const auto &other : m_geneCountFiles)
                {
                    if ((other.filePath == data.filePath) && other.groupName == groupName)
                    {
                        alreadyDuplicated = true;
                        break;
                    }
                }
                if (!alreadyDuplicated)
                {
                    // Duplicate the entry with the new group
                    GeneCountData dup = data;
                    dup.groupName = groupName;
                    dup.sampleName = data.sampleName + "_" + groupName;
                    m_geneCountFiles.append(dup);
                }
                return true;
            }
        }
        return false;
    }

    int GeneCountHandler::autoAssignGroups(const QMap<QString, QString> &patterns)
    {
        int assignedCount = 0;

        for (auto &data : m_geneCountFiles)
        {
            for (auto it = patterns.begin(); it != patterns.end(); ++it)
            {
                QString pattern = it.key();
                QString groupName = it.value();

                if (data.fileName.contains(pattern, Qt::CaseInsensitive))
                {
                    data.groupName = groupName;
                    assignedCount++;
                    break;
                }
            }
        }

        return assignedCount;
    }

    CombinedData GeneCountHandler::generateCountMatrixAndMetadata()
    {
        CombinedData result;
        result.isValid = false;

        if (m_geneCountFiles.isEmpty())
        {
            result.errorMessage = "No gene count files loaded";
            return result;
        }

        // Check if all files have group assignments
        for (const auto &data : m_geneCountFiles)
        {
            if (data.groupName.isEmpty())
            {
                result.errorMessage = QString("File '%1' has no group assignment").arg(data.fileName);
                return result;
            }
        }

        // Collect all unique gene names
        QSet<QString> uniqueGeneNames;
        for (const auto &data : m_geneCountFiles)
        {
            for (auto it = data.geneCounts.begin(); it != data.geneCounts.end(); ++it)
            {
                uniqueGeneNames.insert(it.key());
            }
        }

        result.geneNames = QVector<QString>(uniqueGeneNames.begin(), uniqueGeneNames.end());
        std::sort(result.geneNames.begin(), result.geneNames.end());

        // Collect sample names and group names
        for (const auto &data : m_geneCountFiles)
        {
            result.sampleNames.append(data.sampleName);
            result.groupNames.append(data.groupName);
        }

        // Build count matrix (PPM) and raw count matrix
        result.countMatrix.resize(result.geneNames.size());
        result.rawCountMatrix.resize(result.geneNames.size());
        for (int geneIdx = 0; geneIdx < result.geneNames.size(); ++geneIdx)
        {
            QString geneName = result.geneNames[geneIdx];
            result.countMatrix[geneIdx].resize(m_geneCountFiles.size());
            result.rawCountMatrix[geneIdx].resize(m_geneCountFiles.size());

            for (int sampleIdx = 0; sampleIdx < m_geneCountFiles.size(); ++sampleIdx)
            {
                const auto &data = m_geneCountFiles[sampleIdx];
                double ppmValue = data.geneCounts.value(geneName, 0.0);
                result.countMatrix[geneIdx][sampleIdx] = ppmValue;

                // Use raw counts if available (from xlsx), otherwise back-compute from PPM
                if (!data.geneRawCounts.isEmpty())
                {
                    result.rawCountMatrix[geneIdx][sampleIdx] = data.geneRawCounts.value(geneName, 0);
                }
                else
                {
                    // Back-compute: count = round(ppm * totalReads / 1e6)
                    int backComputed = static_cast<int>(std::round(ppmValue * data.totalReads / 1e6));
                    result.rawCountMatrix[geneIdx][sampleIdx] = std::max(0, backComputed);
                }
            }
        }

        result.isValid = true;
        return result;
    }

    bool GeneCountHandler::exportCountMatrixToCSV(const QString &filePath, const CombinedData &data)
    {
        if (!data.isValid)
        {
            return false;
        }

        QFile file(filePath);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            return false;
        }

        QTextStream out(&file);

        // Write header
        out << "Gene";
        for (const QString &sampleName : data.sampleNames)
        {
            out << "," << sampleName;
        }
        out << "\n";

        // Write data
        for (int geneIdx = 0; geneIdx < data.geneNames.size(); ++geneIdx)
        {
            out << data.geneNames[geneIdx];
            for (int sampleIdx = 0; sampleIdx < data.sampleNames.size(); ++sampleIdx)
            {
                out << "," << QString::number(data.countMatrix[geneIdx][sampleIdx], 'f', 6);
            }
            out << "\n";
        }

        file.close();
        return true;
    }

    bool GeneCountHandler::exportMetadataToCSV(const QString &filePath, const CombinedData &data)
    {
        if (!data.isValid)
        {
            return false;
        }

        QFile file(filePath);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            return false;
        }

        QTextStream out(&file);

        // Write header
        out << "Sample,Group\n";

        // Write data
        for (int i = 0; i < data.sampleNames.size(); ++i)
        {
            out << data.sampleNames[i] << "," << data.groupNames[i] << "\n";
        }

        file.close();
        return true;
    }

    QPair<bool, bool> GeneCountHandler::exportBothFiles(const QString &baseFilePath, const CombinedData &data)
    {
        if (!data.isValid)
        {
            return QPair<bool, bool>(false, false);
        }

        // Export count matrix
        QString countMatrixPath = baseFilePath + "_count_matrix.csv";
        bool countMatrixSuccess = exportCountMatrixToCSV(countMatrixPath, data);

        // Export metadata
        QString metadataPath = baseFilePath + "_metadata.csv";
        bool metadataSuccess = exportMetadataToCSV(metadataPath, data);

        return QPair<bool, bool>(countMatrixSuccess, metadataSuccess);
    }

    void GeneCountHandler::updateGeneCountFilesTable(QTableWidget *tableWidget)
    {
        if (!tableWidget)
        {
            return;
        }

        tableWidget->clear();
        tableWidget->setRowCount(m_geneCountFiles.size());
        tableWidget->setColumnCount(5);

        QStringList headers = {"File", "Sample Name", "Group", "Total Reads", "Total Hits"};
        tableWidget->setHorizontalHeaderLabels(headers);

        for (int row = 0; row < m_geneCountFiles.size(); ++row)
        {
            const auto &data = m_geneCountFiles[row];

            auto *fileItem = new QTableWidgetItem(data.fileName);
            fileItem->setToolTip(data.filePath);
            tableWidget->setItem(row, 0, fileItem);
            tableWidget->setItem(row, 1, new QTableWidgetItem(data.sampleName));
            tableWidget->setItem(row, 2, new QTableWidgetItem(data.groupName));
            tableWidget->setItem(row, 3, new QTableWidgetItem(QString::number(data.totalReads)));
            tableWidget->setItem(row, 4, new QTableWidgetItem(QString::number(data.totalHits)));
        }

        tableWidget->resizeColumnsToContents();
    }

    void GeneCountHandler::updateCountsPreviewTable(QTableWidget *tableWidget, const CombinedData &data)
    {
        if (!tableWidget || !data.isValid)
        {
            return;
        }

        const int maxRows = std::min(20, static_cast<int>(data.geneNames.size()));
        const int maxCols = std::min(10, static_cast<int>(data.sampleNames.size()));

        tableWidget->clear();
        tableWidget->setRowCount(maxRows);
        tableWidget->setColumnCount(maxCols + 1); // +1 for gene names

        QStringList headers = {"Gene"};
        for (int i = 0; i < maxCols; ++i)
        {
            headers.append(data.sampleNames[i]);
        }
        tableWidget->setHorizontalHeaderLabels(headers);

        for (int row = 0; row < maxRows; ++row)
        {
            tableWidget->setItem(row, 0, new QTableWidgetItem(data.geneNames[row]));
            for (int col = 0; col < maxCols; ++col)
            {
                double value = data.countMatrix[row][col];
                tableWidget->setItem(row, col + 1, new QTableWidgetItem(QString::number(value, 'f', 6)));
            }
        }

        tableWidget->resizeColumnsToContents();
    }

    void GeneCountHandler::updateMetadataPreviewTable(QTableWidget *tableWidget, const CombinedData &data)
    {
        if (!tableWidget || !data.isValid)
        {
            return;
        }

        tableWidget->clear();
        tableWidget->setRowCount(data.sampleNames.size());
        tableWidget->setColumnCount(2);

        QStringList headers = {"Sample", "Group"};
        tableWidget->setHorizontalHeaderLabels(headers);

        for (int row = 0; row < data.sampleNames.size(); ++row)
        {
            tableWidget->setItem(row, 0, new QTableWidgetItem(data.sampleNames[row]));
            tableWidget->setItem(row, 1, new QTableWidgetItem(data.groupNames[row]));
        }

        tableWidget->resizeColumnsToContents();
    }

    QString GeneCountHandler::extractSampleName(const QString &fileName)
    {
        // Remove file extension
        QString name = QFileInfo(fileName).baseName();

        // Remove common suffixes
        QStringList suffixes = {"_summary", "_counts", "_gene_counts"};
        for (const QString &suffix : suffixes)
        {
            if (name.endsWith(suffix, Qt::CaseInsensitive))
            {
                name = name.left(name.length() - suffix.length());
            }
        }

        return name;
    }

    bool GeneCountHandler::validateGeneCountData(const GeneCountData &data)
    {
        if (!data.isValid)
        {
            return false;
        }

        if (data.fileName.isEmpty())
        {
            return false;
        }

        if (data.filePath.isEmpty())
        {
            return false;
        }

        if (data.geneCounts.isEmpty())
        {
            return false;
        }

        return true;
    }

    QStringList GeneCountHandler::parseCSVLine(const QString &line)
    {
        QStringList fields;
        QString currentField;
        bool inQuotes = false;

        for (int i = 0; i < line.length(); ++i)
        {
            QChar ch = line[i];

            if (ch == '"')
            {
                inQuotes = !inQuotes;
            }
            else if (ch == ',' && !inQuotes)
            {
                fields.append(currentField.trimmed());
                currentField.clear();
            }
            else
            {
                currentField.append(ch);
            }
        }

        // Add the last field
        fields.append(currentField.trimmed());

        return fields;
    }

    int GeneCountHandler::findDataStartLine(const QStringList &lines)
    {
        for (int i = 0; i < lines.size(); ++i)
        {
            QString line = lines[i].trimmed();
            if (line.startsWith("Chromosome") && line.contains("GeneName"))
            {
                return i + 1; // Return the line after the header
            }
        }
        return -1;
    }

    void GeneCountHandler::extractHeaderInfo(const QStringList &lines, int &totalReads, int &totalHits)
    {
        totalReads = 0;
        totalHits = 0;

        for (const QString &line : lines)
        {
            QString trimmedLine = line.trimmed();
            if (trimmedLine.contains("TotalReads"))
            {
                QStringList parts = trimmedLine.split(",");
                if (parts.size() >= 3)
                {
                    totalReads = parts[2].trimmed().toInt();
                }
            }
            else if (trimmedLine.contains("TotalHits"))
            {
                QStringList parts = trimmedLine.split(",");
                if (parts.size() >= 3)
                {
                    totalHits = parts[2].trimmed().toInt();
                }
            }
        }
    }

} // namespace deseq2

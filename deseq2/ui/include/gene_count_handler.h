#ifndef GENE_COUNT_HANDLER_H
#define GENE_COUNT_HANDLER_H

#include <QString>
#include <QStringList>
#include <QMap>
#include <QVector>
#include <QTableWidget>
#include <memory>
#include <QPair>

namespace deseq2
{

    /**
     * @brief Structure to hold gene count data from a single file
     */
    struct GeneCountData
    {
        QString fileName;                 // Original filename
        QString sampleName;               // Extracted sample name
        QString groupName;                // Assigned group name
        QMap<QString, double> geneCounts; // Gene name -> PPM value mapping
        QMap<QString, int> geneRawCounts; // Gene name -> raw integer count
        int totalReads;                   // Total reads from file
        int totalHits;                    // Total hits from file
        bool isValid;                     // Whether the file was parsed successfully
        QString errorMessage;             // Error message if parsing failed
    };

    /**
     * @brief Structure to hold the combined count matrix and metadata
     */
    struct CombinedData
    {
        QVector<QString> geneNames;           // List of all unique gene names
        QVector<QString> sampleNames;         // List of sample names
        QVector<QString> groupNames;          // List of group names
        QVector<QVector<double>> countMatrix; // Count matrix (genes x samples) - PPM values
        QVector<QVector<int>> rawCountMatrix; // Raw integer count matrix (genes x samples)
        bool isValid;                         // Whether the data is valid
        QString errorMessage;                 // Error message if invalid
    };

    /**
     * @brief Class to handle gene count file operations
     *
     * This class provides functionality to:
     * - Parse gene count CSV files
     * - Extract gene names and PPM values
     * - Combine multiple files into a count matrix
     * - Generate metadata for DESeq2 analysis
     */
    class GeneCountHandler
    {
    public:
        /**
         * @brief Constructor
         */
        GeneCountHandler();

        /**
         * @brief Destructor
         */
        ~GeneCountHandler() = default;

        // Rule of Five - disable copy operations
        GeneCountHandler(const GeneCountHandler &) = delete;
        GeneCountHandler &operator=(const GeneCountHandler &) = delete;
        GeneCountHandler(GeneCountHandler &&) = delete;
        GeneCountHandler &operator=(GeneCountHandler &&) = delete;

        /**
         * @brief Parse a single gene count CSV file
         * @param filePath Path to the CSV file
         * @return GeneCountData structure with parsed data
         */
        GeneCountData parseGeneCountFile(const QString &filePath);

        /**
         * @brief Parse a single gene count XLSX file (from GeneCount++ output)
         * @param filePath Path to the XLSX file
         * @return GeneCountData structure with parsed data including raw counts
         */
        GeneCountData parseGeneCountXlsxFile(const QString &filePath);

        /**
         * @brief Add a gene count file to the collection
         * @param filePath Path to the CSV file
         * @return true if successfully added, false otherwise
         */
        bool addGeneCountFile(const QString &filePath);

        /**
         * @brief Remove a gene count file from the collection
         * @param fileName Name of the file to remove
         * @return true if successfully removed, false otherwise
         */
        bool removeGeneCountFile(const QString &fileName);

        /**
         * @brief Clear all gene count files
         */
        void clearAllFiles();

        /**
         * @brief Get the list of loaded gene count files
         * @return List of GeneCountData structures
         */
        const QVector<GeneCountData> &getGeneCountFiles() const;

        /**
         * @brief Assign a group name to a specific file
         * @param fileName Name of the file
         * @param groupName Group name to assign
         * @return true if successfully assigned, false otherwise
         */
        bool assignGroupToFile(const QString &fileName, const QString &groupName);

        /**
         * @brief Auto-assign groups based on filename patterns
         * @param patterns Map of filename patterns to group names
         * @return Number of files that were auto-assigned
         */
        int autoAssignGroups(const QMap<QString, QString> &patterns);

        /**
         * @brief Generate count matrix and metadata from loaded files
         * @return CombinedData structure with count matrix and metadata
         */
        CombinedData generateCountMatrixAndMetadata();

        /**
         * @brief Export count matrix to CSV file
         * @param filePath Output file path
         * @param data Combined data to export
         * @return true if successfully exported, false otherwise
         */
        bool exportCountMatrixToCSV(const QString &filePath, const CombinedData &data);

        /**
         * @brief Export metadata to CSV file
         * @param filePath Output file path
         * @param data Combined data to export
         * @return true if successfully exported, false otherwise
         */
        bool exportMetadataToCSV(const QString &filePath, const CombinedData &data);

        /**
         * @brief Export both count matrix and metadata to CSV files
         * @param baseFilePath Base file path (without extension)
         * @param data Combined data to export
         * @return Pair of booleans: (countMatrixSuccess, metadataSuccess)
         */
        QPair<bool, bool> exportBothFiles(const QString &baseFilePath, const CombinedData &data);

        /**
         * @brief Update the gene count files table in the UI
         * @param tableWidget Table widget to update
         */
        void updateGeneCountFilesTable(QTableWidget *tableWidget);

        /**
         * @brief Update the counts preview table in the UI
         * @param tableWidget Table widget to update
         * @param data Combined data to display
         */
        void updateCountsPreviewTable(QTableWidget *tableWidget, const CombinedData &data);

        /**
         * @brief Update the metadata preview table in the UI
         * @param tableWidget Table widget to update
         * @param data Combined data to display
         */
        void updateMetadataPreviewTable(QTableWidget *tableWidget, const CombinedData &data);

        /**
         * @brief Extract sample name from filename
         * @param fileName Full filename
         * @return Extracted sample name
         */
        QString extractSampleName(const QString &fileName);

        /**
         * @brief Validate gene count data
         * @param data GeneCountData to validate
         * @return true if valid, false otherwise
         */
        bool validateGeneCountData(const GeneCountData &data);

    private:
        QVector<GeneCountData> m_geneCountFiles; // Collection of parsed gene count files

        /**
         * @brief Parse CSV line and extract fields
         * @param line CSV line to parse
         * @return List of fields
         */
        QStringList parseCSVLine(const QString &line);

        /**
         * @brief Skip header lines in the CSV file
         * @param lines All lines from the file
         * @return Index of the first data line
         */
        int findDataStartLine(const QStringList &lines);

        /**
         * @brief Extract header information from CSV lines
         * @param lines All lines from the file
         * @param totalReads Output parameter for total reads
         * @param totalHits Output parameter for total hits
         */
        void extractHeaderInfo(const QStringList &lines, int &totalReads, int &totalHits);
    };

} // namespace deseq2

#endif // GENE_COUNT_HANDLER_H
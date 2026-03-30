#include "export_engine.h"

#include <QDebug>
#include <QFile>
#include <QFileInfo>
#include <QPainter>
#include <QPixmap>
#include <QPrinter>
#include <QSvgGenerator>
#include <QTextStream>
#include <QWidget>

namespace deepn {

static QString s_lastError;

bool ExportEngine::exportJunctionCSV(const QString& filePath, const GeneJunctionProfile& profile)
{
    s_lastError.clear();

    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        s_lastError = QStringLiteral("Cannot open file for writing: %1").arg(filePath);
        qDebug() << s_lastError;
        return false;
    }

    QTextStream out(&file);
    // Header
    out << "Gene,RefSeq,Position,PositionEnd,QueryStart,QueryEnd,Frame,CDS_Class,PPM,RawCount\n";

    for (const auto& site : profile.sites) {
        out << site.geneName << ","
            << site.refseq << ","
            << site.position << ","
            << site.positionEnd << ","
            << site.queryStart << ","
            << site.queryEnd << ","
            << site.frame << ","
            << site.cdsClass << ","
            << QString::number(site.ppm, 'f', 4) << ","
            << site.rawCount << "\n";
    }

    file.close();
    return true;
}

bool ExportEngine::exportCollapsedCSV(const QString& filePath, const QVector<CollapsedJunction>& collapsed)
{
    s_lastError.clear();

    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        s_lastError = QStringLiteral("Cannot open file for writing: %1").arg(filePath);
        qDebug() << s_lastError;
        return false;
    }

    QTextStream out(&file);
    out << "Position,TotalPPM,VariantCount,DominantFrame,CDS_Class\n";

    for (const auto& cj : collapsed) {
        out << cj.position << ","
            << QString::number(cj.totalPpm, 'f', 4) << ","
            << cj.variantCount << ","
            << cj.dominantFrame << ","
            << cj.cdsClass << "\n";
    }

    file.close();
    return true;
}

bool ExportEngine::exportDepthCSV(const QString& filePath, const DepthProfile& profile)
{
    s_lastError.clear();

    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        s_lastError = QStringLiteral("Cannot open file for writing: %1").arg(filePath);
        qDebug() << s_lastError;
        return false;
    }

    QTextStream out(&file);
    out << "Position,Count,Normalized\n";

    for (const auto& pt : profile.points) {
        out << pt.position << ","
            << pt.count << ","
            << QString::number(pt.normalized, 'f', 4) << "\n";
    }

    file.close();
    return true;
}

bool ExportEngine::exportBatchSummaryCSV(const QString& filePath,
                                          const QVector<GeneJunctionProfile>& profiles)
{
    s_lastError.clear();

    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        s_lastError = QStringLiteral("Cannot open file for writing: %1").arg(filePath);
        qDebug() << s_lastError;
        return false;
    }

    QTextStream out(&file);
    out << "Gene,RefSeq,ORFStart,ORFEnd,mRNALength,TotalSites,TotalReads,SourceFile\n";

    for (const auto& profile : profiles) {
        out << profile.annotation.geneName << ","
            << profile.annotation.refseq << ","
            << profile.annotation.orfStart << ","
            << profile.annotation.orfEnd << ","
            << profile.annotation.mRNALength << ","
            << profile.sites.size() << ","
            << profile.totalReads << ","
            << QFileInfo(profile.sourceFile).fileName() << "\n";
    }

    file.close();
    return true;
}

bool ExportEngine::exportFigureSVG(const QString& filePath, QWidget* widget, QSize size)
{
    s_lastError.clear();

    if (!widget) {
        s_lastError = "Widget is null";
        qDebug() << s_lastError;
        return false;
    }

    QSvgGenerator generator;
    generator.setFileName(filePath);
    generator.setSize(size);
    generator.setViewBox(QRect(0, 0, size.width(), size.height()));
    generator.setTitle(QStringLiteral("DEEPN++ Export"));

    QPainter painter;
    if (!painter.begin(&generator)) {
        s_lastError = QStringLiteral("Cannot begin SVG painter for: %1").arg(filePath);
        qDebug() << s_lastError;
        return false;
    }

    widget->render(&painter);
    painter.end();

    return true;
}

bool ExportEngine::exportFigurePDF(const QString& filePath, QWidget* widget, QSize size)
{
    s_lastError.clear();

    if (!widget) {
        s_lastError = "Widget is null";
        qDebug() << s_lastError;
        return false;
    }

    QPrinter printer(QPrinter::HighResolution);
    printer.setOutputFormat(QPrinter::PdfFormat);
    printer.setOutputFileName(filePath);
    printer.setPageSize(QPageSize(QSizeF(size.width(), size.height()), QPageSize::Point));
    printer.setPageMargins(QMarginsF(0, 0, 0, 0));

    QPainter painter;
    if (!painter.begin(&printer)) {
        s_lastError = QStringLiteral("Cannot begin PDF painter for: %1").arg(filePath);
        qDebug() << s_lastError;
        return false;
    }

    // Scale widget rendering to fit the printer page
    double xScale = printer.pageRect(QPrinter::DevicePixel).width() / static_cast<double>(size.width());
    double yScale = printer.pageRect(QPrinter::DevicePixel).height() / static_cast<double>(size.height());
    double scale = qMin(xScale, yScale);
    painter.scale(scale, scale);

    widget->render(&painter);
    painter.end();

    return true;
}

bool ExportEngine::exportFigurePNG(const QString& filePath, QWidget* widget, QSize size, int dpi)
{
    s_lastError.clear();

    if (!widget) {
        s_lastError = "Widget is null";
        qDebug() << s_lastError;
        return false;
    }

    // Scale size by DPI ratio (base 72 DPI for screen)
    double scaleFactor = dpi / 72.0;
    QSize renderSize(static_cast<int>(size.width() * scaleFactor),
                     static_cast<int>(size.height() * scaleFactor));

    QPixmap pixmap(renderSize);
    pixmap.fill(Qt::white);

    QPainter painter(&pixmap);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setRenderHint(QPainter::TextAntialiasing, true);
    painter.scale(scaleFactor, scaleFactor);
    widget->render(&painter);
    painter.end();

    if (!pixmap.save(filePath, "PNG")) {
        s_lastError = QStringLiteral("Failed to save PNG: %1").arg(filePath);
        qDebug() << s_lastError;
        return false;
    }

    return true;
}

QString ExportEngine::lastError()
{
    return s_lastError;
}

}  // namespace deepn

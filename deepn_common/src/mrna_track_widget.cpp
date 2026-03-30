#include "mrna_track_widget.h"

#include <QDebug>
#include <QFontMetrics>
#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>

#include <cmath>

namespace deepn {

MRNATrackWidget::MRNATrackWidget(QWidget* parent)
    : QWidget(parent)
{
    setMinimumHeight(50);
    setMaximumHeight(70);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    setMouseTracking(true);
}

void MRNATrackWidget::setAnnotation(const GeneAnnotation& annotation)
{
    m_annotation = annotation;
    m_visibleStart = 0;
    m_visibleEnd = annotation.mRNALength;
    update();
}

void MRNATrackWidget::setVisibleRange(int start, int end)
{
    m_visibleStart = start;
    m_visibleEnd = end;
    update();
}

int MRNATrackWidget::positionToX(int pos) const
{
    int range = m_visibleEnd - m_visibleStart;
    if (range <= 0)
        return 0;

    const int margin = 10;
    int drawWidth = width() - 2 * margin;
    return margin + static_cast<int>((static_cast<double>(pos - m_visibleStart) / range) * drawWidth);
}

int MRNATrackWidget::xToPosition(int x) const
{
    const int margin = 10;
    int drawWidth = width() - 2 * margin;
    if (drawWidth <= 0)
        return m_visibleStart;

    double fraction = static_cast<double>(x - margin) / drawWidth;
    fraction = qBound(0.0, fraction, 1.0);
    int range = m_visibleEnd - m_visibleStart;
    return m_visibleStart + static_cast<int>(fraction * range);
}

void MRNATrackWidget::paintEvent(QPaintEvent* event)
{
    Q_UNUSED(event)

    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, true);

    int w = width();
    int h = height();
    int margin = 10;

    // Background
    painter.fillRect(rect(), Qt::white);

    if (!m_annotation.isValid() || m_visibleEnd <= m_visibleStart) {
        painter.setPen(Qt::gray);
        painter.drawText(rect(), Qt::AlignCenter, "No gene annotation loaded");
        return;
    }

    // Track layout: gene track in upper portion, scale bar below
    int trackY = 12;
    int trackH = 14;
    int utrH = 6;
    int utrY = trackY + (trackH - utrH) / 2;
    int scaleY = trackY + trackH + 6;

    // Draw full-length thin line (mRNA backbone)
    painter.setPen(QPen(QColor(CDS_COLOR), 1));
    int x0 = positionToX(qMax(0, m_visibleStart));
    int x1 = positionToX(qMin(m_annotation.mRNALength, m_visibleEnd));
    painter.drawLine(x0, trackY + trackH / 2, x1, trackY + trackH / 2);

    // Draw 5' UTR region (thin, light grey)
    if (m_annotation.orfStart > m_visibleStart) {
        int utrStart = positionToX(qMax(m_visibleStart, 0));
        int utrEnd = positionToX(qMin(m_annotation.orfStart, m_visibleEnd));
        if (utrEnd > utrStart) {
            painter.fillRect(utrStart, utrY, utrEnd - utrStart, utrH, QColor(UTR_COLOR));
        }
    }

    // Draw CDS region (thick, dark grey)
    int cdsDrawStart = qMax(m_annotation.orfStart, m_visibleStart);
    int cdsDrawEnd = qMin(m_annotation.orfEnd, m_visibleEnd);
    if (cdsDrawEnd > cdsDrawStart) {
        int cdsX0 = positionToX(cdsDrawStart);
        int cdsX1 = positionToX(cdsDrawEnd);
        painter.fillRect(cdsX0, trackY, cdsX1 - cdsX0, trackH, QColor(CDS_COLOR));
    }

    // Draw 3' UTR region (thin, light grey)
    if (m_annotation.orfEnd < m_visibleEnd) {
        int utr3Start = positionToX(qMax(m_annotation.orfEnd, m_visibleStart));
        int utr3End = positionToX(qMin(m_annotation.mRNALength, m_visibleEnd));
        if (utr3End > utr3Start) {
            painter.fillRect(utr3Start, utrY, utr3End - utr3Start, utrH, QColor(UTR_COLOR));
        }
    }

    // Draw start codon marker (green triangle pointing down at ATG)
    if (m_annotation.orfStart >= m_visibleStart && m_annotation.orfStart <= m_visibleEnd) {
        int atgX = positionToX(m_annotation.orfStart);
        QPolygonF triangle;
        triangle << QPointF(atgX, trackY - 1)
                 << QPointF(atgX - 4, trackY - 7)
                 << QPointF(atgX + 4, trackY - 7);
        painter.setPen(Qt::NoPen);
        painter.setBrush(QColor(START_CODON_COLOR));
        painter.drawPolygon(triangle);

        QFont font("Helvetica", 7, QFont::Bold);
        painter.setFont(font);
        painter.setPen(QColor(START_CODON_COLOR));
        QFontMetrics fm(font);
        int labelW = fm.horizontalAdvance("ATG");
        painter.drawText(atgX - labelW / 2, trackY - 8, "ATG");
    }

    // Draw stop codon marker (red triangle pointing down at stop)
    if (m_annotation.orfEnd >= m_visibleStart && m_annotation.orfEnd <= m_visibleEnd) {
        int stopX = positionToX(m_annotation.orfEnd);
        QPolygonF triangle;
        triangle << QPointF(stopX, trackY - 1)
                 << QPointF(stopX - 4, trackY - 7)
                 << QPointF(stopX + 4, trackY - 7);
        painter.setPen(Qt::NoPen);
        painter.setBrush(QColor(STOP_CODON_COLOR));
        painter.drawPolygon(triangle);

        QFont font("Helvetica", 7, QFont::Bold);
        painter.setFont(font);
        painter.setPen(QColor(STOP_CODON_COLOR));
        QFontMetrics fm(font);
        int labelW = fm.horizontalAdvance("STOP");
        painter.drawText(stopX - labelW / 2, trackY - 8, "STOP");
    }

    // CDS label (inside the dark band if wide enough)
    if (cdsDrawEnd > cdsDrawStart) {
        int cdsPixelWidth = positionToX(cdsDrawEnd) - positionToX(cdsDrawStart);
        if (cdsPixelWidth > 40) {
            QFont font("Helvetica", 7);
            painter.setFont(font);
            painter.setPen(Qt::white);
            int labelX = (positionToX(cdsDrawStart) + positionToX(cdsDrawEnd)) / 2;
            QFontMetrics fm(font);
            int labelW = fm.horizontalAdvance("CDS");
            painter.drawText(labelX - labelW / 2, trackY + trackH / 2 + 4, "CDS");
        }
    }

    // Draw scale bar with position ticks
    painter.setPen(QPen(Qt::black, 1));
    painter.drawLine(margin, scaleY, w - margin, scaleY);

    // Determine tick interval based on visible range
    int visibleRange = m_visibleEnd - m_visibleStart;
    int tickInterval;
    if (visibleRange > 10000)
        tickInterval = 2000;
    else if (visibleRange > 5000)
        tickInterval = 1000;
    else if (visibleRange > 2000)
        tickInterval = 500;
    else if (visibleRange > 1000)
        tickInterval = 200;
    else if (visibleRange > 500)
        tickInterval = 100;
    else
        tickInterval = 50;

    QFont tickFont("Helvetica", 7);
    painter.setFont(tickFont);
    painter.setPen(Qt::black);
    QFontMetrics fm(tickFont);

    int firstTick = ((m_visibleStart + tickInterval - 1) / tickInterval) * tickInterval;
    if (firstTick < m_visibleStart)
        firstTick += tickInterval;

    for (int pos = firstTick; pos <= m_visibleEnd; pos += tickInterval) {
        int x = positionToX(pos);
        if (x < margin || x > w - margin)
            continue;

        painter.drawLine(x, scaleY - 2, x, scaleY + 2);

        QString label = QString::number(pos);
        int labelW = fm.horizontalAdvance(label);
        int labelX = x - labelW / 2;
        labelX = qMax(margin, qMin(labelX, w - margin - labelW));
        painter.drawText(labelX, scaleY + 12, label);
    }

    // Gene name + info at bottom
    QFont infoFont("Helvetica", 8, QFont::Bold);
    painter.setFont(infoFont);
    painter.setPen(Qt::black);
    QString geneInfo = QStringLiteral("%1 (%2) - %3 bp")
                           .arg(m_annotation.geneName, m_annotation.refseq)
                           .arg(m_annotation.mRNALength);
    painter.drawText(margin, h - 3, geneInfo);
}

void MRNATrackWidget::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton && m_annotation.isValid()) {
        int pos = xToPosition(event->position().toPoint().x());
        pos = qBound(0, pos, m_annotation.mRNALength);
        emit positionClicked(pos);
    }
    QWidget::mousePressEvent(event);
}

}  // namespace deepn

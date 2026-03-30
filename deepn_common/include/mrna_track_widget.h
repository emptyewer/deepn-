#ifndef DEEPN_MRNA_TRACK_WIDGET_H
#define DEEPN_MRNA_TRACK_WIDGET_H

#include "data_structures.h"

#include <QWidget>

namespace deepn {

class MRNATrackWidget : public QWidget {
    Q_OBJECT
public:
    explicit MRNATrackWidget(QWidget* parent = nullptr);

    void setAnnotation(const GeneAnnotation& annotation);
    void setVisibleRange(int start, int end);  // linked to chart zoom

    // Colors
    static constexpr const char* CDS_COLOR = "#374151";         // dark grey
    static constexpr const char* UTR_COLOR = "#D1D5DB";         // light grey
    static constexpr const char* START_CODON_COLOR = "#059669"; // green
    static constexpr const char* STOP_CODON_COLOR = "#DC2626";  // red

signals:
    void positionClicked(int position);

protected:
    void paintEvent(QPaintEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;

private:
    GeneAnnotation m_annotation;
    int m_visibleStart = 0;
    int m_visibleEnd = 0;

    int positionToX(int pos) const;
    int xToPosition(int x) const;
};

}  // namespace deepn

#endif  // DEEPN_MRNA_TRACK_WIDGET_H

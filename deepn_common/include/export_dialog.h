#ifndef DEEPN_EXPORT_DIALOG_H
#define DEEPN_EXPORT_DIALOG_H

#include <QDialog>
#include <QSize>
#include <QString>

class QComboBox;
class QDoubleSpinBox;
class QSpinBox;

namespace deepn {

struct ExportSettings {
    enum Format { PNG, SVG, PDF };
    Format format = PNG;
    int dpi = 300;
    int width = 800;
    int height = 400;

    QString filter() const;
    QString defaultExtension() const;
};

class ExportDialog : public QDialog {
    Q_OBJECT
public:
    explicit ExportDialog(QWidget* parent = nullptr);

    ExportSettings settings() const;
    void setDefaults(const ExportSettings& defaults);

private:
    QComboBox* m_formatCombo;
    QSpinBox* m_dpiSpin;
    QSpinBox* m_widthSpin;
    QSpinBox* m_heightSpin;
};

}  // namespace deepn

#endif  // DEEPN_EXPORT_DIALOG_H

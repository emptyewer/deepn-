#include "export_dialog.h"

#include <QComboBox>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QGroupBox>
#include <QLabel>
#include <QSpinBox>
#include <QVBoxLayout>

namespace deepn {

QString ExportSettings::filter() const
{
    switch (format) {
    case PNG: return "PNG Images (*.png)";
    case SVG: return "SVG Files (*.svg)";
    case PDF: return "PDF Files (*.pdf)";
    }
    return "All Files (*)";
}

QString ExportSettings::defaultExtension() const
{
    switch (format) {
    case PNG: return ".png";
    case SVG: return ".svg";
    case PDF: return ".pdf";
    }
    return ".png";
}

ExportDialog::ExportDialog(QWidget* parent)
    : QDialog(parent)
{
    setWindowTitle("Export Figure");
    setMinimumWidth(320);

    auto* layout = new QVBoxLayout(this);

    // Format group
    auto* formatGroup = new QGroupBox("Format", this);
    auto* formatLayout = new QFormLayout(formatGroup);

    m_formatCombo = new QComboBox(this);
    m_formatCombo->addItem("PNG (raster)", static_cast<int>(ExportSettings::PNG));
    m_formatCombo->addItem("SVG (vector)", static_cast<int>(ExportSettings::SVG));
    m_formatCombo->addItem("PDF (vector)", static_cast<int>(ExportSettings::PDF));
    formatLayout->addRow("Format:", m_formatCombo);

    m_dpiSpin = new QSpinBox(this);
    m_dpiSpin->setRange(72, 1200);
    m_dpiSpin->setValue(300);
    m_dpiSpin->setSuffix(" DPI");
    m_dpiSpin->setSingleStep(50);
    formatLayout->addRow("Resolution:", m_dpiSpin);

    layout->addWidget(formatGroup);

    // Size group
    auto* sizeGroup = new QGroupBox("Dimensions", this);
    auto* sizeLayout = new QFormLayout(sizeGroup);

    m_widthSpin = new QSpinBox(this);
    m_widthSpin->setRange(200, 4000);
    m_widthSpin->setValue(800);
    m_widthSpin->setSuffix(" px");
    m_widthSpin->setSingleStep(100);
    sizeLayout->addRow("Width:", m_widthSpin);

    m_heightSpin = new QSpinBox(this);
    m_heightSpin->setRange(100, 3000);
    m_heightSpin->setValue(400);
    m_heightSpin->setSuffix(" px");
    m_heightSpin->setSingleStep(50);
    sizeLayout->addRow("Height:", m_heightSpin);

    layout->addWidget(sizeGroup);

    // Enable/disable DPI based on format
    connect(m_formatCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, [this](int) {
                auto fmt = static_cast<ExportSettings::Format>(
                    m_formatCombo->currentData().toInt());
                m_dpiSpin->setEnabled(fmt == ExportSettings::PNG);
            });

    // Buttons
    auto* buttons = new QDialogButtonBox(
        QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
    connect(buttons, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
    layout->addWidget(buttons);
}

ExportSettings ExportDialog::settings() const
{
    ExportSettings s;
    s.format = static_cast<ExportSettings::Format>(
        m_formatCombo->currentData().toInt());
    s.dpi = m_dpiSpin->value();
    s.width = m_widthSpin->value();
    s.height = m_heightSpin->value();
    return s;
}

void ExportDialog::setDefaults(const ExportSettings& defaults)
{
    for (int i = 0; i < m_formatCombo->count(); ++i) {
        if (m_formatCombo->itemData(i).toInt() == static_cast<int>(defaults.format)) {
            m_formatCombo->setCurrentIndex(i);
            break;
        }
    }
    m_dpiSpin->setValue(defaults.dpi);
    m_widthSpin->setValue(defaults.width);
    m_heightSpin->setValue(defaults.height);
}

}  // namespace deepn

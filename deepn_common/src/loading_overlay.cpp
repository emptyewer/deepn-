#include "loading_overlay.h"

#include <QApplication>
#include <QEvent>
#include <QScrollBar>

namespace deepn {

LoadingOverlay::LoadingOverlay(QWidget* parent)
    : QWidget(parent)
{
    // 90% opaque black background covering the entire window
    setAttribute(Qt::WA_StyledBackground, true);
    setAttribute(Qt::WA_TransparentForMouseEvents, false);
    setAutoFillBackground(true);
    setStyleSheet("background-color: rgba(0, 0, 0, 230);");

    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    layout->addStretch(1);

    m_titleLabel = new QLabel(this);
    m_titleLabel->setAlignment(Qt::AlignCenter);
    m_titleLabel->setStyleSheet(
        "color: white; font-size: 22px; font-weight: bold; "
        "background: transparent; padding: 20px;");
    layout->addWidget(m_titleLabel);

    m_log = new QPlainTextEdit(this);
    m_log->setReadOnly(true);
    m_log->setMinimumHeight(180);
    m_log->setMaximumHeight(250);
    m_log->setMinimumWidth(500);
    m_log->setMaximumWidth(600);
    m_log->setStyleSheet(
        "color: #8EE89E; background: rgba(0, 0, 0, 180); "
        "border: 1px solid #555; border-radius: 8px; "
        "font-family: monospace; font-size: 13px; padding: 10px;");
    layout->addWidget(m_log, 0, Qt::AlignCenter);

    layout->addStretch(1);

    setVisible(false);

    if (parent) {
        parent->installEventFilter(this);
    }
}

void LoadingOverlay::show(const QString& title)
{
    m_titleLabel->setText(title);
    m_log->clear();
    if (parentWidget()) {
        setGeometry(0, 0, parentWidget()->width(), parentWidget()->height());
    }
    raise();
    setVisible(true);
    QApplication::processEvents();
}

void LoadingOverlay::hide()
{
    setVisible(false);
    QApplication::processEvents();
}

void LoadingOverlay::addMessage(const QString& message)
{
    m_log->appendPlainText(message);
    QScrollBar* sb = m_log->verticalScrollBar();
    sb->setValue(sb->maximum());
    QApplication::processEvents();
}

bool LoadingOverlay::eventFilter(QObject* obj, QEvent* event)
{
    if (obj == parentWidget() && event->type() == QEvent::Resize) {
        setGeometry(0, 0, parentWidget()->width(), parentWidget()->height());
    }
    return QWidget::eventFilter(obj, event);
}

}  // namespace deepn

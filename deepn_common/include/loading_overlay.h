#ifndef DEEPN_LOADING_OVERLAY_H
#define DEEPN_LOADING_OVERLAY_H

#include <QLabel>
#include <QPlainTextEdit>
#include <QVBoxLayout>
#include <QWidget>

namespace deepn {

class LoadingOverlay : public QWidget {
    Q_OBJECT
public:
    explicit LoadingOverlay(QWidget* parent);

    void show(const QString& title = "Loading...");
    void hide();
    void addMessage(const QString& message);

protected:
    bool eventFilter(QObject* obj, QEvent* event) override;

private:
    QLabel* m_titleLabel;
    QPlainTextEdit* m_log;
};

}  // namespace deepn

#endif  // DEEPN_LOADING_OVERLAY_H

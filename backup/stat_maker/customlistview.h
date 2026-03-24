#ifndef CUSTOMLISTVIEW_H
#define CUSTOMLISTVIEW_H

#include <QListWidget>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QMimeData>
#include <QStringListModel>
#include "signals.h"

class CustomListView : public QListWidget {
Q_OBJECT

public:
    explicit CustomListView(QWidget *parent = nullptr);

    ~CustomListView() override;

protected:
    void dragEnterEvent(QDragEnterEvent *event) override;

    void dragMoveEvent(QDragMoveEvent *event) override;

    void dropEvent(QDropEvent *event) override;

    void startDrag(Qt::DropActions supportedActions) override;

    void keyPressEvent(QKeyEvent *event) override;

private:
    QStringListModel *model;
    Signals *sig = Signals::getCommonInstance();

    bool isFileAlreadyInList(const QString &filePath) const;


};

#endif // CUSTOMLISTVIEW_H
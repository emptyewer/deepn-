#include "customlistview.h"
#include <QUrl>
#include <QDebug>
#include <QFileIconProvider>
#include <QFileInfo>
#include <QMimeData>
#include <QDrag>

CustomListView::CustomListView(QWidget *parent)
        : QListWidget(parent) {
    setAcceptDrops(true);
    setDragEnabled(true);
}

void CustomListView::startDrag(Qt::DropActions supportedActions) {
    QList<QListWidgetItem *> items = selectedItems();
    if (items.empty()) {
        return;
    }

    auto *mimeData = new QMimeData;
    QList<QUrl> urlList;
    for (QListWidgetItem *item: items) {
        urlList << item->data(Qt::UserRole).toString();
    }
    mimeData->setUrls(urlList);

    auto *drag = new QDrag(this);
    drag->setMimeData(mimeData);
    drag->exec(supportedActions, Qt::CopyAction);

}

void CustomListView::dragEnterEvent(QDragEnterEvent *event) {
    if (event->mimeData()->hasUrls()) {
        event->acceptProposedAction();
    }
}

void CustomListView::dragMoveEvent(QDragMoveEvent *event) {
    if (event->mimeData()->hasUrls()) {
        event->acceptProposedAction();
    }
}

void CustomListView::dropEvent(QDropEvent *event) {
    if (event->mimeData()->hasUrls()) {
        QList<QUrl> urlList = event->mimeData()->urls();
        QFileIconProvider iconProvider;

        for (const QUrl &url: urlList) {
            QString filepath = url.toLocalFile();
            QFileInfo file(filepath);
            if (!isFileAlreadyInList(filepath)) {
                auto *item = new QListWidgetItem(file.baseName());
                item->setData(Qt::UserRole, url);
                item->setIcon(iconProvider.icon(QFileInfo(url.toLocalFile())));
                addItem(item);
            } else {
                qDebug() << "File already in list";
                emit sig->displayStatus(file.baseName() + " already in list");

            }
        }
        event->acceptProposedAction();
    }
}

bool CustomListView::isFileAlreadyInList(const QString &filePath) const {
    for (int i = 0; i < count(); ++i) {
        QListWidgetItem *item = this->item(i);
        if (item->data(Qt::UserRole).toUrl().toLocalFile() == filePath) {
            return true;
        }
    }
    return false;
}

void CustomListView::keyPressEvent(QKeyEvent *event) {
    if (this->objectName() == "file_list") {
        return;
    }
    if (event->key() == Qt::Key_Delete || event->key() == Qt::Key_Backspace) {
        qDeleteAll(selectedItems());
        event->accept();
    } else {
        QListWidget::keyPressEvent(event);
    }
}

CustomListView::~CustomListView() = default;

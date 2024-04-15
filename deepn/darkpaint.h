#ifndef DARKPAINT_H
#define DARKPAINT_H

#include <QPainter>
#include <QWidget>

class Overlay : public QWidget {
 public:
  Overlay(QWidget *parent) {
    setPalette(Qt::transparent);
    setAttribute(Qt::WA_TransparentForMouseEvents);
  }

 protected:
  void paintEvent(QPaintEvent *event) {
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setBrush(QBrush(QColor(0, 0, 0, 100)));
    painter.setPen(Qt::NoPen);
    painter.setOpacity(0.5);
    painter.drawRect(rect());
  }
};

#endif  // DARKPAINT_H

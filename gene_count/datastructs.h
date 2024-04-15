#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H

#include <QMap>
#include <QString>
#include <algorithm>

struct GCStat {
  QString input;
  int readCount = 0;
  QString readName = "";
  int elapsedTime;
  bool running = true;
};

#endif // DATASTRUCTS_H

#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H

#include <QString>

struct GCStat {
  QString filename;
  int readCount;
  QString currentRefSeq;
  int elapsedTime;
  bool running = true;
};

#endif // DATASTRUCTS_H

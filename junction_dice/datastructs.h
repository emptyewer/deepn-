#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H

#include <QMap>
#include <QString>
#include <algorithm>

enum MapAlgo { blat, blastn, megablast };

struct MapStat {
  MapAlgo algo = blat;
  QString input;
  QString output;
  QString db;
  QString currentRead;
  int elapsedTime;
  float percentComplete;
  QString map_command;
};

struct DiceStat {
  QString input;
  QString output;
  QString summary;
  QString junction;
  int matchLength;
  QString match;
  QString currentRead = "";
  long long fileSize;
  int elapsedTime;
  int totalReads = 0;
  int writtenReads = 0;
  int removedReads = 0;
  int forwardMatches = 0;
  int reverseMatches = 0;
};

struct JDStat {
  MapStat mstat;
  DiceStat dstat;
  bool readdice = false;
  bool dicing = true;
  bool blasting = false;
};

#endif  // DATASTRUCTS_H

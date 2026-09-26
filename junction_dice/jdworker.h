#ifndef JDWORKER_H
#define JDWORKER_H

#include <QDebug>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QHash>
#include <QList>
#include <QMap>
#include <QObject>
#include <QPair>
#include <QProcess>
#include <QRegularExpression>
#include <QString>
#include <QStringList>
#include <QTextStream>
#include <QThread>
#include <QTimer>
#include <QtSql>

#include "datastructs.h"
#include "signals.h"

class JDWorker : public QObject {
  Q_OBJECT

 public:
  explicit JDWorker(JDStat *stat, int fileCount);
  ~JDWorker();

 signals:
  void finished();

 public slots:
  void run();

 private slots:
  void mapCallBack();

 private:
  int fileCount;
  QString jseq_pattern;   // Tier 1: full matchLength bases (primary search)
  QString jseq_pattern2;  // Tier 2: 2/3 of matchLength (recovers 5' junction errors)
  QString jseq_pattern3;  // Tier 3: 1/3 of matchLength (recovers most junction errors)
  QString repeats_sequence;
  QElapsedTimer elapsedTimer;
  QTimer *mapTimer;
  Signals *sig = Signals::getCommonInstance();
  JDStat *stat;
  QProcess process;
  QHash<QString, int> readIndex;
  int mapPollCount = 0;
  bool mapFinished = false;
  QString readDepthFileName;
  QString dbConnectionName;
  QSqlDatabase db;
  QSqlQuery query;
  void doDice();
  void doMapping();
  QString reverseComplement(QString dna_sequence);
  QString translate(QString dna);
  void readDice();
  void writeDiceSummary();
  void readDiceSummary();
  void createDepthDatabase();
};

#endif // JDWORKER_H

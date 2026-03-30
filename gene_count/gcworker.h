#ifndef GCWORKER_H
#define GCWORKER_H

#include <QElapsedTimer>
#include <QObject>
#include <QThread>
#include <QtSql>

#include "datastructs.h"
#include "maphits.h"
#include "signals.h"

class GCWorker : public QObject {
  Q_OBJECT

 public:
  explicit GCWorker(GCStat *stat);

 public slots:
  void run();

 signals:
  void finished();

 private:
  QElapsedTimer elapsedTimer;
  QSqlDatabase db;
  QSqlQuery query;
  Signals *sig = Signals::getCommonInstance();
  GCStat *stat;
  QString mappedOuputDBName;
  QString writeDbConn;
  void setupDB();
  void writeReadHitsToDB(ReadHits& hits);
  void writeToDatabase(QList<ReadHits> *collectedReads);
  void readMapOutput();
  void writeGeneCount();
};

#endif // GCWORKER_H

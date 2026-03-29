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
  QSqlDatabase mem_db;
  QSqlQuery query;
  Signals *sig = Signals::getCommonInstance();
  GCStat *stat;
  QString mappedOuputDBName;
  QString writeDbConn;
  QString fileDbConn;
  QString memDbConn;
  void setupDB();
  void writeToDatabase(QList<ReadHits> *collectedReads);
  void readMapOutput();
  void writeGeneCount();
  void createInMemoryDB();
};

#endif // GCWORKER_H

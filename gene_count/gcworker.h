#ifndef GCWORKER_H
#define GCWORKER_H

#include <QElapsedTimer>
#include <QObject>

#include "signals.h"

class GCWorker : public QObject {
  Q_OBJECT
 public:
  explicit GCWorker(QString file);

 public slots:
  void doWork();
  void finished();
  void stop();

 signals:

 private:
  QElapsedTimer timer;
  QString sam_file;
  Signals *sig = Signals::getCommonInstance();
};

#endif // GCWORKER_H

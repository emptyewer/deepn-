#ifndef SIGNALS_H
#define SIGNALS_H

#include <QObject>

#include "datastructs.h"

class Signals : public QObject {
  Q_OBJECT
 public:
  static Signals *getCommonInstance();

 signals:
  void gc_update_progress_sig(GCStat stat);
  void gc_finished_sig();

 private:
  static Signals *inst_;
  explicit Signals(QObject *parent = nullptr);
  Signals(const Signals &);
};

#endif  // SIGNALS_H

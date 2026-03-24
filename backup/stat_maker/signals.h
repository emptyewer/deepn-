#ifndef SIGNALS_H
#define SIGNALS_H

#include <QObject>

class Signals : public QObject {
Q_OBJECT
public:
    static Signals *getCommonInstance();

signals:

    void displayStatus(QString status);


private:
    static Signals *inst_;

    explicit Signals(QObject *parent = nullptr);

    Signals(const Signals &);
};

#endif  // SIGNALS_H
#include "signals.h"

Signals *Signals::inst_ = NULL;

Signals *Signals::getCommonInstance() {
    if (inst_ == NULL) {
        inst_ = new Signals();
    }
    return (inst_);
}

Signals::Signals(QObject *parent) : QObject(parent) {}
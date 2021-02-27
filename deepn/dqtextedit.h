#ifndef DQTEXTEDIT_H
#define DQTEXTEDIT_H

#include <QPlainTextEdit>

class DQTextEdit : public QPlainTextEdit
{
public:
    DQTextEdit(QWidget* parent);
    ~DQTextEdit() {}
};

#endif // DQTEXTEDIT_H

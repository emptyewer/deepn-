/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.15.13
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPlainTextEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralwidget;
    QVBoxLayout *verticalLayout;
    QPlainTextEdit *jd_output;
    QGridLayout *gridLayout;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QSpinBox *jseq_match_len;
    QPushButton *dice_btn;
    QGroupBox *blastChoices;
    QHBoxLayout *horizontalLayout_4;
    QRadioButton *blat;
    QRadioButton *megablast;
    QRadioButton *blastn;
    QPushButton *exit_btn;
    QButtonGroup *blastGroup;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(708, 791);
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(MainWindow->sizePolicy().hasHeightForWidth());
        MainWindow->setSizePolicy(sizePolicy);
        MainWindow->setMinimumSize(QSize(700, 600));
        QFont font;
        font.setFamily(QString::fromUtf8("Courier New"));
        font.setBold(true);
        MainWindow->setFont(font);
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        verticalLayout = new QVBoxLayout(centralwidget);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        jd_output = new QPlainTextEdit(centralwidget);
        jd_output->setObjectName(QString::fromUtf8("jd_output"));
        jd_output->setReadOnly(true);

        verticalLayout->addWidget(jd_output);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        label = new QLabel(centralwidget);
        label->setObjectName(QString::fromUtf8("label"));

        horizontalLayout->addWidget(label);

        jseq_match_len = new QSpinBox(centralwidget);
        jseq_match_len->setObjectName(QString::fromUtf8("jseq_match_len"));
        jseq_match_len->setMinimum(21);
        jseq_match_len->setMaximum(48);
        jseq_match_len->setSingleStep(3);

        horizontalLayout->addWidget(jseq_match_len);


        gridLayout->addLayout(horizontalLayout, 0, 1, 1, 1);

        dice_btn = new QPushButton(centralwidget);
        dice_btn->setObjectName(QString::fromUtf8("dice_btn"));

        gridLayout->addWidget(dice_btn, 0, 2, 2, 1);

        blastChoices = new QGroupBox(centralwidget);
        blastChoices->setObjectName(QString::fromUtf8("blastChoices"));
        horizontalLayout_4 = new QHBoxLayout(blastChoices);
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        blat = new QRadioButton(blastChoices);
        blastGroup = new QButtonGroup(MainWindow);
        blastGroup->setObjectName(QString::fromUtf8("blastGroup"));
        blastGroup->addButton(blat);
        blat->setObjectName(QString::fromUtf8("blat"));
        blat->setEnabled(true);
        blat->setChecked(false);

        horizontalLayout_4->addWidget(blat);

        megablast = new QRadioButton(blastChoices);
        blastGroup->addButton(megablast);
        megablast->setObjectName(QString::fromUtf8("megablast"));
        megablast->setEnabled(true);

        horizontalLayout_4->addWidget(megablast);

        blastn = new QRadioButton(blastChoices);
        blastGroup->addButton(blastn);
        blastn->setObjectName(QString::fromUtf8("blastn"));
        blastn->setChecked(true);

        horizontalLayout_4->addWidget(blastn);


        gridLayout->addWidget(blastChoices, 1, 1, 1, 1);

        exit_btn = new QPushButton(centralwidget);
        exit_btn->setObjectName(QString::fromUtf8("exit_btn"));

        gridLayout->addWidget(exit_btn, 0, 0, 2, 1);


        verticalLayout->addLayout(gridLayout);

        MainWindow->setCentralWidget(centralwidget);

        retranslateUi(MainWindow);
        QObject::connect(exit_btn, SIGNAL(clicked()), MainWindow, SLOT(close()));

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "Junction Dice++", nullptr));
        label->setText(QCoreApplication::translate("MainWindow", "Junction Sequence Match Length", nullptr));
        dice_btn->setText(QCoreApplication::translate("MainWindow", "Dice and Map", nullptr));
        blastChoices->setTitle(QString());
        blat->setText(QCoreApplication::translate("MainWindow", "BLAT", nullptr));
        megablast->setText(QCoreApplication::translate("MainWindow", "MegaBlast", nullptr));
        blastn->setText(QCoreApplication::translate("MainWindow", "BlastN", nullptr));
        exit_btn->setText(QCoreApplication::translate("MainWindow", "Exit", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H

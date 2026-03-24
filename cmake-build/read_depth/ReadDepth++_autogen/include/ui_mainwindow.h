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
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include "../QCustomplot/qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralwidget;
    QVBoxLayout *verticalLayout_3;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QLineEdit *search_text;
    QFrame *line;
    QLabel *label_2;
    QSpinBox *depth_count_interval;
    QSpacerItem *horizontalSpacer_2;
    QCheckBox *plot_individual;
    QFrame *line_2;
    QPushButton *save;
    QFrame *line_3;
    QHBoxLayout *horizontalLayout_2;
    QVBoxLayout *verticalLayout_2;
    QLabel *label_3;
    QComboBox *file_list;
    QFrame *line_4;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout;
    QTableWidget *depth_table;
    QCustomPlot *plot_view;
    QHBoxLayout *refseq_buttons;
    QTextEdit *sequence_view;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1330, 762);
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        centralwidget->setEnabled(false);
        verticalLayout_3 = new QVBoxLayout(centralwidget);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        label = new QLabel(centralwidget);
        label->setObjectName(QString::fromUtf8("label"));
        QFont font;
        font.setFamily(QString::fromUtf8("Menlo"));
        font.setPointSize(13);
        font.setBold(true);
        label->setFont(font);

        horizontalLayout->addWidget(label);

        search_text = new QLineEdit(centralwidget);
        search_text->setObjectName(QString::fromUtf8("search_text"));
        QFont font1;
        font1.setFamily(QString::fromUtf8("Menlo"));
        search_text->setFont(font1);

        horizontalLayout->addWidget(search_text);

        line = new QFrame(centralwidget);
        line->setObjectName(QString::fromUtf8("line"));
        line->setFrameShape(QFrame::VLine);
        line->setFrameShadow(QFrame::Sunken);

        horizontalLayout->addWidget(line);

        label_2 = new QLabel(centralwidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setFont(font);

        horizontalLayout->addWidget(label_2);

        depth_count_interval = new QSpinBox(centralwidget);
        depth_count_interval->setObjectName(QString::fromUtf8("depth_count_interval"));
        depth_count_interval->setFont(font1);
        depth_count_interval->setMinimum(5);
        depth_count_interval->setMaximum(200);
        depth_count_interval->setSingleStep(10);
        depth_count_interval->setValue(20);

        horizontalLayout->addWidget(depth_count_interval);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_2);

        plot_individual = new QCheckBox(centralwidget);
        plot_individual->setObjectName(QString::fromUtf8("plot_individual"));
        QFont font2;
        font2.setFamily(QString::fromUtf8("Menlo"));
        font2.setBold(true);
        plot_individual->setFont(font2);

        horizontalLayout->addWidget(plot_individual);

        line_2 = new QFrame(centralwidget);
        line_2->setObjectName(QString::fromUtf8("line_2"));
        line_2->setFrameShape(QFrame::VLine);
        line_2->setFrameShadow(QFrame::Sunken);

        horizontalLayout->addWidget(line_2);

        save = new QPushButton(centralwidget);
        save->setObjectName(QString::fromUtf8("save"));
        save->setFont(font1);

        horizontalLayout->addWidget(save);


        verticalLayout_3->addLayout(horizontalLayout);

        line_3 = new QFrame(centralwidget);
        line_3->setObjectName(QString::fromUtf8("line_3"));
        line_3->setFrameShape(QFrame::HLine);
        line_3->setFrameShadow(QFrame::Sunken);

        verticalLayout_3->addWidget(line_3);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        label_3 = new QLabel(centralwidget);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setFont(font);

        verticalLayout_2->addWidget(label_3);

        file_list = new QComboBox(centralwidget);
        file_list->setObjectName(QString::fromUtf8("file_list"));
        file_list->setFont(font1);

        verticalLayout_2->addWidget(file_list);

        line_4 = new QFrame(centralwidget);
        line_4->setObjectName(QString::fromUtf8("line_4"));
        line_4->setFrameShape(QFrame::HLine);
        line_4->setFrameShadow(QFrame::Sunken);

        verticalLayout_2->addWidget(line_4);

        groupBox = new QGroupBox(centralwidget);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        groupBox->setMinimumSize(QSize(320, 0));
        groupBox->setMaximumSize(QSize(640, 16777215));
        verticalLayout = new QVBoxLayout(groupBox);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        depth_table = new QTableWidget(groupBox);
        depth_table->setObjectName(QString::fromUtf8("depth_table"));
        depth_table->setMaximumSize(QSize(640, 16777215));
        depth_table->setFont(font1);
        depth_table->setEditTriggers(QAbstractItemView::NoEditTriggers);
        depth_table->setAlternatingRowColors(true);
        depth_table->setSelectionMode(QAbstractItemView::SingleSelection);
        depth_table->setSelectionBehavior(QAbstractItemView::SelectRows);
        depth_table->setSortingEnabled(true);
        depth_table->horizontalHeader()->setStretchLastSection(true);

        verticalLayout->addWidget(depth_table);


        verticalLayout_2->addWidget(groupBox);


        horizontalLayout_2->addLayout(verticalLayout_2);

        plot_view = new QCustomPlot(centralwidget);
        plot_view->setObjectName(QString::fromUtf8("plot_view"));
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(plot_view->sizePolicy().hasHeightForWidth());
        plot_view->setSizePolicy(sizePolicy);
        plot_view->setStyleSheet(QString::fromUtf8("QWidget {\n"
"    border: 2px solid black;\n"
"    border-radius: 5px;\n"
"}"));

        horizontalLayout_2->addWidget(plot_view);


        verticalLayout_3->addLayout(horizontalLayout_2);

        refseq_buttons = new QHBoxLayout();
        refseq_buttons->setObjectName(QString::fromUtf8("refseq_buttons"));

        verticalLayout_3->addLayout(refseq_buttons);

        sequence_view = new QTextEdit(centralwidget);
        sequence_view->setObjectName(QString::fromUtf8("sequence_view"));
        sequence_view->setMaximumSize(QSize(16777215, 400));
        sequence_view->setFont(font1);
        sequence_view->setReadOnly(true);

        verticalLayout_3->addWidget(sequence_view);

        MainWindow->setCentralWidget(centralwidget);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "Read Depth++", nullptr));
        label->setText(QCoreApplication::translate("MainWindow", "RefSeq Accession Number / Gene Name  ", nullptr));
        search_text->setPlaceholderText(QCoreApplication::translate("MainWindow", "Search by Gene or NM number", nullptr));
        label_2->setText(QCoreApplication::translate("MainWindow", "Depth Count Interval  ", nullptr));
        plot_individual->setText(QCoreApplication::translate("MainWindow", "Plot Individual", nullptr));
        save->setText(QCoreApplication::translate("MainWindow", "Save to CSV", nullptr));
        label_3->setText(QCoreApplication::translate("MainWindow", "Select Analyzed File", nullptr));
        groupBox->setTitle(QCoreApplication::translate("MainWindow", "Results", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H

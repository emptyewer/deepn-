/********************************************************************************
** Form generated from reading UI file 'dockwidget.ui'
**
** Created by: Qt User Interface Compiler version 5.15.13
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DOCKWIDGET_H
#define UI_DOCKWIDGET_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include "../QCustomplot/qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_DockWidget
{
public:
    QWidget *dockWidgetContents;
    QVBoxLayout *verticalLayout_4;
    QVBoxLayout *verticalLayout;
    QLabel *label;
    QHBoxLayout *horizontalLayout;
    QLineEdit *search_text;
    QPushButton *save;
    QVBoxLayout *verticalLayout_3;
    QLabel *label_2;
    QComboBox *dataset_selection;
    QFrame *line;
    QSplitter *splitter;
    QTabWidget *tabWidget;
    QWidget *summary;
    QGridLayout *gridLayout_3;
    QTableWidget *summary_table;
    QWidget *filtered_summary;
    QGridLayout *gridLayout_2;
    QTableWidget *filtered_summary_table;
    QWidget *plot;
    QGridLayout *gridLayout;
    QCustomPlot *plot_view;
    QTabWidget *refseq_info;
    QWidget *gene_tab;
    QGridLayout *gridLayout_4;
    QTextEdit *gene_info;
    QWidget *protein_tab;
    QVBoxLayout *verticalLayout_2;
    QTextEdit *protein_info;

    void setupUi(QDockWidget *DockWidget)
    {
        if (DockWidget->objectName().isEmpty())
            DockWidget->setObjectName(QString::fromUtf8("DockWidget"));
        DockWidget->resize(571, 761);
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(DockWidget->sizePolicy().hasHeightForWidth());
        DockWidget->setSizePolicy(sizePolicy);
        QFont font;
        font.setFamily(QString::fromUtf8("Courier New"));
        DockWidget->setFont(font);
        dockWidgetContents = new QWidget();
        dockWidgetContents->setObjectName(QString::fromUtf8("dockWidgetContents"));
        verticalLayout_4 = new QVBoxLayout(dockWidgetContents);
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        label = new QLabel(dockWidgetContents);
        label->setObjectName(QString::fromUtf8("label"));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Minimum);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy1);
        QFont font1;
        font1.setFamily(QString::fromUtf8("Courier New"));
        font1.setPointSize(13);
        font1.setBold(true);
        label->setFont(font1);

        verticalLayout->addWidget(label);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        search_text = new QLineEdit(dockWidgetContents);
        search_text->setObjectName(QString::fromUtf8("search_text"));
        search_text->setFont(font);

        horizontalLayout->addWidget(search_text);

        save = new QPushButton(dockWidgetContents);
        save->setObjectName(QString::fromUtf8("save"));
        QFont font2;
        font2.setFamily(QString::fromUtf8("Courier New"));
        font2.setBold(true);
        save->setFont(font2);

        horizontalLayout->addWidget(save);


        verticalLayout->addLayout(horizontalLayout);


        verticalLayout_4->addLayout(verticalLayout);

        verticalLayout_3 = new QVBoxLayout();
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        label_2 = new QLabel(dockWidgetContents);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        sizePolicy1.setHeightForWidth(label_2->sizePolicy().hasHeightForWidth());
        label_2->setSizePolicy(sizePolicy1);
        label_2->setFont(font1);

        verticalLayout_3->addWidget(label_2);

        dataset_selection = new QComboBox(dockWidgetContents);
        dataset_selection->setObjectName(QString::fromUtf8("dataset_selection"));
        dataset_selection->setFont(font);

        verticalLayout_3->addWidget(dataset_selection);


        verticalLayout_4->addLayout(verticalLayout_3);

        line = new QFrame(dockWidgetContents);
        line->setObjectName(QString::fromUtf8("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout_4->addWidget(line);

        splitter = new QSplitter(dockWidgetContents);
        splitter->setObjectName(QString::fromUtf8("splitter"));
        splitter->setLineWidth(1);
        splitter->setMidLineWidth(5);
        splitter->setOrientation(Qt::Vertical);
        splitter->setHandleWidth(10);
        tabWidget = new QTabWidget(splitter);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tabWidget->setFont(font);
        summary = new QWidget();
        summary->setObjectName(QString::fromUtf8("summary"));
        gridLayout_3 = new QGridLayout(summary);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        summary_table = new QTableWidget(summary);
        if (summary_table->columnCount() < 5)
            summary_table->setColumnCount(5);
        QTableWidgetItem *__qtablewidgetitem = new QTableWidgetItem();
        summary_table->setHorizontalHeaderItem(0, __qtablewidgetitem);
        QTableWidgetItem *__qtablewidgetitem1 = new QTableWidgetItem();
        summary_table->setHorizontalHeaderItem(1, __qtablewidgetitem1);
        QTableWidgetItem *__qtablewidgetitem2 = new QTableWidgetItem();
        summary_table->setHorizontalHeaderItem(2, __qtablewidgetitem2);
        QTableWidgetItem *__qtablewidgetitem3 = new QTableWidgetItem();
        summary_table->setHorizontalHeaderItem(3, __qtablewidgetitem3);
        QTableWidgetItem *__qtablewidgetitem4 = new QTableWidgetItem();
        summary_table->setHorizontalHeaderItem(4, __qtablewidgetitem4);
        summary_table->setObjectName(QString::fromUtf8("summary_table"));
        summary_table->setSortingEnabled(true);
        summary_table->horizontalHeader()->setStretchLastSection(true);
        summary_table->verticalHeader()->setVisible(false);
        summary_table->verticalHeader()->setStretchLastSection(true);

        gridLayout_3->addWidget(summary_table, 0, 0, 1, 1);

        tabWidget->addTab(summary, QString());
        filtered_summary = new QWidget();
        filtered_summary->setObjectName(QString::fromUtf8("filtered_summary"));
        gridLayout_2 = new QGridLayout(filtered_summary);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        filtered_summary_table = new QTableWidget(filtered_summary);
        filtered_summary_table->setObjectName(QString::fromUtf8("filtered_summary_table"));

        gridLayout_2->addWidget(filtered_summary_table, 0, 0, 1, 1);

        tabWidget->addTab(filtered_summary, QString());
        plot = new QWidget();
        plot->setObjectName(QString::fromUtf8("plot"));
        gridLayout = new QGridLayout(plot);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        plot_view = new QCustomPlot(plot);
        plot_view->setObjectName(QString::fromUtf8("plot_view"));
        QSizePolicy sizePolicy2(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(plot_view->sizePolicy().hasHeightForWidth());
        plot_view->setSizePolicy(sizePolicy2);

        gridLayout->addWidget(plot_view, 0, 0, 1, 1);

        tabWidget->addTab(plot, QString());
        splitter->addWidget(tabWidget);
        refseq_info = new QTabWidget(splitter);
        refseq_info->setObjectName(QString::fromUtf8("refseq_info"));
        gene_tab = new QWidget();
        gene_tab->setObjectName(QString::fromUtf8("gene_tab"));
        gridLayout_4 = new QGridLayout(gene_tab);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        gene_info = new QTextEdit(gene_tab);
        gene_info->setObjectName(QString::fromUtf8("gene_info"));
        QSizePolicy sizePolicy3(QSizePolicy::Expanding, QSizePolicy::Preferred);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(gene_info->sizePolicy().hasHeightForWidth());
        gene_info->setSizePolicy(sizePolicy3);
        gene_info->setFont(font);
        gene_info->setReadOnly(true);

        gridLayout_4->addWidget(gene_info, 0, 0, 1, 1);

        refseq_info->addTab(gene_tab, QString());
        protein_tab = new QWidget();
        protein_tab->setObjectName(QString::fromUtf8("protein_tab"));
        verticalLayout_2 = new QVBoxLayout(protein_tab);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        protein_info = new QTextEdit(protein_tab);
        protein_info->setObjectName(QString::fromUtf8("protein_info"));
        protein_info->setReadOnly(true);

        verticalLayout_2->addWidget(protein_info);

        refseq_info->addTab(protein_tab, QString());
        splitter->addWidget(refseq_info);

        verticalLayout_4->addWidget(splitter);

        DockWidget->setWidget(dockWidgetContents);

        retranslateUi(DockWidget);

        tabWidget->setCurrentIndex(0);
        refseq_info->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(DockWidget);
    } // setupUi

    void retranslateUi(QDockWidget *DockWidget)
    {
        DockWidget->setWindowTitle(QCoreApplication::translate("DockWidget", "Query", nullptr));
        label->setText(QCoreApplication::translate("DockWidget", "RefSeq Accession Number", nullptr));
        save->setText(QCoreApplication::translate("DockWidget", "Save to CSV", nullptr));
        label_2->setText(QCoreApplication::translate("DockWidget", "Select Data", nullptr));
        QTableWidgetItem *___qtablewidgetitem = summary_table->horizontalHeaderItem(0);
        ___qtablewidgetitem->setText(QCoreApplication::translate("DockWidget", "Position", nullptr));
        QTableWidgetItem *___qtablewidgetitem1 = summary_table->horizontalHeaderItem(1);
        ___qtablewidgetitem1->setText(QCoreApplication::translate("DockWidget", "Junctions #", nullptr));
        QTableWidgetItem *___qtablewidgetitem2 = summary_table->horizontalHeaderItem(2);
        ___qtablewidgetitem2->setText(QCoreApplication::translate("DockWidget", "Query Start", nullptr));
        QTableWidgetItem *___qtablewidgetitem3 = summary_table->horizontalHeaderItem(3);
        ___qtablewidgetitem3->setText(QCoreApplication::translate("DockWidget", "CDS", nullptr));
        QTableWidgetItem *___qtablewidgetitem4 = summary_table->horizontalHeaderItem(4);
        ___qtablewidgetitem4->setText(QCoreApplication::translate("DockWidget", "Frame", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(summary), QCoreApplication::translate("DockWidget", "Summary", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(filtered_summary), QCoreApplication::translate("DockWidget", "Filtered Summary", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(plot), QCoreApplication::translate("DockWidget", "Plot", nullptr));
        refseq_info->setTabText(refseq_info->indexOf(gene_tab), QCoreApplication::translate("DockWidget", "Gene", nullptr));
        refseq_info->setTabText(refseq_info->indexOf(protein_tab), QCoreApplication::translate("DockWidget", "Protein", nullptr));
    } // retranslateUi

};

namespace Ui {
    class DockWidget: public Ui_DockWidget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DOCKWIDGET_H

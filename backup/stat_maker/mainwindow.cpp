#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QString>
#include <QMessageBox>
#include <QMutex>
#include <QFileInfo>
#include <QFileIconProvider>
#include <Security/Authorization.h>
#include <Security/AuthorizationTags.h>
#include <unistd.h>
#include <datastructs.h>

bool DEBUG = true;
QMutex mutex;

MainWindow::MainWindow(int argc, char *argv[], QWidget *parent)
        : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    qInstallMessageHandler(MainWindow::customMessageHandler);
    if (DEBUG) {
        files << "/Volumes/Vault/Download.Vault/gc_files/HDPTP_Tail_2_NON_50M_summary.csv"
              << "/Volumes/Vault/Download.Vault/gc_files/HDPTP_Tail_2_SEL_50N_summary.csv"
              << "/Volumes/Vault/Download.Vault/gc_files/HDPTP_V_NON_51M_summary.csv"
              << "/Volumes/Vault/Download.Vault/gc_files/HDPTP_V_SEL_51N_summary.csv"
              << "/Volumes/Vault/Download.Vault/gc_files/Tef_Vector1_Non_50A_summary.csv"
              << "/Volumes/Vault/Download.Vault/gc_files/Tef_Vector1_SEL_50B_summary.csv"
              << "/Volumes/Vault/Download.Vault/gc_files/Tef_Vector2_NON_51A_summary.csv"
              << "/Volumes/Vault/Download.Vault/gc_files/Tef_Vector2_SEL_51B_summary.csv";
    } else {
        for (int i = 1; i < argc; i++) {
            files << argv[i];
        }
    }
    setupSlots();
    loadFiles();
    checkRFramework();
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::setupSlots() {
    connect(sig, &Signals::displayStatus, this, &MainWindow::displayStatus);
}


void MainWindow::loadFiles() {
    QFileIconProvider iconProvider;
    for (const QString &file: files) {
        QFileInfo fileInfo(file);
        auto *item = new QListWidgetItem(fileInfo.baseName());
        item->setIcon(iconProvider.icon(fileInfo));
        item->setData(Qt::UserRole, QUrl::fromLocalFile(file));
        ui->file_list->addItem(item);
    }
}

void MainWindow::customMessageHandler(QtMsgType type, const QMessageLogContext &context, const QString &msg) {
    mutex.lock();
    // Redirect to instance method
    MainWindow *instance = nullptr;
    auto widgets = QApplication::topLevelWidgets();
    for (QWidget *widget: widgets) {
        instance = dynamic_cast<MainWindow *>(widget);
        if (instance) {
            MainWindow::handleMessage(type, context, msg);
            break;
        }
    }
    mutex.unlock();
}

void MainWindow::handleMessage(QtMsgType type, const QMessageLogContext &context, const QString &msg) {
    QString txt;
    switch (type) {
        case QtDebugMsg:
            txt = QString("Debug: %1").arg(msg);
            break;
        case QtWarningMsg:
            txt = QString("Warning: %1").arg(msg);
            break;
        case QtCriticalMsg:
            txt = QString("Critical: %1").arg(msg);
            break;
        case QtFatalMsg:
            txt = QString("Fatal: %1").arg(msg);
            abort();
        case QtInfoMsg:
            break;
    }
    if (DEBUG) {
        qDebug() << txt;
    }
}

void MainWindow::displayStatus(const QString &status) {
    ui->statusbar->showMessage(status, 1000);
}

bool MainWindow::checkRFramework() {
    QStringList paths = {
            "/Library/Frameworks/R.framework",
            "/Library/Frameworks/R.framework/Resources/jags",
    };
    for (const QString &folder: paths) {
        if (pathExists(folder)) {
            ui->install_btn->setEnabled(false);
            return true;  // If any folder exists, return true
        }
    }
    ui->install_btn->setEnabled(true);
    ui->overdisp_btn->setEnabled(false);
    ui->run_btn->setEnabled(false);
    return false;
}

bool MainWindow::pathExists(const QString &path) {
    // Check if the source path exists
    QFileInfo sourceInfo(path);
    if (!sourceInfo.exists()) {
        qDebug() << "Source path does not exist:" << path;
        return false;
    }
    return true;
}

void MainWindow::organizeFiles() {
    // Get all ui elements that are CustomListView
    QList<CustomListView *> listViews = findChildren<CustomListView *>();
    for (CustomListView *listView: listViews) {
        if (listView->objectName() == "file_list") {
            continue;
        }
        if (listView->count() > 0) {
            if (listView->objectName().startsWith("vec")) {
                if (listView->objectName().contains("nonsel")) {
                    if (listView->objectName().contains("_2")) {
                        stat->addGeneCount(
                                new GeneCountFile(listView->item(0)->data(Qt::UserRole).toUrl().toLocalFile(),
                                                  VECTOR,
                                                  NONSELECTED, TWO));
                    } else {
                        stat->addGeneCount(
                                new GeneCountFile(listView->item(0)->data(Qt::UserRole).toUrl().toLocalFile(),
                                                  VECTOR,
                                                  NONSELECTED, ONE));
                    }
                } else {
                    if (listView->objectName().contains("_2")) {
                        stat->addGeneCount(
                                new GeneCountFile(listView->item(0)->data(Qt::UserRole).toUrl().toLocalFile(),
                                                  VECTOR,
                                                  SELECTED, TWO));

                    } else {
                        stat->addGeneCount(
                                new GeneCountFile(listView->item(0)->data(Qt::UserRole).toUrl().toLocalFile(),
                                                  VECTOR,
                                                  SELECTED, ONE));
                    }
                }
            } else {
                if (listView->objectName().contains("bait1")) {
                    if (listView->objectName().contains("nonsel")) {
                        if (listView->objectName().contains("_2")) {
                            stat->addGeneCount(
                                    new GeneCountFile(listView->item(0)->data(Qt::UserRole).toUrl().toLocalFile(),
                                                      BAIT,
                                                      NONSELECTED, TWO));
                        } else {
                            stat->addGeneCount(
                                    new GeneCountFile(listView->item(0)->data(Qt::UserRole).toUrl().toLocalFile(),
                                                      BAIT,
                                                      NONSELECTED, ONE));
                        }
                    } else {
                        if (listView->objectName().contains("_2")) {
                            stat->addGeneCount(
                                    new GeneCountFile(listView->item(0)->data(Qt::UserRole).toUrl().toLocalFile(),
                                                      BAIT,
                                                      SELECTED, TWO));
                        } else {
                            stat->addGeneCount(
                                    new GeneCountFile(listView->item(0)->data(Qt::UserRole).toUrl().toLocalFile(),
                                                      BAIT,
                                                      SELECTED, ONE));
                        }
                    }
                } else if (listView->objectName().contains("bait2")) {
                    if (listView->objectName().contains("nonsel")) {
                        if (listView->objectName().contains("_2")) {
                            stat->addGeneCount(
                                    new GeneCountFile(listView->item(0)->data(Qt::UserRole).toUrl().toLocalFile(),
                                                      BAIT,
                                                      NONSELECTED, TWO));
                        } else {
                            stat->addGeneCount(
                                    new GeneCountFile(listView->item(0)->data(Qt::UserRole).toUrl().toLocalFile(),
                                                      BAIT,
                                                      NONSELECTED, ONE));
                        }
                    } else {
                        if (listView->objectName().contains("_2")) {
                            stat->addGeneCount(
                                    new GeneCountFile(listView->item(0)->data(Qt::UserRole).toUrl().toLocalFile(),
                                                      BAIT,
                                                      SELECTED, TWO));
                        } else {
                            stat->addGeneCount(
                                    new GeneCountFile(listView->item(0)->data(Qt::UserRole).toUrl().toLocalFile(),
                                                      BAIT,
                                                      SELECTED, ONE));
                        }
                    }
                }
            }
        }
//        for (int i = 0; i < listView->count(); ++i) {
//            QListWidgetItem *item = listView->item(i);
//            organizedFiles.insert(item->data(Qt::UserRole).toUrl().toLocalFile(), item->text());
//        }
    }
}

void MainWindow::writeRInputFile(const QString &outputFile) {
    QFile file(outputFile);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qDebug() << "Could not open file for writing:" << outputFile;
        return;
    }
    QTextStream out(&file);
    out << "Vector_Selected_1=" << stat->getGeneCount(VECTOR, SELECTED, ONE)->filepath << "\n";
    out << "Vector_Selected_2=" << stat->getGeneCount(VECTOR, SELECTED, TWO)->filepath << "\n";
//    out << "Vector_Non-Selected_1=" << files["vec1_nonsel"] << "\n";
//    out << "Vector_Non-Selected_2=" << files["vec2_nonsel"] << "\n";
//    out << "Bait1_Selected_1=" << files["bait1_sel_1"] << "\n";
//    out << "Bait1_Selected_2=" << files["bait1_sel_2"] << "\n";
//    out << "Bait1_Non-Selected_1=" << files["bait1_nonsel_1"] << "\n";
//    out << "Bait1_Non-Selected_2=" << files["bait1_nonsel_2"] << "\n";
//    out << "Bait2_Selected_1=" << files["bait2_sel_1"] << "\n";
//    out << "Bait2_Selected_2=" << files["bait2_sel_2"] << "\n";
    file.close();
}

void MainWindow::on_run_btn_clicked() {
    writeRInputFile("/Volumes/Vault/Download.Vault/gc_files/r_input.params");
}

bool MainWindow::copyFolderToFrameworks(const QString &sourcePath, const QString &destinationPath) {

    if (pathExists(QString("/Library/Frameworks/R.Framework"))) {
        qDebug() << "R.Framework exists";
        QMessageBox::StandardButton reply;
        reply = QMessageBox::question(nullptr, "Overwrite Confirmation",
                                      "R.Framework already exists and does not support Stat Maker++ operation. Do you want to overwrite it?",
                                      QMessageBox::Yes | QMessageBox::No);
        if (reply == QMessageBox::Yes) {
            qDebug() << "User chose to overwrite";
            // Add code here to handle the overwriting
        } else {
            qDebug() << "User chose not to overwrite";
            // Add code here to handle the cancellation
        }
    } else {
        qDebug() << "R.Framework does not exist";
    }

    // Ensure the destination path exists
    if (access(destinationPath.toUtf8().constData(), F_OK) == -1) {
        qDebug() << "Destination path does not exist:" << destinationPath;
        return false;
    }

    // Create an authorization reference
    AuthorizationRef authRef;
    OSStatus status = AuthorizationCreate(nullptr, kAuthorizationEmptyEnvironment, kAuthorizationFlagDefaults,
                                          &authRef);
    if (status != errAuthorizationSuccess) {
        qDebug() << "Authorization creation failed";
        return false;
    }

    // Set up the authorization rights
    AuthorizationItem right = {kAuthorizationRightExecute, 0, nullptr, 0};
    AuthorizationRights rights = {1, &right};

    AuthorizationFlags flags = kAuthorizationFlagDefaults | kAuthorizationFlagInteractionAllowed |
                               kAuthorizationFlagPreAuthorize | kAuthorizationFlagExtendRights;

    status = AuthorizationCopyRights(authRef, &rights, nullptr, flags, nullptr);
    if (status != errAuthorizationSuccess) {
        qDebug() << "Authorization rights request failed";
        AuthorizationFree(authRef, kAuthorizationFlagDefaults);
        return false;
    }

    // Convert paths to QByteArray and then to char*
    QByteArray sourcePathBA = sourcePath.toUtf8();
    QByteArray destinationPathBA = destinationPath.toUtf8();

    char *args[] = {const_cast<char *>("-f"),
                    const_cast<char *>("-R"),
                    sourcePathBA.data(),
                    destinationPathBA.data(),
                    nullptr};


    FILE *pipe = nullptr;
    status = AuthorizationExecuteWithPrivileges(authRef, "/bin/cp", kAuthorizationFlagDefaults, args, &pipe);
    if (status != errAuthorizationSuccess) {
        qDebug() << "AuthorizationExecuteWithPrivileges failed with status:" << status;
        AuthorizationFree(authRef, kAuthorizationFlagDefaults);
        return false;
    }

    // Close the pipe and free the authorization reference
    if (pipe) {
        fclose(pipe);
    }
    AuthorizationFree(authRef, kAuthorizationFlagDefaults);

    qDebug() << "Folder copied successfully";
    return true;
}

void MainWindow::on_install_btn_clicked() {
    qDebug() << "Installing R Framework";
}



#include "mainwindow.h"

#include <QApplication>

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  a.setApplicationName("JunctionDice++");
  MainWindow w(argc, argv);
  w.show();
  return a.exec();
}

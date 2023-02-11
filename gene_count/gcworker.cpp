#include "gcworker.h"

#include <SamFile.h>

#include <QApplication>

GCWorker::GCWorker(QString file) { sam_file = file; }

void GCWorker::doWork() {
  timer.start();
  SamFileHeader header;
  // Open the bam file for reading and read the header.
  SamFile samIn(sam_file.toStdString().c_str(), SamFile::READ, &header);
  SamRecord record;
  GCStat stat;
  while (samIn.ReadRecord(header, record)) {
    //    Print the reference positions associated with this read.
    stat.filename = sam_file;
    stat.readCount = samIn.GetCurrentPosition();
    stat.elapsedTime = timer.elapsed() / 1000;
    stat.currentRefSeq = record.getReadName();
    if (stat.readCount % 1000 == 0) {
      emit sig->gc_update_progress_sig(stat);
      QApplication::processEvents();
    }
    Cigar* cigar = record.getCigarInfo();
    for (int i = 0; i < record.getReadLength(); i++) {
      int refPos = cigar->getRefPosition(i, record.get1BasedPosition());
      if (refPos != Cigar::INDEX_NA) {
      }
    }
  }
  stat.running = false;
  emit sig->gc_update_progress_sig(stat);
}

void GCWorker::finished() {}

void GCWorker::stop() {}

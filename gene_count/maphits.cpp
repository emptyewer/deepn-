#include "maphits.h"

#include <QDebug>

MapHit::MapHit(const QString& line)
    : identity_(0.0),
      bitscore_(0.0),
      read_(""),
      gene_(""),
      nm_number_(""),
      orfstart_(-1),
      orfend_(-1),
      length_(0),
      qstart_(-1),
      qend_(-1),
      rstart_(-1),
      rend_(-1) {
  process(line);
}

MapHit::MapHit(const MapHit& other) {
  identity_ = other.identity_;
  bitscore_ = other.bitscore_;
  read_ = other.read_;
  gene_ = other.gene_;
  nm_number_ = other.nm_number_;
  orfstart_ = other.orfstart_;
  orfend_ = other.orfend_;
  length_ = other.length_;
  qstart_ = other.qstart_;
  qend_ = other.qend_;
  rstart_ = other.rstart_;
  rend_ = other.rend_;
}

void MapHit::process(const QString& line) {
  QStringList split_line = line.trimmed().split("\t");
  read_ = split_line.at(0);
  QString map = split_line.at(1);
  identity_ = split_line.at(2).toDouble();
  length_ = split_line.at(3).toInt();
  qstart_ = split_line.at(6).toInt();
  qend_ = split_line.at(7).toInt();
  rstart_ = split_line.at(8).toInt();
  rend_ = split_line.at(9).toInt();
  bitscore_ = split_line.at(11).toDouble();

  QMap<QString, QString> split_map = breakDownMap(map);
  gene_ = split_map["GeneSymbol"];
  nm_number_ = split_map["RefSeqNM"];
  orfstart_ = split_map["ORFstart"].toInt();
  orfend_ = split_map["ORFstop"].toInt();
}

QMap<QString, QString> MapHit::breakDownMap(const QString& inmap) {
  QMap<QString, QString> outmap;
  QStringList com_list = inmap.split("|");
  for (const QString& com : qAsConst(com_list)) {
    if (com.length() > 0) {
      QStringList kv_list = com.split(":");
      outmap.insert(kv_list.at(0), kv_list.at(1));
    }
  }
  return outmap;
}

QString MapHit::frame() const {
  if (rend_ < rstart_) {
    return "backward";
  }
  int rem = qAbs(rstart_ - orfstart_ - qstart_ + 1) % 3;
  if (rem == 0) {
    return "+0_frame";
  } else if (rem == 1) {
    return "+1_frame";
  } else if (rem == 2) {
    return "+2_frame";
  }
  return "";  // Return empty string if none of the above conditions are true
}

QString MapHit::location() const {
  if (rstart_ < orfstart_) {
    return "upstream";
  } else if (rend_ > orfend_) {
    return "downstream";
  }
  return "in_orf";  // If none of the above conditions are true, return "in_orf"
}

QString MapHit::toString() const {
  QString output;
  QTextStream stream(&output);

  stream << "Read: " << read_ << "\n"
         << "Gene: " << gene_ << "\n"
         << "NM Number: " << nm_number_ << "\n"
         << "Identity: " << QString::number(identity_, 'f', 2) << "%\n"
         << "Bit Score: " << QString::number(bitscore_, 'f', 2) << "\n"
         << "ORF Start: " << orfstart_ << "\n"
         << "ORF End: " << orfend_ << "\n"
         << "Length: " << length_ << "\n"
         << "Query Start: " << qstart_ << "\n"
         << "Query End: " << qend_ << "\n"
         << "Ref Start: " << rstart_ << "\n"
         << "Ref End: " << rend_;

  return output;
}

// ReadHits

ReadHits::ReadHits() {}

ReadHits::ReadHits(const ReadHits& other) {
  m_hits = other.m_hits;
  m_read = other.m_read;
}

void ReadHits::add(MapHit hit) {
  MapHit my_hit(hit);
  m_hits.append(my_hit);
  if (!allReadsEqual()) {
    qDebug() << "WARNING : adding read" << my_hit.read_ << "to"
             << m_hits.last().read_;
  } else {
    m_read = m_hits.first().read_;
  }
}

void ReadHits::clear() { m_hits.clear(); }

int ReadHits::size() const { return m_hits.size(); }

bool ReadHits::allReadsEqual() const {
  QStringList readNames;
  foreach (const MapHit& hit, m_hits) {
    readNames.append(hit.read_);
  }
  return readNames.count(readNames.first()) == readNames.size();
}

double ReadHits::lastBitscore() const { return m_hits.last().bitscore_; }

bool ReadHits::operator==(MapHit other) const { return other.read_ == m_read; }

QStringList ReadHits::genes() const {
  QStringList genes;
  foreach (const MapHit& hit, m_hits) {
    if (!genes.contains(hit.gene_)) {
      genes.append(hit.gene_);
    }
  }
  return genes;
}

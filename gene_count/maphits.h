#ifndef MAPHITS_H
#define MAPHITS_H

#include <QMap>
#include <QString>

class MapHit {
 public:
  MapHit(const QString& line);
  MapHit(const MapHit& other);

  void process(const QString& line);
  QString frame() const;
  QString location() const;
  double identity_;
  double bitscore_;
  QString read_;
  QString gene_;
  QString nm_number_;
  int orfstart_;
  int orfend_;
  int length_;
  int qstart_;
  int qend_;
  int rstart_;
  int rend_;

  QString toString() const;

 private:
  QMap<QString, QString> breakDownMap(const QString& inmap);
};

class ReadHitsIterator {
 public:
  ReadHitsIterator(QList<MapHit>::iterator it) : it_(it) {}
  ReadHitsIterator& operator++() {
    ++it_;
    return *this;
  }
  bool operator!=(const ReadHitsIterator& other) const {
    return it_ != other.it_;
  }
  MapHit& operator*() const { return *it_; }
 private:
  QList<MapHit>::iterator it_;
};

class ReadHits {
 public:
  ReadHits();
  ReadHits(const ReadHits& other);
  void add(MapHit hit);
  void clear();
  int size() const;
  bool allReadsEqual() const;
  double lastBitscore() const;
  bool operator==(MapHit other) const;
  QStringList genes() const;
  using iterator = ReadHitsIterator;
  iterator begin() { return iterator(m_hits.begin()); }
  iterator end() { return iterator(m_hits.end()); }

 private:
  QList<MapHit> m_hits;
  QString m_read;
};

#endif // MAPHITS_H

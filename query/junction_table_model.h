#ifndef JUNCTION_TABLE_MODEL_H
#define JUNCTION_TABLE_MODEL_H

#include <data_structures.h>

#include <QAbstractTableModel>
#include <QMap>
#include <QVector>

class JunctionTableModel : public QAbstractTableModel
{
    Q_OBJECT

public:
    enum Column {
        Position = 0,
        PPM,
        QueryStart,
        CDS,
        Frame,
        RawCount,
        RefSeq,
        Enrichment,
        ColumnCount
    };

    explicit JunctionTableModel(QObject* parent = nullptr);

    void setJunctions(const QVector<deepn::JunctionSite>& sites);
    void setCollapsedJunctions(const QVector<deepn::CollapsedJunction>& collapsed);
    void setComparisonData(const QVector<deepn::JunctionSite>& secondarySites);
    void clearComparisonData();
    deepn::JunctionSite junctionAt(int row) const;

    int rowCount(const QModelIndex& parent = {}) const override;
    int columnCount(const QModelIndex& parent = {}) const override;
    QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const override;
    void sort(int column, Qt::SortOrder order) override;

private:
    QVector<deepn::JunctionSite> m_sites;
    QMap<int, double> m_secondaryPpmByPos;  // position -> PPM from secondary dataset
    bool m_isCollapsed = false;
    bool m_hasComparison = false;
};

#endif // JUNCTION_TABLE_MODEL_H

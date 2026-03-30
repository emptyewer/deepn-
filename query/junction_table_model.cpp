#include "junction_table_model.h"

#include <QColor>

#include <algorithm>

using namespace deepn;

JunctionTableModel::JunctionTableModel(QObject* parent)
    : QAbstractTableModel(parent)
{
}

void JunctionTableModel::setJunctions(const QVector<JunctionSite>& sites)
{
    beginResetModel();
    m_sites = sites;
    m_isCollapsed = false;
    endResetModel();
}

void JunctionTableModel::setCollapsedJunctions(const QVector<CollapsedJunction>& collapsed)
{
    beginResetModel();
    m_sites.clear();
    m_sites.reserve(collapsed.size());

    // Convert collapsed entries to JunctionSite for unified display
    for (const auto& cj : collapsed) {
        JunctionSite site;
        site.position = cj.position;
        site.ppm = cj.totalPpm;
        site.queryStart = cj.variantCount;  // repurpose: show variant count
        site.cdsClass = cj.cdsClass;
        site.frame = cj.dominantFrame;
        site.rawCount = cj.variantCount;

        if (!cj.variants.isEmpty()) {
            site.refseq = cj.variants.first().refseq;
            site.geneName = cj.variants.first().geneName;
        }

        m_sites.append(site);
    }

    m_isCollapsed = true;
    endResetModel();
}

JunctionSite JunctionTableModel::junctionAt(int row) const
{
    if (row < 0 || row >= m_sites.size()) return {};
    return m_sites.at(row);
}

int JunctionTableModel::rowCount(const QModelIndex& parent) const
{
    if (parent.isValid()) return 0;
    return m_sites.size();
}

int JunctionTableModel::columnCount(const QModelIndex& parent) const
{
    if (parent.isValid()) return 0;
    return ColumnCount;
}

QVariant JunctionTableModel::data(const QModelIndex& index, int role) const
{
    if (!index.isValid()) return {};
    if (index.row() < 0 || index.row() >= m_sites.size()) return {};

    const auto& site = m_sites.at(index.row());

    if (role == Qt::DisplayRole) {
        switch (index.column()) {
        case Position:   return site.position;
        case PPM:        return QString::number(site.ppm, 'f', 2);
        case QueryStart:
            if (m_isCollapsed) return QString("%1 variants").arg(site.rawCount);
            return site.queryStart;
        case CDS:        return site.cdsLabel();
        case Frame:      return site.frameLabel();
        case RawCount:   return site.rawCount;
        case RefSeq:     return site.refseq;
        default:         return {};
        }
    }

    if (role == Qt::TextAlignmentRole) {
        switch (index.column()) {
        case Position:
        case PPM:
        case QueryStart:
        case RawCount:
            return QVariant(Qt::AlignRight | Qt::AlignVCenter);
        default:
            return QVariant(Qt::AlignLeft | Qt::AlignVCenter);
        }
    }

    if (role == Qt::ForegroundRole) {
        // Color-code by junction category
        JunctionCategory cat = classifyJunction(site);
        QColor color;
        switch (cat) {
        case JunctionCategory::InOrfInFrame:     color = QColor("#1E40AF"); break;
        case JunctionCategory::UpstreamInFrame:  color = QColor("#0D9488"); break;
        case JunctionCategory::InOrfOutOfFrame:  color = QColor("#D97706"); break;
        case JunctionCategory::DownstreamOrBack: color = QColor("#6B7280"); break;
        }
        // Use a slightly darker version for readability on white background
        return color.darker(130);
    }

    if (role == Qt::ToolTipRole) {
        QString tip = QString("Position: %1\nPPM: %2\nFrame: %3\nCDS: %4")
                          .arg(site.position)
                          .arg(site.ppm, 0, 'f', 2)
                          .arg(site.frameLabel())
                          .arg(site.cdsLabel());
        if (m_isCollapsed) {
            tip += QString("\nVariants: %1").arg(site.rawCount);
        }
        return tip;
    }

    // For sorting, provide raw numeric data
    if (role == Qt::UserRole) {
        switch (index.column()) {
        case Position:   return site.position;
        case PPM:        return site.ppm;
        case QueryStart: return site.queryStart;
        case RawCount:   return site.rawCount;
        default:         return {};
        }
    }

    return {};
}

QVariant JunctionTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (orientation != Qt::Horizontal || role != Qt::DisplayRole) return {};

    switch (section) {
    case Position:   return "Position";
    case PPM:        return "PPM";
    case QueryStart: return m_isCollapsed ? "Variants" : "QueryStart";
    case CDS:        return "CDS";
    case Frame:      return "Frame";
    case RawCount:   return m_isCollapsed ? "Variants" : "Raw Count";
    case RefSeq:     return "RefSeq";
    default:         return {};
    }
}

void JunctionTableModel::sort(int column, Qt::SortOrder order)
{
    if (m_sites.isEmpty()) return;

    beginResetModel();

    auto comparator = [column, order](const JunctionSite& a, const JunctionSite& b) -> bool {
        bool lessThan = false;

        switch (column) {
        case Position:
            lessThan = a.position < b.position;
            break;
        case PPM:
            lessThan = a.ppm < b.ppm;
            break;
        case QueryStart:
            lessThan = a.queryStart < b.queryStart;
            break;
        case CDS:
            lessThan = a.cdsClass < b.cdsClass;
            break;
        case Frame:
            lessThan = a.frame < b.frame;
            break;
        case RawCount:
            lessThan = a.rawCount < b.rawCount;
            break;
        case RefSeq:
            lessThan = a.refseq < b.refseq;
            break;
        default:
            lessThan = a.position < b.position;
            break;
        }

        return (order == Qt::AscendingOrder) ? lessThan : !lessThan;
    };

    std::stable_sort(m_sites.begin(), m_sites.end(), comparator);

    endResetModel();
}

//
// Created by Venky Krishnamani on 2024/7/12.
//

#ifndef STATMAKER_DATASTRUCTS_H
#define STATMAKER_DATASTRUCTS_H

#include <QString>
#include <QMap>
#include <utility>

enum FileType {
    BAIT,
    VECTOR
};

enum SelectionType {
    NONSELECTED,
    SELECTED
};

enum GroupType {
    ONE,
    TWO
};

struct GeneCountFile {
    QString name;
    QString filepath;
    QMap<QString, int> counts;
    int total;
    FileType type;
    SelectionType selection;
    GroupType group;

    explicit GeneCountFile(QString filepath, FileType type, SelectionType selection, GroupType group) {
        this->name = std::move(name);
        this->filepath = std::move(filepath);
        this->counts = std::move(counts);
        this->total = total;
        this->type = type;
        this->selection = selection;
        this->group = group;
    }
};

struct Stat {
    QList<GeneCountFile *> *geneCounts;
    int threshold;

    void addGeneCount(GeneCountFile *geneCount) const {
        geneCounts->append(geneCount);
    }

    GeneCountFile *getGeneCount(const FileType ftype, const SelectionType stype, const GroupType gtype) {
        for (auto geneCount: *geneCounts) {
            if (geneCount->type == ftype && geneCount->selection == stype && geneCount->group == gtype) {
                return geneCount;
            }
        }
        return nullptr;
    }
};

#endif //STATMAKER_DATASTRUCTS_H

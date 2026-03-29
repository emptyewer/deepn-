#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace deseq2 {

struct EnrichmentResult {
    std::string gene;
    std::string bait;
    double wald_stat = 0;
    double pvalue = 1;
    double log2FC = 0;
    double rank_score = 0;
    double kde_contribution = 0;
    double total_score = 0;
};

struct SpecificityResult {
    std::string gene;
    std::string bait;
    double mean_spec_score = 0;
    double kde_contribution = 0;
    double total_score = 0;
};

struct JunctionData {
    std::string transcript_id;
    std::string gene;
    int in_frame_reads_selected = 0;
    int total_reads_selected = 0;
    int in_frame_reads_nonselected = 0;
    int total_reads_nonselected = 0;
};

struct InFrameResult {
    std::string gene;
    std::string bait;
    std::string transcripts;
    double in_frame_prop_s = 0;
    double in_frame_prop_ns = 0;
    double z_statistic = 0;
    double freq_score = 0;
};

struct Y2HScore {
    std::string gene;
    std::string bait;
    double enrichment_score = 0;
    double specificity_score = 0;
    double in_frame_score = 0;
    std::string in_frame_transcripts;
    double sum_scores = 0;
    double borda_score = 0;
};

// Enrichment Score: DESeq2 S vs N → rank by Wald stat → KDE weighting
class EnrichmentScorer {
public:
    std::vector<EnrichmentResult> compute(
        const Eigen::MatrixXd& deseq_results, // 6-col: baseMean, lfc, se, stat, pval, padj
        const std::vector<std::string>& gene_names,
        const std::string& bait_name,
        double p_threshold = 1.0,
        double fc_threshold = 0.0,
        int n_kde_bins = 100
    );
};

// Specificity Score: pairwise bait contrasts → rank → KDE weighting
// Input: map of "baitA_vs_baitB" → DESeq2 6-column results matrix
struct PairwiseContrast {
    std::string bait_numerator;
    std::string bait_denominator;
    Eigen::MatrixXd results; // 6-col DESeq2 results
    std::vector<std::string> gene_names;
};

class SpecificityScorer {
public:
    std::vector<SpecificityResult> compute(
        const std::vector<PairwiseContrast>& contrasts,
        int n_baits_in_group,
        double p_threshold = 1.0,
        double fc_threshold = 0.0,
        int n_kde_bins = 100
    );
};

// In-Frame Score: junction read frame analysis (two-proportion z-test)
class InFrameScorer {
public:
    std::vector<InFrameResult> compute(
        const std::vector<JunctionData>& junction_data,
        const std::string& bait_name
    );
};

// Borda rank aggregation
class BordaAggregator {
public:
    std::vector<Y2HScore> aggregate(
        const std::vector<EnrichmentResult>& enrichment,
        const std::vector<SpecificityResult>& specificity,
        const std::vector<InFrameResult>& in_frame
    );
};

} // namespace deseq2

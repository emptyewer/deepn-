#include "y2h_scores.h"
#include <iostream>
#include <set>
#include <sstream>

namespace deseq2 {

// Helper: apply KDE-based relative fold-change weighting
// Replicates the R `get_rel_fc_score` function from Y2H-SCORES
static void applyKdeWeighting(
    std::vector<double>& rank_scores,
    const std::vector<double>& fc_values,
    std::vector<double>& total_scores,
    int n_bins)
{
    if (rank_scores.empty()) return;

    // Find bin edges (evenly spaced across rank_score range)
    double min_score = *std::min_element(rank_scores.begin(), rank_scores.end());
    double max_score = *std::max_element(rank_scores.begin(), rank_scores.end());
    double bin_width = (max_score - min_score + 1e-10) / n_bins;

    total_scores.resize(rank_scores.size());

    // Assign each element to a bin
    std::vector<int> bin_indices(rank_scores.size());
    for (size_t i = 0; i < rank_scores.size(); i++) {
        bin_indices[i] = std::min(n_bins - 1,
            static_cast<int>((rank_scores[i] - min_score) / bin_width));
    }

    // For each bin, rank by fold change and compute relative contribution
    for (int b = 0; b < n_bins; b++) {
        // Collect indices in this bin
        std::vector<size_t> in_bin;
        for (size_t i = 0; i < bin_indices.size(); i++) {
            if (bin_indices[i] == b) in_bin.push_back(i);
        }
        if (in_bin.empty()) continue;

        // Rank by FC within bin
        std::vector<std::pair<double, size_t>> fc_rank;
        for (size_t idx : in_bin) {
            fc_rank.push_back({fc_values[idx], idx});
        }
        std::sort(fc_rank.begin(), fc_rank.end());

        double max_rank = static_cast<double>(fc_rank.size());
        double bin_min_score = rank_scores[in_bin[0]];
        double bin_max_score = rank_scores[in_bin[0]];
        for (size_t idx : in_bin) {
            bin_min_score = std::min(bin_min_score, rank_scores[idx]);
            bin_max_score = std::max(bin_max_score, rank_scores[idx]);
        }

        double bin_count = static_cast<double>(in_bin.size());

        for (size_t r = 0; r < fc_rank.size(); r++) {
            size_t idx = fc_rank[r].second;
            double rel_fc_score = (r + 1.0) / max_rank;
            double contribution = rel_fc_score *
                (bin_max_score - bin_min_score + 0.0001) / bin_count;
            total_scores[idx] = rank_scores[idx] + contribution;
        }
    }
}

// Helper: compute ranks (1 = smallest, ties get average rank)
static std::vector<double> computeRanks(const std::vector<double>& values) {
    int n = values.size();
    std::vector<double> ranks(n);
    std::vector<std::pair<double, int>> sorted;
    for (int i = 0; i < n; i++) {
        sorted.push_back({values[i], i});
    }
    std::sort(sorted.begin(), sorted.end());

    int i = 0;
    while (i < n) {
        int j = i;
        while (j < n && sorted[j].first == sorted[i].first) j++;
        double avg_rank = (i + j + 1.0) / 2.0; // average rank for ties
        for (int k = i; k < j; k++) {
            ranks[sorted[k].second] = avg_rank;
        }
        i = j;
    }
    return ranks;
}

// ============================================================================
// EnrichmentScorer
// ============================================================================
std::vector<EnrichmentResult> EnrichmentScorer::compute(
    const Eigen::MatrixXd& deseq_results,
    const std::vector<std::string>& gene_names,
    const std::string& bait_name,
    double p_threshold,
    double fc_threshold,
    int n_kde_bins)
{
    int n_genes = deseq_results.rows();
    std::vector<EnrichmentResult> results;

    // Filter: pvalue < threshold AND log2FC > fc_threshold
    for (int i = 0; i < n_genes; i++) {
        double pval = deseq_results(i, 4);  // raw p-value
        double lfc = deseq_results(i, 1);   // log2FC
        double stat = deseq_results(i, 3);  // Wald statistic

        if (std::isnan(pval) || std::isnan(lfc) || std::isnan(stat)) continue;
        if (pval >= p_threshold || lfc <= fc_threshold) continue;

        EnrichmentResult r;
        r.gene = gene_names[i];
        r.bait = bait_name;
        r.wald_stat = stat;
        r.pvalue = pval;
        r.log2FC = lfc;
        results.push_back(r);
    }

    if (results.empty()) return results;

    // Rank by -stat (higher stat = lower rank number = higher score)
    std::vector<double> neg_stats;
    for (auto& r : results) neg_stats.push_back(-r.wald_stat);
    std::vector<double> ranks = computeRanks(neg_stats);

    double max_rank = *std::max_element(ranks.begin(), ranks.end());
    std::vector<double> rank_scores(results.size());
    std::vector<double> fc_values(results.size());

    for (size_t i = 0; i < results.size(); i++) {
        rank_scores[i] = (max_rank - ranks[i]) / max_rank;
        results[i].rank_score = rank_scores[i];
        fc_values[i] = results[i].log2FC;
    }

    // Apply KDE weighting
    std::vector<double> total_scores;
    applyKdeWeighting(rank_scores, fc_values, total_scores, n_kde_bins);

    for (size_t i = 0; i < results.size(); i++) {
        results[i].kde_contribution = total_scores[i] - results[i].rank_score;
        results[i].total_score = total_scores[i];
    }

    return results;
}

// ============================================================================
// SpecificityScorer
// ============================================================================
std::vector<SpecificityResult> SpecificityScorer::compute(
    const std::vector<PairwiseContrast>& contrasts,
    int n_baits_in_group,
    double p_threshold,
    double fc_threshold,
    int n_kde_bins)
{
    // Collect per-bait spec scores from pairwise contrasts
    // Key: "gene\tbait" → vector of bait_spec_scores
    struct SpecEntry {
        double stat;
        double log2FC;
        double pvalue;
    };
    std::map<std::string, std::vector<SpecEntry>> bait_entries;

    for (const auto& contrast : contrasts) {
        int n = contrast.results.rows();
        for (int i = 0; i < n; i++) {
            double lfc = contrast.results(i, 1);
            double pval = contrast.results(i, 4);
            double stat = contrast.results(i, 3);
            if (std::isnan(pval) || std::isnan(lfc) || std::isnan(stat)) continue;

            // Positive LFC → numerator bait is enriched
            // Negative LFC → denominator bait is enriched
            std::string bait;
            if (lfc > 0) {
                bait = contrast.bait_numerator;
            } else if (lfc < 0) {
                bait = contrast.bait_denominator;
            } else {
                continue;
            }

            std::string key = contrast.gene_names[i] + "\t" + bait;
            bait_entries[key].push_back({stat, lfc, pval});
        }
    }

    // Filter by thresholds, collect all entries for ranking
    struct RankEntry {
        std::string gene;
        std::string bait;
        double abs_stat;
        double abs_lfc;
        double pvalue;
    };
    std::vector<RankEntry> all_entries;

    for (auto& [key, entries] : bait_entries) {
        for (auto& e : entries) {
            if (e.pvalue < p_threshold && std::abs(e.log2FC) > fc_threshold) {
                auto tab_pos = key.find('\t');
                all_entries.push_back({
                    key.substr(0, tab_pos),
                    key.substr(tab_pos + 1),
                    std::abs(e.stat),
                    std::abs(e.log2FC),
                    e.pvalue
                });
            }
        }
    }

    if (all_entries.empty()) return {};

    // Rank by -|stat|
    std::vector<double> neg_abs_stats;
    for (auto& e : all_entries) neg_abs_stats.push_back(-e.abs_stat);
    std::vector<double> ranks = computeRanks(neg_abs_stats);

    double max_rank = *std::max_element(ranks.begin(), ranks.end());

    // Compute bait_spec_score per entry
    std::vector<double> bait_spec_scores(all_entries.size());
    for (size_t i = 0; i < all_entries.size(); i++) {
        bait_spec_scores[i] = (max_rank - ranks[i]) / max_rank;
    }

    // Group by gene+bait, average spec_score, collect LFCs
    struct GroupData {
        std::vector<double> scores;
        std::vector<double> lfcs;
    };
    std::map<std::string, GroupData> grouped;

    for (size_t i = 0; i < all_entries.size(); i++) {
        std::string key = all_entries[i].gene + "\t" + all_entries[i].bait;
        grouped[key].scores.push_back(bait_spec_scores[i]);
        grouped[key].lfcs.push_back(all_entries[i].abs_lfc);
    }

    // Compute mean spec_score per gene-bait
    std::vector<double> mean_scores;
    std::vector<double> mean_lfcs;
    std::vector<std::string> keys;

    for (auto& [key, data] : grouped) {
        double mean_s = 0;
        for (double s : data.scores) mean_s += s;
        mean_s /= n_baits_in_group;

        double mean_l = 0;
        for (double l : data.lfcs) mean_l += l;
        mean_l /= data.lfcs.size();

        keys.push_back(key);
        mean_scores.push_back(mean_s);
        mean_lfcs.push_back(mean_l);
    }

    // Apply KDE weighting
    std::vector<double> total_scores;
    applyKdeWeighting(mean_scores, mean_lfcs, total_scores, n_kde_bins);

    // Build results
    std::vector<SpecificityResult> results;
    for (size_t i = 0; i < keys.size(); i++) {
        auto tab_pos = keys[i].find('\t');
        SpecificityResult r;
        r.gene = keys[i].substr(0, tab_pos);
        r.bait = keys[i].substr(tab_pos + 1);
        r.mean_spec_score = mean_scores[i];
        r.kde_contribution = total_scores[i] - mean_scores[i];
        r.total_score = total_scores[i];
        results.push_back(r);
    }

    return results;
}

// ============================================================================
// InFrameScorer
// ============================================================================
std::vector<InFrameResult> InFrameScorer::compute(
    const std::vector<JunctionData>& junction_data,
    const std::string& bait_name)
{
    // Group by transcript, summing reads across replicates
    struct TranscriptSums {
        std::string gene;
        int in_frame_s = 0;
        int total_s = 0;
        int in_frame_ns = 0;
        int total_ns = 0;
    };
    std::map<std::string, TranscriptSums> by_transcript;

    for (const auto& jd : junction_data) {
        auto& ts = by_transcript[jd.transcript_id];
        ts.gene = jd.gene;
        ts.in_frame_s += jd.in_frame_reads_selected;
        ts.total_s += jd.total_reads_selected;
        ts.in_frame_ns += jd.in_frame_reads_nonselected;
        ts.total_ns += jd.total_reads_nonselected;
    }

    // Filter: require selected reads > 0 and non-selected reads > 0
    // (or handle missing NS data)
    struct ScoredTranscript {
        std::string transcript_id;
        std::string gene;
        double in_frame_prop_s;
        double in_frame_prop_ns;
        double z_stat;
        int in_frame_s;
    };
    std::vector<ScoredTranscript> scored;

    for (auto& [tid, ts] : by_transcript) {
        if (ts.total_s <= 0) continue;

        bool has_ns = (ts.total_ns > 0);
        // Filter: require both S and NS reads > 0 when NS data exists
        if (has_ns && (ts.total_s <= 0 || ts.total_ns <= 0)) continue;
        if (has_ns && ts.total_s < ts.total_ns) continue;

        double prop_s = static_cast<double>(ts.in_frame_s) / ts.total_s;
        double prop_ns;

        if (has_ns && ts.total_ns > 0) {
            prop_ns = static_cast<double>(ts.in_frame_ns) / ts.total_ns;
        } else {
            prop_ns = 1.0 / 3.0; // assume random 3-frame baseline
        }

        // Two-proportion z-test
        double n_s = ts.total_s;
        double n_ns = has_ns ? ts.total_ns : n_s; // if no NS, use S count
        double pi_hat = (prop_s * n_s + prop_ns * n_ns) / (n_s + n_ns);

        double denom = std::sqrt(pi_hat * (1.0 - pi_hat) * (1.0 / n_s + 1.0 / n_ns));
        if (denom <= 0 || std::isnan(denom) || pi_hat <= 0 || pi_hat >= 1) continue;

        double z = (prop_s - prop_ns) / denom;
        if (std::isnan(z)) continue;

        scored.push_back({tid, ts.gene, prop_s, prop_ns, z, ts.in_frame_s});
    }

    if (scored.empty()) return {};

    // Rank by z_statistic
    std::vector<double> z_vals;
    for (auto& s : scored) z_vals.push_back(s.z_stat);
    std::vector<double> ranks = computeRanks(z_vals);
    double max_rank = *std::max_element(ranks.begin(), ranks.end());

    // freq_score = rank / max_rank
    for (size_t i = 0; i < scored.size(); i++) {
        scored[i].z_stat = ranks[i] / max_rank; // reuse as freq_score temp
    }

    // Per gene: filter transcripts with >= 50% of max in_frame reads, return max freq_score
    std::map<std::string, std::vector<size_t>> by_gene;
    for (size_t i = 0; i < scored.size(); i++) {
        by_gene[scored[i].gene].push_back(i);
    }

    std::vector<InFrameResult> results;
    for (auto& [gene, indices] : by_gene) {
        int max_in_frame = 0;
        for (size_t idx : indices) {
            max_in_frame = std::max(max_in_frame, scored[idx].in_frame_s);
        }

        double max_freq = 0;
        std::vector<std::string> best_transcripts;

        for (size_t idx : indices) {
            // Filter: >= 50% of max in-frame reads
            if (scored[idx].in_frame_s < max_in_frame / 2) continue;

            double freq = scored[idx].z_stat; // stored freq_score above
            if (freq > max_freq) {
                max_freq = freq;
                best_transcripts.clear();
                best_transcripts.push_back(scored[idx].transcript_id);
            } else if (freq == max_freq) {
                best_transcripts.push_back(scored[idx].transcript_id);
            }
        }

        std::string transcript_str;
        for (size_t i = 0; i < best_transcripts.size(); i++) {
            if (i > 0) transcript_str += ",";
            transcript_str += best_transcripts[i];
        }

        InFrameResult r;
        r.gene = gene;
        r.bait = bait_name;
        r.freq_score = max_freq;
        r.transcripts = transcript_str;
        results.push_back(r);
    }

    return results;
}

// ============================================================================
// BordaAggregator
// ============================================================================
std::vector<Y2HScore> BordaAggregator::aggregate(
    const std::vector<EnrichmentResult>& enrichment,
    const std::vector<SpecificityResult>& specificity,
    const std::vector<InFrameResult>& in_frame)
{
    // Collect all gene-bait pairs
    std::map<std::string, Y2HScore> scores;

    for (const auto& e : enrichment) {
        std::string key = e.gene + "\t" + e.bait;
        scores[key].gene = e.gene;
        scores[key].bait = e.bait;
        scores[key].enrichment_score = e.total_score;
    }
    for (const auto& s : specificity) {
        std::string key = s.gene + "\t" + s.bait;
        scores[key].gene = s.gene;
        scores[key].bait = s.bait;
        scores[key].specificity_score = s.total_score;
    }
    for (const auto& f : in_frame) {
        std::string key = f.gene + "\t" + f.bait;
        scores[key].gene = f.gene;
        scores[key].bait = f.bait;
        scores[key].in_frame_score = f.freq_score;
        scores[key].in_frame_transcripts = f.transcripts;
    }

    // Build result vector
    std::vector<Y2HScore> result;
    for (auto& [key, s] : scores) {
        s.sum_scores = s.enrichment_score + s.specificity_score + s.in_frame_score;
        result.push_back(s);
    }

    int n = result.size();
    if (n == 0) return result;

    // Rank each score independently (higher score = rank 1)
    auto rankDescending = [](const std::vector<double>& vals) {
        std::vector<double> neg_vals(vals.size());
        for (size_t i = 0; i < vals.size(); i++) neg_vals[i] = -vals[i];
        return computeRanks(neg_vals);
    };

    std::vector<double> e_vals, s_vals, f_vals;
    for (auto& r : result) {
        e_vals.push_back(r.enrichment_score);
        s_vals.push_back(r.specificity_score);
        f_vals.push_back(r.in_frame_score);
    }

    auto e_ranks = rankDescending(e_vals);
    auto s_ranks = rankDescending(s_vals);
    auto f_ranks = rankDescending(f_vals);

    // Count how many score types are available
    bool has_e = !enrichment.empty();
    bool has_s = !specificity.empty();
    bool has_f = !in_frame.empty();
    int n_metrics = (has_e ? 1 : 0) + (has_s ? 1 : 0) + (has_f ? 1 : 0);
    if (n_metrics == 0) return result;

    // Borda: average rank, then score = n / avg_rank
    for (int i = 0; i < n; i++) {
        double sum_ranks = 0;
        if (has_e) sum_ranks += e_ranks[i];
        if (has_s) sum_ranks += s_ranks[i];
        if (has_f) sum_ranks += f_ranks[i];
        double avg_rank = sum_ranks / n_metrics;
        result[i].borda_score = (avg_rank > 0) ? n / avg_rank : 0;
    }

    // Sort by borda descending
    std::sort(result.begin(), result.end(),
        [](const Y2HScore& a, const Y2HScore& b) {
            return a.borda_score > b.borda_score;
        });

    return result;
}

} // namespace deseq2

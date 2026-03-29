#include "y2h_scores.h"
#include "deseq_dataset.h"
#include "deseq_stats.h"
#include <iostream>
#include <cmath>
#include <random>
#include <cassert>

using namespace deseq2;

static int passed = 0, failed = 0;

void check(bool cond, const std::string& name) {
    if (cond) { passed++; std::cout << "  PASS: " << name << std::endl; }
    else { failed++; std::cout << "  FAIL: " << name << std::endl; }
}

// ============================================================================
// Test 1: Enrichment Score
// ============================================================================
void testEnrichmentScore() {
    std::cout << "\n=== Test 1: Enrichment Score ===" << std::endl;

    // Create synthetic DESeq2 results (6 columns: baseMean, lfc, se, stat, pval, padj)
    // 20 genes: 5 strongly enriched, 5 moderately enriched, 10 not enriched
    int n = 20;
    Eigen::MatrixXd results(n, 6);
    std::vector<std::string> genes;

    for (int i = 0; i < n; i++) {
        genes.push_back("gene_" + std::to_string(i));
        if (i < 5) {
            // Strongly enriched: high stat, low pvalue, high LFC
            results(i, 0) = 1000;        // baseMean
            results(i, 1) = 4.0 - i*0.2; // log2FC
            results(i, 2) = 0.3;         // SE
            results(i, 3) = 10.0 - i;    // stat
            results(i, 4) = 1e-10;       // pvalue
            results(i, 5) = 1e-8;        // padj
        } else if (i < 10) {
            // Moderately enriched
            results(i, 0) = 500;
            results(i, 1) = 2.0 - (i-5)*0.2;
            results(i, 2) = 0.5;
            results(i, 3) = 3.0 - (i-5)*0.3;
            results(i, 4) = 0.01;
            results(i, 5) = 0.05;
        } else {
            // Not enriched (negative LFC or high pvalue)
            results(i, 0) = 200;
            results(i, 1) = -0.5;
            results(i, 2) = 1.0;
            results(i, 3) = -0.5;
            results(i, 4) = 0.8;
            results(i, 5) = 0.95;
        }
    }

    EnrichmentScorer scorer;

    // With default thresholds (p<1, fc>0): should include all positive LFC genes
    auto scores = scorer.compute(results, genes, "bait1", 1.0, 0.0);
    check(scores.size() == 10, "10 genes with positive LFC pass filter (got " +
          std::to_string(scores.size()) + ")");

    // With strict threshold
    auto strict = scorer.compute(results, genes, "bait1", 0.05, 1.0);
    check(strict.size() == 10, "10 genes pass p<0.05 and lfc>1.0 (got " +
          std::to_string(strict.size()) + ")");

    // Highest stat gene should have highest score
    if (!scores.empty()) {
        double max_score = 0;
        std::string best_gene;
        for (auto& s : scores) {
            if (s.total_score > max_score) {
                max_score = s.total_score;
                best_gene = s.gene;
            }
        }
        check(best_gene == "gene_0",
              "Highest Wald stat gene has highest enrichment score (got " + best_gene + ")");
    }

    // Scores should be in [0, ~1+epsilon] range
    bool valid_range = true;
    for (auto& s : scores) {
        if (s.total_score < 0 || s.total_score > 2.0) valid_range = false;
    }
    check(valid_range, "All enrichment scores in valid range [0, 2)");

    // rank_score should be (max_rank - rank) / max_rank, so in [0, 1)
    bool valid_ranks = true;
    for (auto& s : scores) {
        if (s.rank_score < 0 || s.rank_score > 1.0) valid_ranks = false;
    }
    check(valid_ranks, "All rank scores in [0, 1]");
}

// ============================================================================
// Test 2: In-Frame Score
// ============================================================================
void testInFrameScore() {
    std::cout << "\n=== Test 2: In-Frame Score ===" << std::endl;

    std::vector<JunctionData> data;

    // Gene A: 3 replicates, strongly in-frame in selected, random in non-selected
    for (int r = 0; r < 3; r++) {
        data.push_back({"geneA.1", "geneA", 90, 100, 33, 100});
    }

    // Gene B: in-frame in both (no enrichment)
    for (int r = 0; r < 3; r++) {
        data.push_back({"geneB.1", "geneB", 50, 100, 50, 100});
    }

    // Gene C: strongly in-frame selected, no reads non-selected
    for (int r = 0; r < 3; r++) {
        data.push_back({"geneC.1", "geneC", 80, 100, 0, 0});
    }

    // Gene D: not in-frame at all
    for (int r = 0; r < 3; r++) {
        data.push_back({"geneD.1", "geneD", 10, 100, 10, 100});
    }

    InFrameScorer scorer;
    auto results = scorer.compute(data, "bait1");

    check(results.size() >= 2, "At least 2 genes scored (got " +
          std::to_string(results.size()) + ")");

    // Gene A should have highest score (90% in-frame selected vs 33% non-selected)
    double geneA_score = 0, geneB_score = 0;
    for (auto& r : results) {
        if (r.gene == "geneA") geneA_score = r.freq_score;
        if (r.gene == "geneB") geneB_score = r.freq_score;
    }
    check(geneA_score > geneB_score,
          "Gene with higher in-frame enrichment scores higher (" +
          std::to_string(geneA_score) + " > " + std::to_string(geneB_score) + ")");

    // All scores should be in [0, 1]
    bool valid = true;
    for (auto& r : results) {
        if (r.freq_score < 0 || r.freq_score > 1.0) valid = false;
    }
    check(valid, "All in-frame scores in [0, 1]");
}

// ============================================================================
// Test 3: Borda Aggregation
// ============================================================================
void testBordaAggregation() {
    std::cout << "\n=== Test 3: Borda Aggregation ===" << std::endl;

    // Create 5 genes with known scores
    std::vector<EnrichmentResult> enrichment;
    std::vector<SpecificityResult> specificity;
    std::vector<InFrameResult> in_frame;

    // Gene A: best in all three → should have highest Borda
    enrichment.push_back({"geneA", "bait1", 10, 1e-10, 4, 0.95, 0.04, 0.99});
    specificity.push_back({"geneA", "bait1", 0.9, 0.05, 0.95});
    in_frame.push_back({"geneA", "bait1", "t1", 0.9, 0.33, 5.0, 0.99});

    // Gene B: second in all three
    enrichment.push_back({"geneB", "bait1", 8, 1e-8, 3, 0.8, 0.03, 0.83});
    specificity.push_back({"geneB", "bait1", 0.7, 0.04, 0.74});
    in_frame.push_back({"geneB", "bait1", "t2", 0.8, 0.33, 3.0, 0.8});

    // Gene C: good enrichment, bad specificity, medium in-frame
    enrichment.push_back({"geneC", "bait1", 7, 1e-6, 2.5, 0.7, 0.02, 0.72});
    specificity.push_back({"geneC", "bait1", 0.1, 0.01, 0.11});
    in_frame.push_back({"geneC", "bait1", "t3", 0.6, 0.33, 1.0, 0.5});

    // Gene D: only in enrichment
    enrichment.push_back({"geneD", "bait1", 5, 0.01, 1.5, 0.5, 0.01, 0.51});

    // Gene E: only in in-frame
    in_frame.push_back({"geneE", "bait1", "t5", 0.95, 0.33, 6.0, 0.6});

    BordaAggregator aggregator;
    auto scores = aggregator.aggregate(enrichment, specificity, in_frame);

    check(scores.size() == 5, "All 5 gene-bait pairs in output (got " +
          std::to_string(scores.size()) + ")");

    // Gene A should have highest Borda (best in all three)
    double geneA_borda = 0, geneB_borda = 0;
    for (auto& s : scores) {
        if (s.gene == "geneA") geneA_borda = s.borda_score;
        if (s.gene == "geneB") geneB_borda = s.borda_score;
    }
    check(geneA_borda > geneB_borda,
          "Gene best in all 3 metrics has highest Borda (" +
          std::to_string(geneA_borda) + " > " + std::to_string(geneB_borda) + ")");

    // Sum scores should be correct
    for (auto& s : scores) {
        double expected_sum = s.enrichment_score + s.specificity_score + s.in_frame_score;
        check(std::abs(s.sum_scores - expected_sum) < 1e-10,
              s.gene + ": sum_scores correct (" + std::to_string(s.sum_scores) + ")");
    }

    // All Borda scores should be positive
    bool all_positive = true;
    for (auto& s : scores) {
        if (s.borda_score <= 0) all_positive = false;
    }
    check(all_positive, "All Borda scores positive");

    // Gene with only one metric should have lower Borda than consistent genes
    double geneD_borda = 0;
    for (auto& s : scores) {
        if (s.gene == "geneD") geneD_borda = s.borda_score;
    }
    check(geneA_borda > geneD_borda,
          "Consistent gene outranks single-metric gene in Borda (" +
          std::to_string(geneA_borda) + " > " + std::to_string(geneD_borda) + ")");
}

// ============================================================================
// Test 4: Full Pipeline with DESeq2
// ============================================================================
void testFullPipeline() {
    std::cout << "\n=== Test 4: Full DESeq2 → Y2H-SCORES Pipeline ===" << std::endl;

    // Create synthetic Y2H data: 2 conditions (selected vs non-selected), 3 reps each
    // 50 genes: 10 true interactors (enriched in selection), 40 background
    std::mt19937 rng(42);
    const int nGenes = 50;
    const int nPerGroup = 3;

    Eigen::MatrixXd counts(nPerGroup * 2, nGenes);
    std::vector<std::string> gene_names;

    for (int g = 0; g < nGenes; g++) {
        gene_names.push_back("prey_" + std::to_string(g));
        double baseMu = 100 + g * 5;
        double fold = (g < 10) ? 5.0 : 1.0; // 10 true interactors at 5x enrichment

        for (int s = 0; s < nPerGroup; s++) {
            std::poisson_distribution<int> dN(baseMu);
            std::poisson_distribution<int> dS(baseMu * fold);
            counts(s, g) = std::max(1, dN(rng));              // non-selected
            counts(s + nPerGroup, g) = std::max(1, dS(rng));  // selected
        }
    }

    Eigen::MatrixXd metadata(nPerGroup * 2, 1);
    for (int i = 0; i < nPerGroup; i++) metadata(i, 0) = 0; // non-selected
    for (int i = nPerGroup; i < nPerGroup * 2; i++) metadata(i, 0) = 1; // selected

    // Run DESeq2
    DeseqDataSet dds(counts, metadata);
    dds.fitSizeFactors();
    dds.fitGenewiseDispersions();
    dds.fitDispersionTrend();
    dds.fitDispersionPrior();
    dds.fitMAPDispersions();
    dds.fitLFC();
    dds.calculateCooks();
    dds.refit();

    Eigen::VectorXd contrast(2);
    contrast << 0.0, 1.0;
    DeseqStats ds(dds, contrast, 0.05, true, true);
    ds.runWaldTest();
    ds.cooksFiltering();
    ds.independentFiltering();
    Eigen::MatrixXd deseq_results = ds.summary();

    // Compute enrichment score
    EnrichmentScorer escorer;
    auto enrichment = escorer.compute(deseq_results, gene_names, "testBait", 1.0, 0.0);

    check(!enrichment.empty(), "Enrichment scores computed (" +
          std::to_string(enrichment.size()) + " genes)");

    // True interactors should dominate top scores
    int top10_true = 0;
    std::sort(enrichment.begin(), enrichment.end(),
        [](const EnrichmentResult& a, const EnrichmentResult& b) {
            return a.total_score > b.total_score;
        });

    int check_top = std::min(10, static_cast<int>(enrichment.size()));
    for (int i = 0; i < check_top; i++) {
        // True interactors are prey_0 through prey_9
        int gene_num = std::stoi(enrichment[i].gene.substr(5));
        if (gene_num < 10) top10_true++;
    }
    check(top10_true >= 7,
          ">=70% of top 10 enrichment scores are true interactors (" +
          std::to_string(top10_true) + "/10)");

    // Create fake junction data for in-frame scoring
    std::vector<JunctionData> junctions;
    for (int g = 0; g < nGenes; g++) {
        int in_frame_s = (g < 10) ? 80 + rng() % 20 : 30 + rng() % 10;
        int total_s = 100;
        int in_frame_ns = 33;
        int total_ns = 100;
        for (int r = 0; r < 3; r++) {
            junctions.push_back({
                gene_names[g] + ".1", gene_names[g],
                in_frame_s, total_s, in_frame_ns, total_ns
            });
        }
    }

    InFrameScorer ifscorer;
    auto in_frame = ifscorer.compute(junctions, "testBait");

    check(!in_frame.empty(), "In-frame scores computed (" +
          std::to_string(in_frame.size()) + " genes)");

    // Borda aggregation (no specificity in this test)
    BordaAggregator borda;
    auto y2h_scores = borda.aggregate(enrichment, {}, in_frame);

    check(!y2h_scores.empty(), "Y2H-SCORES computed (" +
          std::to_string(y2h_scores.size()) + " gene-bait pairs)");

    // Top Borda genes should be true interactors
    int top10_borda_true = 0;
    for (int i = 0; i < std::min(10, static_cast<int>(y2h_scores.size())); i++) {
        int gene_num = std::stoi(y2h_scores[i].gene.substr(5));
        if (gene_num < 10) top10_borda_true++;
    }
    check(top10_borda_true >= 7,
          ">=70% of top 10 Borda scores are true interactors (" +
          std::to_string(top10_borda_true) + "/10)");
}

// ============================================================================
// Test 5: Edge Cases
// ============================================================================
void testEdgeCases() {
    std::cout << "\n=== Test 5: Edge Cases ===" << std::endl;

    // Empty input
    EnrichmentScorer escorer;
    Eigen::MatrixXd empty_results(0, 6);
    auto empty = escorer.compute(empty_results, {}, "bait", 1.0, 0.0);
    check(empty.empty(), "Empty DESeq2 results → empty enrichment scores");

    // All NaN input
    Eigen::MatrixXd nan_results(5, 6);
    nan_results.setConstant(std::nan(""));
    std::vector<std::string> genes = {"a", "b", "c", "d", "e"};
    auto nan_scores = escorer.compute(nan_results, genes, "bait", 1.0, 0.0);
    check(nan_scores.empty(), "All-NaN DESeq2 results → empty enrichment scores");

    // Single gene
    Eigen::MatrixXd single(1, 6);
    single << 500, 3.0, 0.5, 6.0, 1e-8, 1e-6;
    auto single_score = escorer.compute(single, {"gene1"}, "bait", 1.0, 0.0);
    check(single_score.size() == 1, "Single gene scores correctly");
    if (!single_score.empty()) {
        check(single_score[0].rank_score == 0.0, "Single gene rank_score = 0");
    }

    // In-frame with no non-selected data
    std::vector<JunctionData> no_ns;
    no_ns.push_back({"t1", "g1", 80, 100, 0, 0});
    no_ns.push_back({"t1", "g1", 85, 100, 0, 0});
    InFrameScorer ifscorer;
    auto no_ns_result = ifscorer.compute(no_ns, "bait");
    // Should still produce a result using 1/3 baseline
    check(!no_ns_result.empty(), "In-frame works with no non-selected data");
}

// ============================================================================
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << " Y2H-SCORES Validation Suite" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Reference: Velásquez-Zapata et al. (2021)" << std::endl;
    std::cout << "PLoS Computational Biology 17(4):e1008890" << std::endl;

    testEnrichmentScore();
    testInFrameScore();
    testBordaAggregation();
    testFullPipeline();
    testEdgeCases();

    std::cout << "\n========================================" << std::endl;
    std::cout << "TOTAL: " << passed << " passed, " << failed << " failed" << std::endl;
    std::cout << "========================================" << std::endl;

    if (failed > 0) {
        std::cout << "*** " << failed << " TESTS FAILED ***" << std::endl;
        return 1;
    }
    std::cout << "*** ALL TESTS PASSED ***" << std::endl;
    return 0;
}

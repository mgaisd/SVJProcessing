/**
 * create_scouting_jec_map.cpp
 * 
 * Fast C++ implementation for creating scouting JEC correction maps.
 * Uses ROOT TTree for I/O and explicit event loops for unique jet matching.
 * 
 * Performance: ~50-100x faster than Python implementation
 * 
 * Outputs:
 *   - ROOT file with TH2D histograms (sum, sum_squared, counts) per bin
 *   - JSON file with correction factors and uncertainties
 * 
 * Author: SVJ Scouting Analysis
 * Date: February 2026
 */

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2D.h"
#include "TError.h"
#include "TObjArray.h"

// JSON output helper
#include <sstream>

// Jet matching result
struct JetMatch {
    float scout_pt;
    float scout_eta;
    float offline_pt;
    float deltaR;
};

// Configuration
struct Config {
    std::vector<std::string> input_files;
    std::string output_dir;
    std::string year;
    std::string jet_type;  // "Jet" or "FatJet"
    float dr_threshold;
    float pt_min;
    int max_jets;  // -1 for all jets
    
    // Binning
    int n_pt_bins;
    float pt_min_bin, pt_max_bin;
    int n_eta_bins;
    float eta_min_bin, eta_max_bin;
};

// Calculate deltaR
inline float deltaR(float eta1, float phi1, float eta2, float phi2) {
    float deta = eta1 - eta2;
    float dphi = phi1 - phi2;
    // Wrap phi to [-pi, pi]
    while (dphi > M_PI) dphi -= 2*M_PI;
    while (dphi < -M_PI) dphi += 2*M_PI;
    return std::sqrt(deta*deta + dphi*dphi);
}

// Create directory (recursive)
void createDirectory(const std::string& path) {
    // Use system mkdir -p for simplicity
    std::string cmd = "mkdir -p " + path;
    int ret = system(cmd.c_str());
    if (ret != 0) {
        std::cerr << "Warning: Could not create directory " << path << std::endl;
    }
}

// Match jets with unique matching (each offline jet matched to at most one scout)
std::vector<JetMatch> matchJets(
    const std::vector<float>& scout_pt,
    const std::vector<float>& scout_eta,
    const std::vector<float>& scout_phi,
    const std::vector<float>& offline_pt,
    const std::vector<float>& offline_eta,
    const std::vector<float>& offline_phi,
    float dr_threshold,
    int& n_scouts,
    int& n_matched)
{
    std::vector<JetMatch> matches;
    
    int n_scout = scout_pt.size();
    int n_offline = offline_pt.size();
    
    n_scouts += n_scout;
    
    if (n_scout == 0 || n_offline == 0) {
        return matches;
    }
    
    // Store all candidate matches: (scout_idx, offline_idx, deltaR)
    struct Candidate {
        int scout_idx;
        int offline_idx;
        float dr;
    };
    std::vector<Candidate> candidates;
    
    // Find best match for each scout jet
    for (int i_scout = 0; i_scout < n_scout; ++i_scout) {
        float best_dr = dr_threshold;
        int best_offline_idx = -1;
        
        for (int i_offline = 0; i_offline < n_offline; ++i_offline) {
            float dr = deltaR(scout_eta[i_scout], scout_phi[i_scout],
                            offline_eta[i_offline], offline_phi[i_offline]);
            
            if (dr < best_dr) {
                best_dr = dr;
                best_offline_idx = i_offline;
            }
        }
        
        if (best_offline_idx >= 0) {
            candidates.push_back({i_scout, best_offline_idx, best_dr});
        }
    }
    
    if (candidates.empty()) {
        return matches;
    }
    
    // Resolve conflicts: sort by offline index, then by deltaR
    std::sort(candidates.begin(), candidates.end(),
        [](const Candidate& a, const Candidate& b) {
            if (a.offline_idx != b.offline_idx) return a.offline_idx < b.offline_idx;
            return a.dr < b.dr;
        });
    
    // Keep only first occurrence of each unique offline jet (smallest deltaR)
    std::vector<bool> offline_used(n_offline, false);
    
    for (const auto& cand : candidates) {
        if (!offline_used[cand.offline_idx]) {
            offline_used[cand.offline_idx] = true;
            
            // Sanity checks
            float s_pt = scout_pt[cand.scout_idx];
            float o_pt = offline_pt[cand.offline_idx];
            
            if (s_pt > 0 && o_pt > 0 && std::isfinite(s_pt) && std::isfinite(o_pt)) {
                JetMatch match;
                match.scout_pt = s_pt;
                match.scout_eta = scout_eta[cand.scout_idx];
                match.offline_pt = o_pt;
                match.deltaR = cand.dr;
                matches.push_back(match);
                n_matched++;
            }
        }
    }
    
    return matches;
}

// Process one ROOT file for BOTH jet types simultaneously
void processFileBoth(
    const std::string& filename,
    const Config& config_jet, TH2D* h_sum_jet, TH2D* h_sum_sq_jet, TH2D* h_counts_jet,
    long& n_total_jets, long& n_rejected_neg_jet, long& n_rejected_outliers_jet,
    const Config& config_fatjet, TH2D* h_sum_fatjet, TH2D* h_sum_sq_fatjet, TH2D* h_counts_fatjet,
    long& n_total_fatjets, long& n_rejected_neg_fatjet, long& n_rejected_outliers_fatjet)
{
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "  Warning: Cannot open file: " << filename << std::endl;
        return;
    }
    
    TTree* tree = (TTree*)file->Get("mmtree/Events");
    if (!tree) {
        std::cerr << "  Warning: Cannot find tree mmtree/Events in " << filename << std::endl;
        file->Close();
        return;
    }
    
    // Enable TTreeCache for faster reading (30 MB cache)
    tree->SetCacheSize(30000000);
    tree->AddBranchToCache("*", kTRUE);
    tree->StopCacheLearningPhase();
    
    // Set up branches for BOTH jet types
    std::vector<float> *jet_pt = nullptr, *jet_eta = nullptr, *jet_phi = nullptr;
    std::vector<float> *offline_jet_pt = nullptr, *offline_jet_eta = nullptr, *offline_jet_phi = nullptr;
    std::vector<float> *fatjet_pt = nullptr, *fatjet_eta = nullptr, *fatjet_phi = nullptr;
    std::vector<float> *offline_fatjet_pt = nullptr, *offline_fatjet_eta = nullptr, *offline_fatjet_phi = nullptr;
    
    tree->SetBranchAddress("Jet_pt", &jet_pt);
    tree->SetBranchAddress("Jet_eta", &jet_eta);
    tree->SetBranchAddress("Jet_phi", &jet_phi);
    tree->SetBranchAddress("OfflineJet_pt", &offline_jet_pt);
    tree->SetBranchAddress("OfflineJet_eta", &offline_jet_eta);
    tree->SetBranchAddress("OfflineJet_phi", &offline_jet_phi);
    
    tree->SetBranchAddress("FatJet_pt", &fatjet_pt);
    tree->SetBranchAddress("FatJet_eta", &fatjet_eta);
    tree->SetBranchAddress("FatJet_phi", &fatjet_phi);
    tree->SetBranchAddress("OfflineFatJet_pt", &offline_fatjet_pt);
    tree->SetBranchAddress("OfflineFatJet_eta", &offline_fatjet_eta);
    tree->SetBranchAddress("OfflineFatJet_phi", &offline_fatjet_phi);
    
    int n_events = tree->GetEntries();
    
    // Process events - fill BOTH jet types simultaneously
    for (int i_evt = 0; i_evt < n_events; ++i_evt) {
        tree->GetEntry(i_evt);
        
        // ===== Process AK4 (Jet) =====
        if (jet_pt && offline_jet_pt) {
            std::vector<float> s_pt, s_eta, s_phi;
            std::vector<float> o_pt, o_eta, o_phi;
            
            int max_scout = (config_jet.max_jets > 0) ? std::min(config_jet.max_jets, (int)jet_pt->size()) : jet_pt->size();
            
            for (int i = 0; i < max_scout; ++i) {
                if ((*jet_pt)[i] > config_jet.pt_min) {
                    s_pt.push_back((*jet_pt)[i]);
                    s_eta.push_back((*jet_eta)[i]);
                    s_phi.push_back((*jet_phi)[i]);
                }
            }
            
            for (size_t i = 0; i < offline_jet_pt->size(); ++i) {
                if ((*offline_jet_pt)[i] > config_jet.pt_min) {
                    o_pt.push_back((*offline_jet_pt)[i]);
                    o_eta.push_back((*offline_jet_eta)[i]);
                    o_phi.push_back((*offline_jet_phi)[i]);
                }
            }
            
            int n_scouts_temp = 0, n_matched_temp = 0;
            std::vector<JetMatch> matches = matchJets(s_pt, s_eta, s_phi, o_pt, o_eta, o_phi,
                                                      config_jet.dr_threshold, n_scouts_temp, n_matched_temp);
            
            for (const auto& match : matches) {
                n_total_jets++;
                float ratio = match.offline_pt / match.scout_pt;
                
                if (ratio <= 0.1 || ratio >= 10.0) {
                    n_rejected_outliers_jet++;
                    continue;
                }
                
                int bin_x = h_sum_jet->GetXaxis()->FindBin(match.scout_pt);
                int bin_y = h_sum_jet->GetYaxis()->FindBin(match.scout_eta);
                
                if (bin_x >= 1 && bin_x <= h_sum_jet->GetNbinsX() &&
                    bin_y >= 1 && bin_y <= h_sum_jet->GetNbinsY()) {
                    h_sum_jet->Fill(match.scout_pt, match.scout_eta, ratio);
                    h_sum_sq_jet->Fill(match.scout_pt, match.scout_eta, ratio * ratio);
                    h_counts_jet->Fill(match.scout_pt, match.scout_eta);
                }
            }
        }
        
        // ===== Process AK8 (FatJet) =====
        if (fatjet_pt && offline_fatjet_pt) {
            std::vector<float> s_pt, s_eta, s_phi;
            std::vector<float> o_pt, o_eta, o_phi;
            
            int max_scout = (config_fatjet.max_jets > 0) ? std::min(config_fatjet.max_jets, (int)fatjet_pt->size()) : fatjet_pt->size();
            
            for (int i = 0; i < max_scout; ++i) {
                if ((*fatjet_pt)[i] > config_fatjet.pt_min) {
                    s_pt.push_back((*fatjet_pt)[i]);
                    s_eta.push_back((*fatjet_eta)[i]);
                    s_phi.push_back((*fatjet_phi)[i]);
                }
            }
            
            for (size_t i = 0; i < offline_fatjet_pt->size(); ++i) {
                if ((*offline_fatjet_pt)[i] > config_fatjet.pt_min) {
                    o_pt.push_back((*offline_fatjet_pt)[i]);
                    o_eta.push_back((*offline_fatjet_eta)[i]);
                    o_phi.push_back((*offline_fatjet_phi)[i]);
                }
            }
            
            int n_scouts_temp = 0, n_matched_temp = 0;
            std::vector<JetMatch> matches = matchJets(s_pt, s_eta, s_phi, o_pt, o_eta, o_phi,
                                                      config_fatjet.dr_threshold, n_scouts_temp, n_matched_temp);
            
            for (const auto& match : matches) {
                n_total_fatjets++;
                float ratio = match.offline_pt / match.scout_pt;
                
                if (ratio <= 0.1 || ratio >= 10.0) {
                    n_rejected_outliers_fatjet++;
                    continue;
                }
                
                int bin_x = h_sum_fatjet->GetXaxis()->FindBin(match.scout_pt);
                int bin_y = h_sum_fatjet->GetYaxis()->FindBin(match.scout_eta);
                
                if (bin_x >= 1 && bin_x <= h_sum_fatjet->GetNbinsX() &&
                    bin_y >= 1 && bin_y <= h_sum_fatjet->GetNbinsY()) {
                    h_sum_fatjet->Fill(match.scout_pt, match.scout_eta, ratio);
                    h_sum_sq_fatjet->Fill(match.scout_pt, match.scout_eta, ratio * ratio);
                    h_counts_fatjet->Fill(match.scout_pt, match.scout_eta);
                }
            }
        }
    }
    
    file->Close();
}

// Save results as correctionlib schema v2 JSON with nested binning nodes.
//
// Output format (loadable directly with correctionlib.CorrectionSet.from_file):
//
//   {
//     "schema_version": 2,
//     "corrections": [{
//       "name": "scouting_jec_<JetType>",
//       "inputs": [{"name":"JetPt","type":"real"}, {"name":"JetEta","type":"real"}],
//       "output": {"name":"factor","type":"real"},
//       "data": {
//         "nodetype": "binning",   // outer: bins on JetPt
//         "input": "JetPt",
//         "edges": [...],
//         "flow": "clamp",
//         "content": [             // one inner node per pt bin
//           {
//             "nodetype": "binning",
//             "input": "JetEta",
//             "edges": [...],
//             "flow": "clamp",
//             "content": [1.05, 1.03, ...]   // one value per eta bin
//           },
//           ...
//         ]
//       }
//     }]
//   }
//
// Bins with fewer than min_count entries are set to 1.0 (no correction).
void saveJSON(const Config& config, TH2D* h_sum, TH2D* h_sum_sq, TH2D* h_counts,
              long n_total, long n_rejected_neg, long n_rejected_outliers)
{
    // Output filename matches the correctionlib convention used by
    // corrections_pfnano_coffea.py and skimmer_utils.py.
    std::string json_file = config.output_dir + "/scouting_jec_corrections_"
                            + config.year + "_" + config.jet_type + "_correctionlib.json";
    std::ofstream out(json_file);

    if (!out.is_open()) {
        std::cerr << "Error: Cannot create JSON file: " << json_file << std::endl;
        return;
    }

    out << std::fixed << std::setprecision(6);

    const int    min_count    = 10;
    const int    n_pt         = h_sum->GetNbinsX();
    const int    n_eta        = h_sum->GetNbinsY();
    std::string  corr_name    = "scouting_jec_" + config.jet_type;

    // ------------------------------------------------------------------
    // Pre-compute correction values so the writing loop below is clean.
    // corrections[ix][iy] corresponds to ROOT histogram bin (ix+1, iy+1).
    // Bins with insufficient stats are left at 1.0 (no correction).
    // ------------------------------------------------------------------
    std::vector<std::vector<double>> corrections(n_pt, std::vector<double>(n_eta, 1.0));
    for (int ix = 0; ix < n_pt; ++ix) {
        for (int iy = 0; iy < n_eta; ++iy) {
            double sum   = h_sum->GetBinContent(ix + 1, iy + 1);
            double count = h_counts->GetBinContent(ix + 1, iy + 1);
            if (count >= min_count) {
                double corr = sum / count;
                // Sanity clamp: treat extreme values as no-correction
                corrections[ix][iy] = (corr > 0.1 && corr < 10.0) ? corr : 1.0;
            }
        }
    }

    // ------------------------------------------------------------------
    // Helper: write a ROOT axis's bin edges as a JSON array inline.
    // ------------------------------------------------------------------
    auto write_edges = [&](TAxis* axis, int n_bins) {
        out << "[";
        for (int i = 1; i <= n_bins + 1; ++i) {
            out << axis->GetBinLowEdge(i);
            if (i <= n_bins) out << ", ";
        }
        out << "]";
    };

    // ------------------------------------------------------------------
    // Write correctionlib schema v2 JSON
    // ------------------------------------------------------------------
    out << "{\n";
    out << "  \"schema_version\": 2,\n";
    out << "  \"description\": \"Scouting JEC residual corrections for "
        << config.jet_type << ", derived from scouting/offline jet matching.\",\n";
    out << "  \"corrections\": [\n";
    out << "    {\n";
    out << "      \"name\": \"" << corr_name << "\",\n";
    out << "      \"version\": 1,\n";
    out << "      \"description\": \"Scouting JEC correction for " << config.jet_type
        << ". dR<" << config.dr_threshold
        << ", pt_min=" << config.pt_min << " GeV.\",\n";
    out << "      \"inputs\": [\n";
    out << "        {\"name\": \"JetPt\",  \"type\": \"real\","
        << " \"description\": \"Scouting jet pT (GeV)\"},\n";
    out << "        {\"name\": \"JetEta\", \"type\": \"real\","
        << " \"description\": \"Scouting jet eta\"}\n";
    out << "      ],\n";
    out << "      \"output\": {\"name\": \"factor\", \"type\": \"real\","
        << " \"description\": \"JEC scale factor (offline/scouting)\"},\n";

    // Outer binning node: bins on JetPt
    out << "      \"data\": {\n";
    out << "        \"nodetype\": \"binning\",\n";
    out << "        \"input\": \"JetPt\",\n";
    out << "        \"edges\": ";
    write_edges(h_sum->GetXaxis(), n_pt);
    out << ",\n";
    out << "        \"flow\": \"clamp\",\n";

    // content: one inner binning node per pt bin, each binning on JetEta
    out << "        \"content\": [\n";
    for (int ix = 0; ix < n_pt; ++ix) {
        out << "          {\n";
        out << "            \"nodetype\": \"binning\",\n";
        out << "            \"input\": \"JetEta\",\n";
        out << "            \"edges\": ";
        write_edges(h_sum->GetYaxis(), n_eta);
        out << ",\n";
        out << "            \"flow\": \"clamp\",\n";
        out << "            \"content\": [";
        for (int iy = 0; iy < n_eta; ++iy) {
            out << corrections[ix][iy];
            if (iy < n_eta - 1) out << ", ";
        }
        out << "]\n";
        out << "          }";
        if (ix < n_pt - 1) out << ",";
        out << "\n";
    }
    out << "        ]\n";   // end content
    out << "      }\n";     // end data
    out << "    }\n";       // end correction object
    out << "  ]\n";         // end corrections array
    out << "}\n";

    out.close();

    long n_accepted = n_total - n_rejected_neg - n_rejected_outliers;
    std::cout << "Saved correctionlib JSON: " << json_file << "\n";
    std::cout << "  Correction name: " << corr_name << "\n";
    std::cout << "  Bins: " << n_pt << " pt x " << n_eta << " eta\n";
    std::cout << "  Jets accepted: " << n_accepted << " / " << n_total << "\n";

    // ------------------------------------------------------------------
    // Write flat diagnostics JSON for the plotting script
    // ------------------------------------------------------------------
    std::string diag_file = config.output_dir + "/scouting_jec_diagnostics_"
                            + config.year + "_" + config.jet_type + ".json";
    std::ofstream diag(diag_file);
    if (!diag.is_open()) {
        std::cerr << "Warning: Cannot create diagnostics file: " << diag_file << "\n";
        return;
    }
    diag << std::fixed << std::setprecision(6);

    // Collect bin edges
    std::vector<double> pt_edges, eta_edges;
    for (int i = 1; i <= n_pt + 1; ++i) pt_edges.push_back(h_sum->GetXaxis()->GetBinLowEdge(i));
    for (int i = 1; i <= n_eta + 1; ++i) eta_edges.push_back(h_sum->GetYaxis()->GetBinLowEdge(i));

    // Compute uncertainties and stat_uncertainties
    std::vector<std::vector<double>> uncertainties(n_pt, std::vector<double>(n_eta, 0.0));
    std::vector<std::vector<double>> stat_uncertainties(n_pt, std::vector<double>(n_eta, 0.0));
    std::vector<std::vector<double>> counts_arr(n_pt, std::vector<double>(n_eta, 0.0));
    for (int ix = 0; ix < n_pt; ++ix) {
        for (int iy = 0; iy < n_eta; ++iy) {
            double sum    = h_sum->GetBinContent(ix + 1, iy + 1);
            double sum_sq = h_sum_sq->GetBinContent(ix + 1, iy + 1);
            double count  = h_counts->GetBinContent(ix + 1, iy + 1);
            counts_arr[ix][iy] = count;
            if (count >= 2) {
                double mean   = sum / count;
                double var    = (sum_sq / count) - mean * mean;
                double stddev = (var > 0) ? std::sqrt(var) : 0.0;
                uncertainties[ix][iy]      = stddev;
                stat_uncertainties[ix][iy] = stddev / std::sqrt(count);
            }
        }
    }

    auto write_1d = [&](const std::vector<double>& v) {
        diag << "[";
        for (size_t k = 0; k < v.size(); ++k) {
            diag << v[k];
            if (k + 1 < v.size()) diag << ", ";
        }
        diag << "]";
    };
    auto write_2d = [&](const std::vector<std::vector<double>>& arr) {
        diag << "[\n";
        for (int ix = 0; ix < n_pt; ++ix) {
            diag << "    [";
            for (int iy = 0; iy < n_eta; ++iy) {
                diag << arr[ix][iy];
                if (iy < n_eta - 1) diag << ", ";
            }
            diag << "]";
            if (ix < n_pt - 1) diag << ",";
            diag << "\n";
        }
        diag << "  ]";
    };

    diag << "{\n";
    diag << "  \"jet_type\": \"" << config.jet_type << "\",\n";
    diag << "  \"pt_bins\": ";  write_1d(pt_edges);  diag << ",\n";
    diag << "  \"eta_bins\": "; write_1d(eta_edges); diag << ",\n";
    diag << "  \"corrections\": ";        write_2d(corrections);        diag << ",\n";
    diag << "  \"uncertainties\": ";      write_2d(uncertainties);      diag << ",\n";
    diag << "  \"stat_uncertainties\": "; write_2d(stat_uncertainties); diag << ",\n";
    diag << "  \"counts\": ";             write_2d(counts_arr);         diag << ",\n";
    diag << "  \"statistics\": {\n";
    diag << "    \"total_matched\": "        << n_total            << ",\n";
    diag << "    \"rejected_negative_pt\": " << n_rejected_neg     << ",\n";
    diag << "    \"rejected_outliers\": "    << n_rejected_outliers << ",\n";
    diag << "    \"accepted\": "             << n_accepted          << "\n";
    diag << "  }\n";
    diag << "}\n";
    diag.close();
    std::cout << "Saved diagnostics JSON: " << diag_file << "\n";
}

int main(int argc, char** argv) {
    // Suppress ROOT warnings about TClassTable
    gErrorIgnoreLevel = kError;
    
    std::cout << "=============================================================\n";
    std::cout << "Scouting JEC Map Creator (C++ version)\n";
    std::cout << "=============================================================\n";
    
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " [--output-dir <dir>] <file1.root> [file2.root ...]\n";
        std::cerr << "       " << argv[0] << " [--output-dir <dir>] --filelist <filelist.txt>\n";
        std::cerr << "  Processes both Jet (AK4) and FatJet (AK8) in a single pass\n";
        return 1;
    }
    
    // Parse arguments
    std::string output_dir = "test";
    std::string year = "";
    std::vector<std::string> input_files;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--output-dir") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --output-dir requires a directory argument\n";
                return 1;
            }
            output_dir = argv[++i];
        } else if (std::string(argv[i]) == "--year") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --year requires an argument\n";
                return 1;
            }
            year = argv[++i];
        } else if (std::string(argv[i]) == "--filelist") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --filelist requires a file path argument\n";
                return 1;
            }
            std::ifstream filelist(argv[++i]);
            if (!filelist.is_open()) {
                std::cerr << "Error: Cannot open filelist: " << argv[i] << "\n";
                return 1;
            }
            std::string line;
            while (std::getline(filelist, line)) {
                if (!line.empty()) input_files.push_back(line);
            }
        } else {
            input_files.push_back(argv[i]);
        }
    }
    
    if (year.empty()) {
        std::cerr << "Warning: --year not specified, output files will not include year\n";
    }

    // Configuration for AK4 (Jet)
    Config config_jet;
    config_jet.jet_type = "Jet";
    config_jet.year = year;
    config_jet.dr_threshold = 0.15;
    config_jet.output_dir = output_dir;
    config_jet.pt_min = 15.0;
    config_jet.pt_min_bin = 15.0;
    config_jet.pt_max_bin = 2000.0;
    config_jet.n_pt_bins = 30;
    config_jet.eta_min_bin = -5.0;
    config_jet.eta_max_bin = 5.0;
    config_jet.n_eta_bins = 24;
    config_jet.max_jets = -1;  // All jets
    config_jet.input_files = input_files;
    
    // Configuration for AK8 (FatJet)
    Config config_fatjet;
    config_fatjet.jet_type = "FatJet";
    config_fatjet.year = year;
    config_fatjet.dr_threshold = 0.15;
    config_fatjet.output_dir = output_dir;
    config_fatjet.pt_min = 150.0;
    config_fatjet.pt_min_bin = 150.0;
    config_fatjet.pt_max_bin = 2000.0;
    config_fatjet.n_pt_bins = 30;
    config_fatjet.eta_min_bin = -2.4;
    config_fatjet.eta_max_bin = 2.4;
    config_fatjet.n_eta_bins = 24;
    config_fatjet.max_jets = 2;  // Leading 2 jets
    config_fatjet.input_files = input_files;
    
    std::cout << "Configuration:\n";
    std::cout << "  Processing: Jet (AK4) + FatJet (AK8) simultaneously\n";
    std::cout << "  deltaR threshold: " << config_jet.dr_threshold << "\n";
    std::cout << "  Input files: " << input_files.size() << "\n";
    std::cout << "  Output directory: " << config_jet.output_dir << "\n\n";
    
    // Create output directory
    createDirectory(config_jet.output_dir);
    
    // Create histograms for AK4 (Jet)
    TH2D* h_sum_jet = new TH2D("h_sum_jet", "Sum of ratios (Jet)",
                          config_jet.n_pt_bins, config_jet.pt_min_bin, config_jet.pt_max_bin,
                          config_jet.n_eta_bins, config_jet.eta_min_bin, config_jet.eta_max_bin);
    TH2D* h_sum_sq_jet = new TH2D("h_sum_sq_jet", "Sum of ratios squared (Jet)",
                             config_jet.n_pt_bins, config_jet.pt_min_bin, config_jet.pt_max_bin,
                             config_jet.n_eta_bins, config_jet.eta_min_bin, config_jet.eta_max_bin);
    TH2D* h_counts_jet = new TH2D("h_counts_jet", "Number of jets (Jet)",
                             config_jet.n_pt_bins, config_jet.pt_min_bin, config_jet.pt_max_bin,
                             config_jet.n_eta_bins, config_jet.eta_min_bin, config_jet.eta_max_bin);
    
    h_sum_jet->Sumw2();
    h_sum_sq_jet->Sumw2();
    h_counts_jet->Sumw2();
    
    // Create histograms for AK8 (FatJet)
    TH2D* h_sum_fatjet = new TH2D("h_sum_fatjet", "Sum of ratios (FatJet)",
                          config_fatjet.n_pt_bins, config_fatjet.pt_min_bin, config_fatjet.pt_max_bin,
                          config_fatjet.n_eta_bins, config_fatjet.eta_min_bin, config_fatjet.eta_max_bin);
    TH2D* h_sum_sq_fatjet = new TH2D("h_sum_sq_fatjet", "Sum of ratios squared (FatJet)",
                             config_fatjet.n_pt_bins, config_fatjet.pt_min_bin, config_fatjet.pt_max_bin,
                             config_fatjet.n_eta_bins, config_fatjet.eta_min_bin, config_fatjet.eta_max_bin);
    TH2D* h_counts_fatjet = new TH2D("h_counts_fatjet", "Number of jets (FatJet)",
                             config_fatjet.n_pt_bins, config_fatjet.pt_min_bin, config_fatjet.pt_max_bin,
                             config_fatjet.n_eta_bins, config_fatjet.eta_min_bin, config_fatjet.eta_max_bin);
    
    h_sum_fatjet->Sumw2();
    h_sum_sq_fatjet->Sumw2();
    h_counts_fatjet->Sumw2();
    
    // Statistics for both jet types
    long n_total_jets = 0, n_rejected_neg = 0, n_rejected_outliers = 0;
    long n_total_fatjets = 0, n_rejected_neg_fatjet = 0, n_rejected_outliers_fatjet = 0;
    
    // Process files - BOTH jet types simultaneously
    std::cout << "Processing files (both Jet and FatJet)...\n";
    for (size_t i = 0; i < input_files.size(); ++i) {
        double pct = 100.0 * (i + 1) / input_files.size();
        std::cout << "\r  File " << (i+1) << "/" << input_files.size()
                  << " (" << std::fixed << std::setprecision(1) << pct << "%)   " << std::flush;
        
        // Process BOTH jet types from same file in one pass
        processFileBoth(input_files[i],
                       config_jet, h_sum_jet, h_sum_sq_jet, h_counts_jet,
                       n_total_jets, n_rejected_neg, n_rejected_outliers,
                       config_fatjet, h_sum_fatjet, h_sum_sq_fatjet, h_counts_fatjet,
                       n_total_fatjets, n_rejected_neg_fatjet, n_rejected_outliers_fatjet);
    }
    std::cout << "\n";
    
    // Print statistics for AK4 (Jet)
    std::cout << "\n============================================\n";
    std::cout << "AK4 (Jet) Summary:\n";
    std::cout << "  Total matched jets: " << n_total_jets << "\n";
    std::cout << "  Rejected (outliers): " << n_rejected_outliers 
              << " (" << (n_total_jets > 0 ? 100.0*n_rejected_outliers/n_total_jets : 0) << "%)\n";
    std::cout << "  Accepted: " << (n_total_jets - n_rejected_neg - n_rejected_outliers)
              << " (" << (n_total_jets > 0 ? 100.0*(n_total_jets-n_rejected_neg-n_rejected_outliers)/n_total_jets : 0) << "%)\n";


    // Print statistics for AK8 (FatJet)
    std::cout << "\nAK8 (FatJet) Summary:\n";
    std::cout << "  Total matched jets: " << n_total_fatjets << "\n";
    std::cout << "  Rejected (outliers): " << n_rejected_outliers_fatjet 
              << " (" << (n_total_fatjets > 0 ? 100.0*n_rejected_outliers_fatjet/n_total_fatjets : 0) << "%)\n";
    std::cout << "  Accepted: " << (n_total_fatjets - n_rejected_neg_fatjet - n_rejected_outliers_fatjet)
              << " (" << (n_total_fatjets > 0 ? 100.0*(n_total_fatjets-n_rejected_neg_fatjet-n_rejected_outliers_fatjet)/n_total_fatjets : 0) << "%)\n";

    
    // Save ROOT histograms for AK4
    std::string root_file_jet = config_jet.output_dir + "/scouting_jec_histograms_"
                                 + config_jet.year + "_Jet.root";
    TFile* out_file_jet = new TFile(root_file_jet.c_str(), "RECREATE");
    h_sum_jet->Write("h_sum");
    h_sum_sq_jet->Write("h_sum_sq");
    h_counts_jet->Write("h_counts");
    out_file_jet->Close();
    std::cout << "\nSaved ROOT file: " << root_file_jet << std::endl;
    
    // Save ROOT histograms for AK8
    std::string root_file_fatjet = config_fatjet.output_dir + "/scouting_jec_histograms_"
                                    + config_fatjet.year + "_FatJet.root";
    TFile* out_file_fatjet = new TFile(root_file_fatjet.c_str(), "RECREATE");
    h_sum_fatjet->Write("h_sum");
    h_sum_sq_fatjet->Write("h_sum_sq");
    h_counts_fatjet->Write("h_counts");
    out_file_fatjet->Close();
    std::cout << "Saved ROOT file: " << root_file_fatjet << std::endl;
    
    // Save JSON for AK4
    saveJSON(config_jet, h_sum_jet, h_sum_sq_jet, h_counts_jet,
             n_total_jets, n_rejected_neg, n_rejected_outliers);
    
    // Save JSON for AK8
    saveJSON(config_fatjet, h_sum_fatjet, h_sum_sq_fatjet, h_counts_fatjet,
             n_total_fatjets, n_rejected_neg_fatjet, n_rejected_outliers_fatjet);
    
    // Cleanup
    delete h_sum_jet;
    delete h_sum_sq_jet;
    delete h_counts_jet;
    delete h_sum_fatjet;
    delete h_sum_sq_fatjet;
    delete h_counts_fatjet;
    
    std::cout << "\n=============================================================\n";
    std::cout << "Done! Run plot_jec_corrections.py to create plots.\n";
    std::cout << "=============================================================\n";
    
    return 0;
}

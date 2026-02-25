#!/usr/bin/env python3
"""
Script to compare corrected vs uncorrected jet pT and MET distributions
across QCD bins and signal samples with proper cross-section weighting.

JEC DIAGNOSTIC METRICS AND BEST PRACTICES:
==========================================

1. RATIO PLOTS (Reco/Gen):
   - PRIMARY DIAGNOSTIC: Compare corrected vs uncorrected ratio to Gen
   - Corrected jets should be closer to 1.0 across all pT ranges
   - Look for pT-dependent trends (corrections should flatten response)
   - Typical uncorrected jets are ~10-20% low at low pT due to detector effects
   
2. RESPONSE METRICS:
   - Mean response: <pT_reco / pT_gen> in bins of pT
   - RMS/width: spread of response (resolution)
   - Bias: deviation of mean response from 1.0
   
3. MET CLOSURE:
   - MET is very sensitive to JECs (vectorial sum of all objects)
   - Better JECs → MET_reco closer to MET_gen
   - Compare MET resolution: RMS(MET_reco - MET_gen)
   - Watch for systematic shifts in MET_reco/MET_gen ratio
   
4. PT-DEPENDENT EFFECTS:
   - Low pT (<100 GeV): largest corrections, sensitive to pileup
   - High pT (>500 GeV): smaller corrections, radiation effects
   - Plot response vs pT in bins to check correction quality
   
5. JET MATCHING CONSIDERATIONS:
   
   NO MATCHING NEEDED FOR:
   - Inclusive jet pT distributions (as done here)
   - Average response checks across full sample
   - Statistical shape comparisons
   
   MATCHING NEEDED FOR:
   - Event-by-event response calculations
   - Individual jet pT_reco/pT_gen ratios
   - Detailed jet-by-jet JEC validation
   - Resolution studies (requires deltaR < 0.4 matching typically)
   
   For this diagnostic script:
   - We look at DISTRIBUTIONS, not per-jet correlations
   - Statistical comparison of shapes is sufficient
   - Matching would be needed for per-jet response studies
   
6. ADDITIONAL USEFUL METRICS (not yet implemented):
   - Delta(pT)/pT vs eta: spatial dependence of corrections
   - Response vs number of vertices (pileup dependence)
   - Mean jet multiplicity: corrected vs uncorrected
   - Jet mass resolution (for fat jets)
   - Alpha/MPF methods for data-MC JEC validation

"""

import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import awkward as ak
from pathlib import Path
import glob

# Set CMS style
plt.style.use(hep.style.CMS)

# Base path for the datasets
BASE_PATH = "root://cmsdcache-kit-disk.gridka.de:1094//store/user/mgaisdor/SVJScouting_skims_JEC/2018/s_channel_scouting_pre_selection/nominal"

# Define the datasets
DATASETS = {
    "QCD_HT300to500": f"{BASE_PATH}/QCD_HT300to500",
    "QCD_HT500to700": f"{BASE_PATH}/QCD_HT500to700",
    "QCD_HT700to1000": f"{BASE_PATH}/QCD_HT700to1000",
    "QCD_HT1000to1500": f"{BASE_PATH}/QCD_HT1000to1500",
    "QCD_HT1500to2000": f"{BASE_PATH}/QCD_HT1500to2000",
    "QCD_HT2000toInf": f"{BASE_PATH}/QCD_HT2000toInf",
    "Signal_mMed-1100_mDark-20_rinv-0.3": f"{BASE_PATH}/s-channel_mMed-1100_mDark-20_rinv-0.3",
}

# Integrated luminosity for 2018 (fb^-1)
LUMI = 59.83  # fb^-1


def get_file_list(dataset_path):
    """Get list of ROOT files in the dataset directory via xrootd."""
    # For xrootd paths, we need to list files using uproot or subprocess
    # This assumes files follow a pattern like *.root
    import subprocess
    
    try:
        # Try to list files using xrdfs
        cmd = f"xrdfs cmsdcache-kit-disk.gridka.de:1094 ls {dataset_path.replace('root://cmsdcache-kit-disk.gridka.de:1094/', '')}"
        result = subprocess.run(cmd.split(), capture_output=True, text=True)
        files = [f"root://cmsdcache-kit-disk.gridka.de:1094/{line.strip()}" 
                 for line in result.stdout.split('\n') 
                 if line.strip().endswith('.root')]
        return files
    except Exception as e:
        print(f"Warning: Could not list files for {dataset_path}: {e}")
        return []


def read_dataset(dataset_name, dataset_path):
    """
    Read a dataset and return weighted histograms.
    
    Returns:
        dict with keys 'pt_corr', 'pt_uncorr', 'weight_sum'
    """
    print(f"\nProcessing {dataset_name}...")
    
    files = get_file_list(dataset_path)
    if not files:
        print(f"  No files found for {dataset_name}")
        return None
    
    print(f"  Found {len(files)} files")
    
    pt_corr_all = []
    pt_uncorr_all = []
    pt_gen_all = []
    met_corr_all = []
    met_uncorr_all = []
    met_gen_all = []
    weights_all = []
    
    total_initial_events = 0
    total_xsec = None
    
    for i, file_path in enumerate(files):
        try:
            print(f"  Reading file {i+1}/{len(files)}: {Path(file_path).name}")
            
            with uproot.open(file_path) as f:
                # Read metadata
                try:
                    cutflow = f["CutFlow"]
                    initial_events = cutflow["Initial"].array(library="np")[0]
                    total_initial_events += initial_events
                except Exception as e:
                    print(f"    Warning: Could not read CutFlow: {e}")
                    continue
                
                try:
                    metadata = f["Metadata"]
                    xsec = metadata["GenCrossSection"].array(library="np")[0]
                    if total_xsec is None:
                        total_xsec = xsec
                except Exception as e:
                    print(f"    Warning: Could not read GenCrossSection: {e}")
                    xsec = 1.0
                
                # Read Events tree
                try:
                    events = f["Events"]
                    pt_corr = events["FatJet_pt"].array(library="ak")
                    pt_uncorr = events["FatJet_pt_uncorr"].array(library="ak")
                    pt_gen = events["GenFatJet_pt"].array(library="ak")
                    
                    # Read MET branches
                    met_corr = events["ScoutMET_pt"].array(library="ak")
                    met_uncorr = events["ScoutMET_pt_uncorr"].array(library="ak")
                    met_gen = events["genMET_pt"].array(library="ak")
                    
                    # Take only the two leading jets per event
                    pt_corr = pt_corr[:, :2]
                    pt_uncorr = pt_uncorr[:, :2]
                    pt_gen = pt_gen[:, :2]
                    
                    # Flatten jets from all events
                    pt_corr_flat = ak.flatten(pt_corr)
                    pt_uncorr_flat = ak.flatten(pt_uncorr)
                    pt_gen_flat = ak.flatten(pt_gen)
                    
                    pt_corr_all.append(pt_corr_flat)
                    pt_uncorr_all.append(pt_uncorr_flat)
                    pt_gen_all.append(pt_gen_flat)
                    
                    # MET is per-event, not per-jet
                    met_corr_all.append(met_corr)
                    met_uncorr_all.append(met_uncorr)
                    met_gen_all.append(met_gen)
                    
                    # Store number of jets for weight calculation
                    weights_all.append(len(pt_corr_flat))
                    
                except Exception as e:
                    print(f"    Warning: Could not read Events: {e}")
                    continue
                    
        except Exception as e:
            print(f"    Error opening file: {e}")
            continue
    
    if not pt_corr_all:
        print(f"  No data read for {dataset_name}")
        return None
    
    # Concatenate all arrays
    pt_corr_concat = ak.concatenate(pt_corr_all)
    pt_uncorr_concat = ak.concatenate(pt_uncorr_all)
    pt_gen_concat = ak.concatenate(pt_gen_all)
    met_corr_concat = ak.concatenate(met_corr_all)
    met_uncorr_concat = ak.concatenate(met_uncorr_all)
    met_gen_concat = ak.concatenate(met_gen_all)
    
    # Calculate per-jet weight (cross-section weighted, no luminosity scaling)
    if total_xsec is not None and total_initial_events > 0:
        # Weight by cross-section only (no luminosity factor)
        weight_per_event = total_xsec / total_initial_events
        print(f"  Total initial events: {total_initial_events}")
        print(f"  Cross section: {total_xsec} pb")
        print(f"  Weight per event: {weight_per_event:.6e} pb")
    else:
        weight_per_event = 1.0
        print(f"  Using unit weight")
    
    # Create uniform weight arrays for each jet collection
    # Note: GenFatJet and FatJet can have different multiplicities
    weights_reco = np.full(len(pt_corr_concat), weight_per_event)
    weights_gen = np.full(len(pt_gen_concat), weight_per_event)
    weights_met = np.full(len(met_corr_concat), weight_per_event)
    
    return {
        'pt_corr': np.array(pt_corr_concat),
        'pt_uncorr': np.array(pt_uncorr_concat),
        'pt_gen': np.array(pt_gen_concat),
        'met_corr': np.array(met_corr_concat),
        'met_uncorr': np.array(met_uncorr_concat),
        'met_gen': np.array(met_gen_concat),
        'weights_reco': weights_reco,
        'weights_gen': weights_gen,
        'weights_met': weights_met,
        'n_jets_reco': len(pt_corr_concat),
        'n_jets_gen': len(pt_gen_concat),
        'n_events': len(met_corr_concat),
        'weight_sum_reco': np.sum(weights_reco),
        'weight_sum_gen': np.sum(weights_gen),
        'weight_sum_met': np.sum(weights_met)
    }


def make_comparison_plots(data_dict, output_prefix="fatjet_pt_comparison"):
    """Create comparison plots for corrected/uncorrected vs gen FatJet_pt."""
    
    # Define histogram bins
    bins = np.linspace(0, 2000, 51)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    
    # Colors
    qcd_color = '#3b5998'
    signal_color = '#d62728'
    gen_color = '#2ca02c'
    
    # Separate QCD and signal data
    qcd_data = []
    signal_data = None
    
    for name, data in data_dict.items():
        if data is None:
            continue
        is_signal = 'Signal' in name or 'mMed' in name
        if is_signal:
            signal_data = data
        else:
            qcd_data.append(data)
    
    # Plot 1: QCD
    if qcd_data:
        # Combine all QCD bins
        total_qcd_corr = None
        total_qcd_uncorr = None
        total_qcd_gen = None
        
        for data in qcd_data:
            hist_corr, _ = np.histogram(data['pt_corr'], bins=bins, weights=data['weights_reco'])
            hist_uncorr, _ = np.histogram(data['pt_uncorr'], bins=bins, weights=data['weights_reco'])
            hist_gen, _ = np.histogram(data['pt_gen'], bins=bins, weights=data['weights_gen'])
            
            if total_qcd_corr is None:
                total_qcd_corr = hist_corr
                total_qcd_uncorr = hist_uncorr
                total_qcd_gen = hist_gen
            else:
                total_qcd_corr += hist_corr
                total_qcd_uncorr += hist_uncorr
                total_qcd_gen += hist_gen
        
        # Calculate Poisson errors properly for weighted histograms
        # Error = sqrt(sum of weights^2) for each bin
        qcd_corr_err_sq = None
        qcd_uncorr_err_sq = None
        qcd_gen_err_sq = None
        
        for data in qcd_data:
            hist_corr_w2, _ = np.histogram(data['pt_corr'], bins=bins, weights=data['weights_reco']**2)
            hist_uncorr_w2, _ = np.histogram(data['pt_uncorr'], bins=bins, weights=data['weights_reco']**2)
            hist_gen_w2, _ = np.histogram(data['pt_gen'], bins=bins, weights=data['weights_gen']**2)
            
            if qcd_corr_err_sq is None:
                qcd_corr_err_sq = hist_corr_w2
                qcd_uncorr_err_sq = hist_uncorr_w2
                qcd_gen_err_sq = hist_gen_w2
            else:
                qcd_corr_err_sq += hist_corr_w2
                qcd_uncorr_err_sq += hist_uncorr_w2
                qcd_gen_err_sq += hist_gen_w2
        
        qcd_corr_err = np.sqrt(qcd_corr_err_sq)
        qcd_uncorr_err = np.sqrt(qcd_uncorr_err_sq)
        qcd_gen_err = np.sqrt(qcd_gen_err_sq)
        
        # Calculate ratios
        ratio_corr_gen = np.divide(total_qcd_corr, total_qcd_gen,
                                   where=total_qcd_gen>0, out=np.ones_like(total_qcd_corr))
        ratio_uncorr_gen = np.divide(total_qcd_uncorr, total_qcd_gen,
                                     where=total_qcd_gen>0, out=np.ones_like(total_qcd_uncorr))
        
        # Error propagation for ratios
        ratio_corr_gen_err = ratio_corr_gen * np.sqrt(
            np.divide(qcd_corr_err**2, total_qcd_corr**2, where=total_qcd_corr>0, out=np.zeros_like(total_qcd_corr)) +
            np.divide(qcd_gen_err**2, total_qcd_gen**2, where=total_qcd_gen>0, out=np.zeros_like(total_qcd_gen))
        )
        ratio_uncorr_gen_err = ratio_uncorr_gen * np.sqrt(
            np.divide(qcd_uncorr_err**2, total_qcd_uncorr**2, where=total_qcd_uncorr>0, out=np.zeros_like(total_qcd_uncorr)) +
            np.divide(qcd_gen_err**2, total_qcd_gen**2, where=total_qcd_gen>0, out=np.zeros_like(total_qcd_gen))
        )
        
        # Create figure
        fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True,
                                 gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.0})
        ax_main = axes[0]
        ax_ratio = axes[1]
        
        # Main plot
        ax_main.hist(bins[:-1], bins=bins, weights=total_qcd_corr,
                    histtype='step', linewidth=2, color=qcd_color, linestyle='-',
                    label='Corrected', alpha=0.7)
        ax_main.hist(bins[:-1], bins=bins, weights=total_qcd_uncorr,
                    histtype='step', linewidth=2, color=qcd_color, linestyle='--',
                    label='Uncorrected', alpha=0.7)
        ax_main.hist(bins[:-1], bins=bins, weights=total_qcd_gen,
                    histtype='step', linewidth=2, color=gen_color, linestyle='-',
                    label='Gen')
        
        ax_main.errorbar(bin_centers, total_qcd_corr, yerr=qcd_corr_err,
                        fmt='none', ecolor=qcd_color, elinewidth=2, alpha=0.6)
        ax_main.errorbar(bin_centers, total_qcd_uncorr, yerr=qcd_uncorr_err,
                        fmt='none', ecolor=qcd_color, elinewidth=2, alpha=0.6)
        ax_main.errorbar(bin_centers, total_qcd_gen, yerr=qcd_gen_err,
                        fmt='none', ecolor=gen_color, elinewidth=2, alpha=0.6)
        
        # Ratio plot
        ax_ratio.hist(bins[:-1], bins=bins, weights=ratio_corr_gen,
                     histtype='step', linewidth=2, color=qcd_color, linestyle='-',
                     label='Corr./Gen')
        ax_ratio.hist(bins[:-1], bins=bins, weights=ratio_uncorr_gen,
                     histtype='step', linewidth=2, color=qcd_color, linestyle='--',
                     label='Uncorr./Gen')
        
        ax_ratio.errorbar(bin_centers, ratio_corr_gen, yerr=ratio_corr_gen_err,
                         fmt='none', ecolor=qcd_color, elinewidth=2, alpha=0.6)
        ax_ratio.errorbar(bin_centers, ratio_uncorr_gen, yerr=ratio_uncorr_gen_err,
                         fmt='none', ecolor=qcd_color, elinewidth=2, alpha=0.6)
        
        # Calculate dynamic ratio limits
        valid_ratios = np.concatenate([ratio_corr_gen[total_qcd_gen > 0],
                                       ratio_uncorr_gen[total_qcd_gen > 0]])
        ratio_ylim_min = 0.3
        ratio_ylim_max = 2.0
        
        # Formatting
        ax_main.set_ylabel('Arbitrary Units', ha='right', y=1.0, fontsize=14)
        ax_main.set_yscale('log')
        
        # Add vertical line at 150 GeV
        ax_main.axvline(150, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
        
        # Add title
        ax_main.set_title(r"$\mathit{Private\ Work}$ $\mathbf{(CMS\ Simulation)}$", loc='left', fontsize=20, pad=0, x=0.0)
        
        # Create legend with QCD title
        legend = ax_main.legend(loc='upper right', frameon=True, fontsize=10, title='QCD')
        legend.get_title().set_fontsize(11)
        legend.get_title().set_fontweight('bold')
        
        ax_main.grid(True, alpha=0.2, which='both')
        
        # Remove bottom tick from top plot
        ax_main.tick_params(axis='y', which='both', labelbottom=False, labelsize=12)
        yticks = ax_main.get_yticks()
        ax_main.set_yticks(yticks[1:])  # Remove first (bottom) tick
        
        ax_ratio.set_xlabel('FatJet $p_T$ [GeV]', ha='right', x=1.0, fontsize=14)
        ax_ratio.set_ylabel('Reco / Gen', ha='right', y=1.0, fontsize=14)
        ax_ratio.set_ylim(ratio_ylim_min, ratio_ylim_max)
        ax_ratio.axhline(1.0, color='gray', linestyle='--', linewidth=1)
        ax_ratio.axvline(150, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
        ax_ratio.legend(loc='upper right', frameon=True, fontsize=9)
        ax_ratio.grid(True, alpha=0.2, which='both')
        ax_ratio.tick_params(axis='both', labelsize=12)
        
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_qcd.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_prefix}_qcd.pdf", bbox_inches='tight')
        print(f"\nSaved QCD plots to {output_prefix}_qcd.png and .pdf")
        plt.close()
    
    # Plot 2: Signal
    if signal_data is not None:
        # Create histograms
        signal_corr, _ = np.histogram(signal_data['pt_corr'], bins=bins, weights=signal_data['weights_reco'])
        signal_uncorr, _ = np.histogram(signal_data['pt_uncorr'], bins=bins, weights=signal_data['weights_reco'])
        signal_gen, _ = np.histogram(signal_data['pt_gen'], bins=bins, weights=signal_data['weights_gen'])
        
        # Calculate Poisson errors properly for weighted histograms
        # Error = sqrt(sum of weights^2) for each bin
        signal_corr_w2, _ = np.histogram(signal_data['pt_corr'], bins=bins, weights=signal_data['weights_reco']**2)
        signal_uncorr_w2, _ = np.histogram(signal_data['pt_uncorr'], bins=bins, weights=signal_data['weights_reco']**2)
        signal_gen_w2, _ = np.histogram(signal_data['pt_gen'], bins=bins, weights=signal_data['weights_gen']**2)
        
        signal_corr_err = np.sqrt(signal_corr_w2)
        signal_uncorr_err = np.sqrt(signal_uncorr_w2)
        signal_gen_err = np.sqrt(signal_gen_w2)
        
        # Calculate ratios
        ratio_corr_gen = np.divide(signal_corr, signal_gen,
                                   where=signal_gen>0, out=np.ones_like(signal_corr))
        ratio_uncorr_gen = np.divide(signal_uncorr, signal_gen,
                                     where=signal_gen>0, out=np.ones_like(signal_uncorr))
        
        # Error propagation for ratios
        ratio_corr_gen_err = ratio_corr_gen * np.sqrt(
            np.divide(signal_corr_err**2, signal_corr**2, where=signal_corr>0, out=np.zeros_like(signal_corr)) +
            np.divide(signal_gen_err**2, signal_gen**2, where=signal_gen>0, out=np.zeros_like(signal_gen))
        )
        ratio_uncorr_gen_err = ratio_uncorr_gen * np.sqrt(
            np.divide(signal_uncorr_err**2, signal_uncorr**2, where=signal_uncorr>0, out=np.zeros_like(signal_uncorr)) +
            np.divide(signal_gen_err**2, signal_gen**2, where=signal_gen>0, out=np.zeros_like(signal_gen))
        )
        
        # Create figure
        fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True,
                                 gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.0})
        ax_main = axes[0]
        ax_ratio = axes[1]
        
        # Main plot
        ax_main.hist(bins[:-1], bins=bins, weights=signal_corr,
                    histtype='step', linewidth=2, color=signal_color, linestyle='-',
                    label='Corrected', alpha=0.7)
        ax_main.hist(bins[:-1], bins=bins, weights=signal_uncorr,
                    histtype='step', linewidth=2, color=signal_color, linestyle='--',
                    label='Uncorrected', alpha=0.7)
        ax_main.hist(bins[:-1], bins=bins, weights=signal_gen,
                    histtype='step', linewidth=2, color=gen_color, linestyle='-',
                    label='Gen')
        
        ax_main.errorbar(bin_centers, signal_corr, yerr=signal_corr_err,
                        fmt='none', ecolor=signal_color, elinewidth=2, alpha=0.6)
        ax_main.errorbar(bin_centers, signal_uncorr, yerr=signal_uncorr_err,
                        fmt='none', ecolor=signal_color, elinewidth=2, alpha=0.6)
        ax_main.errorbar(bin_centers, signal_gen, yerr=signal_gen_err,
                        fmt='none', ecolor=gen_color, elinewidth=2, alpha=0.6)
        
        # Ratio plot
        ax_ratio.hist(bins[:-1], bins=bins, weights=ratio_corr_gen,
                     histtype='step', linewidth=2, color=signal_color, linestyle='-',
                     label='Corr./Gen')
        ax_ratio.hist(bins[:-1], bins=bins, weights=ratio_uncorr_gen,
                     histtype='step', linewidth=2, color=signal_color, linestyle='--',
                     label='Uncorr./Gen')
        
        ax_ratio.errorbar(bin_centers, ratio_corr_gen, yerr=ratio_corr_gen_err,
                         fmt='none', ecolor=signal_color, elinewidth=2, alpha=0.6)
        ax_ratio.errorbar(bin_centers, ratio_uncorr_gen, yerr=ratio_uncorr_gen_err,
                         fmt='none', ecolor=signal_color, elinewidth=2, alpha=0.6)
        
        # Calculate dynamic ratio limits
        valid_ratios = np.concatenate([ratio_corr_gen[signal_gen > 0],
                                       ratio_uncorr_gen[signal_gen > 0]])
        ratio_ylim_min = 0.3
        ratio_ylim_max = 2.0
        
        # Formatting
        ax_main.set_ylabel('Arbitrary Units', ha='right', y=1.0, fontsize=14)
        ax_main.set_yscale('log')
        
        # Add vertical line at 150 GeV
        ax_main.axvline(150, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
        
        # Add title
        ax_main.set_title(r"$\mathit{Private\ Work}$ $\mathbf{(CMS\ Simulation)}$", loc='left', fontsize=20, pad=0, x=0.0)
        
        # Create legend with signal info as title
        legend_title = r"SVJ $m_{Z'}$ = 1.1 TeV" + '\n' + r"$m_{\mathrm{Dark}}$ = 20 GeV, $r_{\mathrm{inv}}$ = 0.3"
        legend = ax_main.legend(loc='upper right', frameon=True, fontsize=10, title=legend_title)
        legend.get_title().set_fontsize(10)
        legend.get_title().set_fontweight('bold')
        
        ax_main.grid(True, alpha=0.2, which='both')
        
        # Remove bottom tick from top plot
        ax_main.tick_params(axis='y', which='both', labelbottom=False, labelsize=12)
        yticks = ax_main.get_yticks()
        ax_main.set_yticks(yticks[1:])  # Remove first (bottom) tick
        
        ax_ratio.set_xlabel('FatJet $p_T$ [GeV]', ha='right', x=1.0, fontsize=14)
        ax_ratio.set_ylabel('Reco / Gen', ha='right', y=1.0, fontsize=14)
        ax_ratio.set_ylim(ratio_ylim_min, ratio_ylim_max)
        ax_ratio.axhline(1.0, color='gray', linestyle='--', linewidth=1)
        ax_ratio.axvline(150, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
        ax_ratio.legend(loc='upper right', frameon=True, fontsize=9)
        ax_ratio.grid(True, alpha=0.2, which='both')
        ax_ratio.tick_params(axis='both', labelsize=12)
        
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_signal.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_prefix}_signal.pdf", bbox_inches='tight')
        print(f"Saved Signal plots to {output_prefix}_signal.png and .pdf")
        plt.close()


def make_met_comparison_plots(data_dict, output_prefix="met_comparison"):
    """Create comparison plots for corrected/uncorrected vs gen MET."""
    
    # Define histogram bins for MET
    bins = np.linspace(0, 1000, 51)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    
    # Colors
    qcd_color = '#3b5998'
    signal_color = '#d62728'
    gen_color = '#2ca02c'
    
    # Separate QCD and signal data
    qcd_data = []
    signal_data = None
    
    for name, data in data_dict.items():
        if data is None:
            continue
        is_signal = 'Signal' in name or 'mMed' in name
        if is_signal:
            signal_data = data
        else:
            qcd_data.append(data)
    
    # Plot 1: QCD MET
    if qcd_data:
        # Combine all QCD bins
        total_qcd_met_corr = None
        total_qcd_met_uncorr = None
        total_qcd_met_gen = None
        
        for data in qcd_data:
            hist_met_corr, _ = np.histogram(data['met_corr'], bins=bins, weights=data['weights_met'])
            hist_met_uncorr, _ = np.histogram(data['met_uncorr'], bins=bins, weights=data['weights_met'])
            hist_met_gen, _ = np.histogram(data['met_gen'], bins=bins, weights=data['weights_met'])
            
            if total_qcd_met_corr is None:
                total_qcd_met_corr = hist_met_corr
                total_qcd_met_uncorr = hist_met_uncorr
                total_qcd_met_gen = hist_met_gen
            else:
                total_qcd_met_corr += hist_met_corr
                total_qcd_met_uncorr += hist_met_uncorr
                total_qcd_met_gen += hist_met_gen
        
        # Calculate errors for weighted histograms
        qcd_met_corr_err_sq = None
        qcd_met_uncorr_err_sq = None
        qcd_met_gen_err_sq = None
        
        for data in qcd_data:
            hist_met_corr_w2, _ = np.histogram(data['met_corr'], bins=bins, weights=data['weights_met']**2)
            hist_met_uncorr_w2, _ = np.histogram(data['met_uncorr'], bins=bins, weights=data['weights_met']**2)
            hist_met_gen_w2, _ = np.histogram(data['met_gen'], bins=bins, weights=data['weights_met']**2)
            
            if qcd_met_corr_err_sq is None:
                qcd_met_corr_err_sq = hist_met_corr_w2
                qcd_met_uncorr_err_sq = hist_met_uncorr_w2
                qcd_met_gen_err_sq = hist_met_gen_w2
            else:
                qcd_met_corr_err_sq += hist_met_corr_w2
                qcd_met_uncorr_err_sq += hist_met_uncorr_w2
                qcd_met_gen_err_sq += hist_met_gen_w2
        
        qcd_met_corr_err = np.sqrt(qcd_met_corr_err_sq)
        qcd_met_uncorr_err = np.sqrt(qcd_met_uncorr_err_sq)
        qcd_met_gen_err = np.sqrt(qcd_met_gen_err_sq)
        
        # Calculate ratios
        ratio_met_corr_gen = np.divide(total_qcd_met_corr, total_qcd_met_gen,
                                        where=total_qcd_met_gen>0, out=np.ones_like(total_qcd_met_corr))
        ratio_met_uncorr_gen = np.divide(total_qcd_met_uncorr, total_qcd_met_gen,
                                          where=total_qcd_met_gen>0, out=np.ones_like(total_qcd_met_uncorr))
        
        # Error propagation for ratios
        ratio_met_corr_gen_err = ratio_met_corr_gen * np.sqrt(
            np.divide(qcd_met_corr_err**2, total_qcd_met_corr**2, where=total_qcd_met_corr>0, out=np.zeros_like(total_qcd_met_corr)) +
            np.divide(qcd_met_gen_err**2, total_qcd_met_gen**2, where=total_qcd_met_gen>0, out=np.zeros_like(total_qcd_met_gen))
        )
        ratio_met_uncorr_gen_err = ratio_met_uncorr_gen * np.sqrt(
            np.divide(qcd_met_uncorr_err**2, total_qcd_met_uncorr**2, where=total_qcd_met_uncorr>0, out=np.zeros_like(total_qcd_met_uncorr)) +
            np.divide(qcd_met_gen_err**2, total_qcd_met_gen**2, where=total_qcd_met_gen>0, out=np.zeros_like(total_qcd_met_gen))
        )
        
        # Create figure
        fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True,
                                 gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.0})
        ax_main = axes[0]
        ax_ratio = axes[1]
        
        # Main plot
        ax_main.hist(bins[:-1], bins=bins, weights=total_qcd_met_corr,
                    histtype='step', linewidth=2, color=qcd_color, linestyle='-',
                    label='Corrected (ScoutMET)', alpha=0.7)
        ax_main.hist(bins[:-1], bins=bins, weights=total_qcd_met_uncorr,
                    histtype='step', linewidth=2, color=qcd_color, linestyle='--',
                    label='Uncorrected (ScoutMET)', alpha=0.7)
        ax_main.hist(bins[:-1], bins=bins, weights=total_qcd_met_gen,
                    histtype='step', linewidth=2, color=gen_color, linestyle='-',
                    label='Gen MET')
        
        ax_main.errorbar(bin_centers, total_qcd_met_corr, yerr=qcd_met_corr_err,
                        fmt='none', ecolor=qcd_color, elinewidth=2, alpha=0.6)
        ax_main.errorbar(bin_centers, total_qcd_met_uncorr, yerr=qcd_met_uncorr_err,
                        fmt='none', ecolor=qcd_color, elinewidth=2, alpha=0.6)
        ax_main.errorbar(bin_centers, total_qcd_met_gen, yerr=qcd_met_gen_err,
                        fmt='none', ecolor=gen_color, elinewidth=2, alpha=0.6)
        
        # Ratio plot
        ax_ratio.hist(bins[:-1], bins=bins, weights=ratio_met_corr_gen,
                     histtype='step', linewidth=2, color=qcd_color, linestyle='-',
                     label='Corr./Gen')
        ax_ratio.hist(bins[:-1], bins=bins, weights=ratio_met_uncorr_gen,
                     histtype='step', linewidth=2, color=qcd_color, linestyle='--',
                     label='Uncorr./Gen')
        
        ax_ratio.errorbar(bin_centers, ratio_met_corr_gen, yerr=ratio_met_corr_gen_err,
                         fmt='none', ecolor=qcd_color, elinewidth=2, alpha=0.6)
        ax_ratio.errorbar(bin_centers, ratio_met_uncorr_gen, yerr=ratio_met_uncorr_gen_err,
                         fmt='none', ecolor=qcd_color, elinewidth=2, alpha=0.6)
        
        # Formatting
        ax_main.set_ylabel('Arbitrary Units', ha='right', y=1.0, fontsize=14)
        ax_main.set_yscale('log')
        
        # Add title
        ax_main.set_title(r"$\mathit{Private\ Work}$ $\mathbf{(CMS\ Simulation)}$", loc='left', fontsize=20, pad=0, x=0.0)
        
        # Create legend with QCD title
        legend = ax_main.legend(loc='upper right', frameon=True, fontsize=10, title='QCD')
        legend.get_title().set_fontsize(11)
        legend.get_title().set_fontweight('bold')
        
        ax_main.grid(True, alpha=0.2, which='both')
        
        # Remove bottom tick from top plot
        ax_main.tick_params(axis='y', which='both', labelbottom=False, labelsize=12)
        yticks = ax_main.get_yticks()
        ax_main.set_yticks(yticks[1:])
        
        ax_ratio.set_xlabel(r'$E_T^{\mathrm{miss}}$ [GeV]', ha='right', x=1.0, fontsize=14)
        ax_ratio.set_ylabel('Reco / Gen', ha='right', y=1.0, fontsize=14)
        ax_ratio.set_ylim(0.5, 1.5)
        ax_ratio.axhline(1.0, color='gray', linestyle='--', linewidth=1)
        ax_ratio.legend(loc='upper right', frameon=True, fontsize=9)
        ax_ratio.grid(True, alpha=0.2, which='both')
        ax_ratio.tick_params(axis='both', labelsize=12)
        
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_qcd.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_prefix}_qcd.pdf", bbox_inches='tight')
        print(f"\nSaved QCD MET plots to {output_prefix}_qcd.png and .pdf")
        plt.close()
    
    # Plot 2: Signal MET
    if signal_data is not None:
        # Create histograms
        signal_met_corr, _ = np.histogram(signal_data['met_corr'], bins=bins, weights=signal_data['weights_met'])
        signal_met_uncorr, _ = np.histogram(signal_data['met_uncorr'], bins=bins, weights=signal_data['weights_met'])
        signal_met_gen, _ = np.histogram(signal_data['met_gen'], bins=bins, weights=signal_data['weights_met'])
        
        # Calculate errors for weighted histograms
        signal_met_corr_w2, _ = np.histogram(signal_data['met_corr'], bins=bins, weights=signal_data['weights_met']**2)
        signal_met_uncorr_w2, _ = np.histogram(signal_data['met_uncorr'], bins=bins, weights=signal_data['weights_met']**2)
        signal_met_gen_w2, _ = np.histogram(signal_data['met_gen'], bins=bins, weights=signal_data['weights_met']**2)
        
        signal_met_corr_err = np.sqrt(signal_met_corr_w2)
        signal_met_uncorr_err = np.sqrt(signal_met_uncorr_w2)
        signal_met_gen_err = np.sqrt(signal_met_gen_w2)
        
        # Calculate ratios
        ratio_met_corr_gen = np.divide(signal_met_corr, signal_met_gen,
                                        where=signal_met_gen>0, out=np.ones_like(signal_met_corr))
        ratio_met_uncorr_gen = np.divide(signal_met_uncorr, signal_met_gen,
                                          where=signal_met_gen>0, out=np.ones_like(signal_met_uncorr))
        
        # Error propagation for ratios
        ratio_met_corr_gen_err = ratio_met_corr_gen * np.sqrt(
            np.divide(signal_met_corr_err**2, signal_met_corr**2, where=signal_met_corr>0, out=np.zeros_like(signal_met_corr)) +
            np.divide(signal_met_gen_err**2, signal_met_gen**2, where=signal_met_gen>0, out=np.zeros_like(signal_met_gen))
        )
        ratio_met_uncorr_gen_err = ratio_met_uncorr_gen * np.sqrt(
            np.divide(signal_met_uncorr_err**2, signal_met_uncorr**2, where=signal_met_uncorr>0, out=np.zeros_like(signal_met_uncorr)) +
            np.divide(signal_met_gen_err**2, signal_met_gen**2, where=signal_met_gen>0, out=np.zeros_like(signal_met_gen))
        )
        
        # Create figure
        fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True,
                                 gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.0})
        ax_main = axes[0]
        ax_ratio = axes[1]
        
        # Main plot
        ax_main.hist(bins[:-1], bins=bins, weights=signal_met_corr,
                    histtype='step', linewidth=2, color=signal_color, linestyle='-',
                    label='Corrected (ScoutMET)', alpha=0.7)
        ax_main.hist(bins[:-1], bins=bins, weights=signal_met_uncorr,
                    histtype='step', linewidth=2, color=signal_color, linestyle='--',
                    label='Uncorrected (ScoutMET)', alpha=0.7)
        ax_main.hist(bins[:-1], bins=bins, weights=signal_met_gen,
                    histtype='step', linewidth=2, color=gen_color, linestyle='-',
                    label='Gen MET')
        
        ax_main.errorbar(bin_centers, signal_met_corr, yerr=signal_met_corr_err,
                        fmt='none', ecolor=signal_color, elinewidth=2, alpha=0.6)
        ax_main.errorbar(bin_centers, signal_met_uncorr, yerr=signal_met_uncorr_err,
                        fmt='none', ecolor=signal_color, elinewidth=2, alpha=0.6)
        ax_main.errorbar(bin_centers, signal_met_gen, yerr=signal_met_gen_err,
                        fmt='none', ecolor=gen_color, elinewidth=2, alpha=0.6)
        
        # Ratio plot
        ax_ratio.hist(bins[:-1], bins=bins, weights=ratio_met_corr_gen,
                     histtype='step', linewidth=2, color=signal_color, linestyle='-',
                     label='Corr./Gen')
        ax_ratio.hist(bins[:-1], bins=bins, weights=ratio_met_uncorr_gen,
                     histtype='step', linewidth=2, color=signal_color, linestyle='--',
                     label='Uncorr./Gen')
        
        ax_ratio.errorbar(bin_centers, ratio_met_corr_gen, yerr=ratio_met_corr_gen_err,
                         fmt='none', ecolor=signal_color, elinewidth=2, alpha=0.6)
        ax_ratio.errorbar(bin_centers, ratio_met_uncorr_gen, yerr=ratio_met_uncorr_gen_err,
                         fmt='none', ecolor=signal_color, elinewidth=2, alpha=0.6)
        
        # Formatting
        ax_main.set_ylabel('Arbitrary Units', ha='right', y=1.0, fontsize=14)
        ax_main.set_yscale('log')
        
        # Add title
        ax_main.set_title(r"$\mathit{Private\ Work}$ $\mathbf{(CMS\ Simulation)}$", loc='left', fontsize=20, pad=0, x=0.0)
        
        # Create legend with signal info as title
        legend_title = r"SVJ $m_{Z'}$ = 1.1 TeV" + '\n' + r"$m_{\mathrm{Dark}}$ = 20 GeV, $r_{\mathrm{inv}}$ = 0.3"
        legend = ax_main.legend(loc='upper right', frameon=True, fontsize=10, title=legend_title)
        legend.get_title().set_fontsize(10)
        legend.get_title().set_fontweight('bold')
        
        ax_main.grid(True, alpha=0.2, which='both')
        
        # Remove bottom tick from top plot
        ax_main.tick_params(axis='y', which='both', labelbottom=False, labelsize=12)
        yticks = ax_main.get_yticks()
        ax_main.set_yticks(yticks[1:])
        
        ax_ratio.set_xlabel(r'$E_T^{\mathrm{miss}}$ [GeV]', ha='right', x=1.0, fontsize=14)
        ax_ratio.set_ylabel('Reco / Gen', ha='right', y=1.0, fontsize=14)
        ax_ratio.set_ylim(0.5, 1.5)
        ax_ratio.axhline(1.0, color='gray', linestyle='--', linewidth=1)
        ax_ratio.legend(loc='upper right', frameon=True, fontsize=9)
        ax_ratio.grid(True, alpha=0.2, which='both')
        ax_ratio.tick_params(axis='both', labelsize=12)
        
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_signal.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_prefix}_signal.pdf", bbox_inches='tight')
        print(f"Saved Signal MET plots to {output_prefix}_signal.png and .pdf")
        plt.close()


def make_response_and_resolution_plots(data_dict, output_prefix="jet_response_resolution"):
    """
    Create JEC and JER diagnostic plots.
    
    JEC diagnostic: Mean response vs pT (should be close to 1.0)
    JER diagnostic: Resolution (RMS/width) vs pT (should match between corrected and data)
    """
    
    # Define pT bins for response and resolution
    pt_bins = np.array([100, 150, 200, 300, 400, 600, 800, 1200, 2000])
    pt_centers = (pt_bins[:-1] + pt_bins[1:]) / 2
    
    # Separate QCD and signal
    qcd_data = []
    signal_data = None
    
    for name, data in data_dict.items():
        if data is None:
            continue
        is_signal = 'Signal' in name or 'mMed' in name
        if is_signal:
            signal_data = data
        else:
            qcd_data.append(data)
    
    # Process QCD
    if qcd_data:
        # Calculate response (reco/gen) in pT bins
        mean_response_corr = []
        mean_response_uncorr = []
        rms_response_corr = []
        rms_response_uncorr = []
        mean_response_corr_err = []
        mean_response_uncorr_err = []
        rms_response_corr_err = []
        rms_response_uncorr_err = []
        
        for i in range(len(pt_bins) - 1):
            pt_min, pt_max = pt_bins[i], pt_bins[i+1]
            
            # Collect all jets in this pT bin across all QCD samples
            response_corr_list = []
            response_uncorr_list = []
            weights_list = []
            
            for data in qcd_data:
                # Use gen pT for binning (standard practice)
                mask = (data['pt_gen'] >= pt_min) & (data['pt_gen'] < pt_max)
                if np.sum(mask) > 0:
                    response_corr = data['pt_corr'][mask] / data['pt_gen'][mask]
                    response_uncorr = data['pt_uncorr'][mask] / data['pt_gen'][mask]
                    weights = data['weights_reco'][mask]
                    
                    response_corr_list.extend(response_corr)
                    response_uncorr_list.extend(response_uncorr)
                    weights_list.extend(weights)
            
            if len(response_corr_list) > 0:
                response_corr_arr = np.array(response_corr_list)
                response_uncorr_arr = np.array(response_uncorr_list)
                weights_arr = np.array(weights_list)
                
                # Weighted mean
                mean_corr = np.average(response_corr_arr, weights=weights_arr)
                mean_uncorr = np.average(response_uncorr_arr, weights=weights_arr)
                
                # Weighted standard deviation (resolution)
                variance_corr = np.average((response_corr_arr - mean_corr)**2, weights=weights_arr)
                variance_uncorr = np.average((response_uncorr_arr - mean_uncorr)**2, weights=weights_arr)
                rms_corr = np.sqrt(variance_corr)
                rms_uncorr = np.sqrt(variance_uncorr)
                
                # Error on mean (uncertainty on the mean)
                n_eff = np.sum(weights_arr)**2 / np.sum(weights_arr**2)  # Effective N
                err_mean_corr = rms_corr / np.sqrt(n_eff)
                err_mean_uncorr = rms_uncorr / np.sqrt(n_eff)
                
                # Error on RMS (approximate)
                err_rms_corr = rms_corr / np.sqrt(2 * n_eff)
                err_rms_uncorr = rms_uncorr / np.sqrt(2 * n_eff)
                
                mean_response_corr.append(mean_corr)
                mean_response_uncorr.append(mean_uncorr)
                rms_response_corr.append(rms_corr)
                rms_response_uncorr.append(rms_uncorr)
                mean_response_corr_err.append(err_mean_corr)
                mean_response_uncorr_err.append(err_mean_uncorr)
                rms_response_corr_err.append(err_rms_corr)
                rms_response_uncorr_err.append(err_rms_uncorr)
            else:
                mean_response_corr.append(np.nan)
                mean_response_uncorr.append(np.nan)
                rms_response_corr.append(np.nan)
                rms_response_uncorr.append(np.nan)
                mean_response_corr_err.append(0)
                mean_response_uncorr_err.append(0)
                rms_response_corr_err.append(0)
                rms_response_uncorr_err.append(0)
        
        # Create figure with 2 subplots
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        # Plot 1: Mean Response (JEC diagnostic)
        ax1 = axes[0]
        ax1.errorbar(pt_centers, mean_response_corr, yerr=mean_response_corr_err,
                    marker='o', markersize=8, linewidth=2, capsize=4,
                    label='Corrected (JEC+JER)', color='#2ca02c')
        ax1.errorbar(pt_centers, mean_response_uncorr, yerr=mean_response_uncorr_err,
                    marker='s', markersize=8, linewidth=2, capsize=4,
                    label='Uncorrected (raw)', color='#d62728')
        
        ax1.axhline(1.0, color='gray', linestyle='--', linewidth=2, alpha=0.5, label='Perfect response')
        ax1.fill_between(pt_centers, 0.98, 1.02, alpha=0.2, color='gray', label='±2% band')
        
        ax1.set_xlabel(r'Gen Jet $p_T$ [GeV]', fontsize=14)
        ax1.set_ylabel(r'Mean Response $\langle p_T^{\mathrm{reco}} / p_T^{\mathrm{gen}} \rangle$', fontsize=14)
        ax1.set_title('JEC Diagnostic: Mean Response vs $p_T$ (QCD)', fontsize=16, pad=10)
        ax1.set_xscale('log')
        ax1.set_ylim(0.85, 1.15)
        ax1.legend(loc='best', fontsize=11, frameon=True)
        ax1.grid(True, alpha=0.3, which='both')
        
        hep.cms.label(ax=ax1, data=False, label='Simulation Preliminary', 
                     year='2018', lumi=f'{LUMI:.1f}')
        
        # Plot 2: Resolution (JER diagnostic)
        ax2 = axes[1]
        ax2.errorbar(pt_centers, rms_response_corr, yerr=rms_response_corr_err,
                    marker='o', markersize=8, linewidth=2, capsize=4,
                    label='Corrected (JEC+JER)', color='#2ca02c')
        ax2.errorbar(pt_centers, rms_response_uncorr, yerr=rms_response_uncorr_err,
                    marker='s', markersize=8, linewidth=2, capsize=4,
                    label='Uncorrected (raw)', color='#d62728')
        
        ax2.set_xlabel(r'Gen Jet $p_T$ [GeV]', fontsize=14)
        ax2.set_ylabel(r'Resolution (RMS of $p_T^{\mathrm{reco}} / p_T^{\mathrm{gen}}$)', fontsize=14)
        ax2.set_title('JER Diagnostic: Resolution vs $p_T$ (QCD)', fontsize=16, pad=10)
        ax2.set_xscale('log')
        ax2.set_ylim(0, max(max(rms_response_uncorr), max(rms_response_corr)) * 1.3)
        ax2.legend(loc='best', fontsize=11, frameon=True)
        ax2.grid(True, alpha=0.3, which='both')
        
        # Add text annotation about JER expectations
        ax2.text(0.05, 0.95, 'JER should INCREASE resolution\n(make distribution wider)',
                transform=ax2.transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        
        hep.cms.label(ax=ax2, data=False, label='Simulation Preliminary', 
                     year='2018', lumi=f'{LUMI:.1f}')
        
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_qcd.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_prefix}_qcd.pdf", bbox_inches='tight')
        print(f"Saved QCD response/resolution plots to {output_prefix}_qcd.png and .pdf")
        plt.close()
        
        # Print quantitative summary
        print("\n" + "="*80)
        print("JEC/JER Performance Summary (QCD):")
        print("="*80)
        print(f"{'pT Range [GeV]':<20} {'Mean Resp (Corr)':<18} {'Mean Resp (Uncorr)':<18} {'RMS (Corr)':<15} {'RMS (Uncorr)':<15}")
        print("-"*90)
        for i in range(len(pt_bins) - 1):
            if not np.isnan(mean_response_corr[i]):
                print(f"{pt_bins[i]:.0f}-{pt_bins[i+1]:.0f}".ljust(20),
                      f"{mean_response_corr[i]:.4f}±{mean_response_corr_err[i]:.4f}".ljust(18),
                      f"{mean_response_uncorr[i]:.4f}±{mean_response_uncorr_err[i]:.4f}".ljust(18),
                      f"{rms_response_corr[i]:.4f}".ljust(15),
                      f"{rms_response_uncorr[i]:.4f}".ljust(15))
        print("="*80)
        print("Note: RMS shows width of response distribution (JER effect)")
        print("      Corrected RMS should be LARGER than uncorr (JER widens distribution)")
        print("="*80 + "\n")


def make_response_distribution_plot(data_dict, output_prefix="jet_response_distributions"):
    """
    Plot response distributions in different pT bins to visualize resolution.
    Shows how JER increases the width of the distribution.
    """
    
    # Separate QCD data
    qcd_data = []
    for name, data in data_dict.items():
        if data is not None and 'Signal' not in name and 'mMed' not in name:
            qcd_data.append(data)
    
    if not qcd_data:
        return
    
    # Collect all data
    all_pt_gen = []
    all_response_corr = []
    all_response_uncorr = []
    all_weights = []
    
    for data in qcd_data:
        min_len = min(len(data['pt_gen']), len(data['pt_corr']))
        all_pt_gen.extend(data['pt_gen'][:min_len])
        all_response_corr.extend(data['pt_corr'][:min_len] / data['pt_gen'][:min_len])
        all_response_uncorr.extend(data['pt_uncorr'][:min_len] / data['pt_gen'][:min_len])
        all_weights.extend(data['weights_reco'][:min_len])
    
    all_pt_gen = np.array(all_pt_gen)
    all_response_corr = np.array(all_response_corr)
    all_response_uncorr = np.array(all_response_uncorr)
    
    # Create plot with 3 different pT bins to show resolution evolution
    pt_bins_to_plot = [(200, 300), (400, 600), (800, 1200)]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for ax, (pt_min, pt_max) in zip(axes, pt_bins_to_plot):
        mask = (all_pt_gen >= pt_min) & (all_pt_gen < pt_max)
        
        if np.sum(mask) > 10:
            response_corr_bin = all_response_corr[mask]
            response_uncorr_bin = all_response_uncorr[mask]
            
            # Calculate statistics
            mean_corr = np.mean(response_corr_bin)
            std_corr = np.std(response_corr_bin)
            mean_uncorr = np.mean(response_uncorr_bin)
            std_uncorr = np.std(response_uncorr_bin)
            
            # Create histogram
            bins = np.linspace(0.5, 1.5, 40)
            ax.hist(response_uncorr_bin, bins=bins, histtype='step', linewidth=2.5,
                   color='#d62728', label=f'Uncorrected\nμ={mean_uncorr:.3f}, σ={std_uncorr:.3f}',
                   density=True, alpha=0.8)
            ax.hist(response_corr_bin, bins=bins, histtype='step', linewidth=2.5,
                   color='#2ca02c', label=f'Corrected (JEC+JER)\nμ={mean_corr:.3f}, σ={std_corr:.3f}',
                   density=True, alpha=0.8)
            
            ax.axvline(1.0, color='gray', linestyle='--', linewidth=2, alpha=0.5)
            ax.set_xlabel(r'Response ($p_T^{\mathrm{reco}} / p_T^{\mathrm{gen}}$)', fontsize=12)
            ax.set_ylabel('Normalized entries', fontsize=12)
            ax.set_title(f'{pt_min}-{pt_max} GeV', fontsize=14, pad=10)
            ax.legend(loc='upper right', fontsize=10, frameon=True)
            ax.grid(True, alpha=0.3)
            ax.set_ylim(0, None)
            
            # Add annotation about JER effect
            if std_corr > std_uncorr * 1.05:
                ax.text(0.05, 0.95, f'JER increases σ\nby {(std_corr/std_uncorr-1)*100:.1f}%',
                       transform=ax.transAxes, fontsize=9, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
            elif std_corr < std_uncorr * 0.95:
                ax.text(0.05, 0.95, f'Warning: σ decreased\nby {(1-std_corr/std_uncorr)*100:.1f}%',
                       transform=ax.transAxes, fontsize=9, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='orange', alpha=0.5))
    
    fig.suptitle('Response Distributions in $p_T$ Bins (QCD)', fontsize=16, y=1.02)
    
    # Add CMS label to first plot
    hep.cms.label(ax=axes[0], data=False, label='Simulation Preliminary', 
                 year='2018', lumi=f'{LUMI:.1f}')
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}.pdf", bbox_inches='tight')
    print(f"Saved response distribution plots to {output_prefix}.png and .pdf")
    plt.close()


def main():
    """Main function to process all datasets and create plots."""
    
    print("="*80)
    print("FatJet pT Comparison: Corrected vs Uncorrected")
    print("="*80)
    
    # Read all datasets
    data_dict = {}
    for name, path in DATASETS.items():
        result = read_dataset(name, path)
        if result is not None:
            data_dict[name] = result
            print(f"  Successfully processed {name}")
            print(f"    N reco jets: {result['n_jets_reco']}")
            print(f"    N gen jets: {result['n_jets_gen']}")
            print(f"    Weighted sum (reco): {result['weight_sum_reco']:.2f}")
    
    if not data_dict:
        print("\nError: No data was successfully read!")
        return
    
    print("\n" + "="*80)
    print("Creating comparison plots...")
    print("="*80)
    
    make_comparison_plots(data_dict)
    make_met_comparison_plots(data_dict)
    make_response_and_resolution_plots(data_dict)
    make_response_distribution_plot(data_dict)
    

    print("\n" + "="*80)
    print("Done!")
    print("="*80)


if __name__ == "__main__":
    main()

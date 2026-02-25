#!/usr/bin/env python3
"""
Plot JEC correction maps from C++ output.

Reads JSON file created by create_scouting_jec_map C++ program
and creates matplotlib plots in CMS style.

Author: SVJ Scouting Analysis
Date: February 2026
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import argparse
import os

# Set CMS style
plt.style.use(hep.style.CMS)


def plot_correction_map(data, output_dir):
    """Create 2D correction map plot."""
    jet_type = data['jet_type']
    pt_bins = np.array(data['pt_bins'])
    eta_bins = np.array(data['eta_bins'])
    corrections = np.array(data['corrections'])
    counts = np.array(data['counts'])
    
    # Check if there's any data
    if counts.max() == 0:
        print(f"Skipping correction map - no data for {jet_type}")
        return
    
    # Mask bins with insufficient statistics
    corrections_masked = np.ma.masked_where(counts < 10, corrections)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    mesh = ax.pcolormesh(pt_bins, eta_bins, corrections_masked.T,
                         cmap='viridis', vmin=0.5, vmax=1.5,
                         shading='flat', edgecolors='none')
    
    ax.set_xlabel("Jet $p_T$ [GeV]", fontsize=18)
    ax.set_ylabel("Jet $\\eta$", fontsize=18)
    ax.set_xlim(pt_bins[0], pt_bins[-1])
    ax.set_ylim(eta_bins[0], eta_bins[-1])
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    cbar = plt.colorbar(mesh, ax=ax)
    cbar.set_label("Correction Factor (Offline/Scouting)", fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    
    hep.cms.label(loc=0, data=False, lumi=None, year=2017, ax=ax)
    
    jet_label = "AK4 Jets" if jet_type == "Jet" else "AK8 Jets"
    ax.text(0.95, 0.95, jet_label, transform=ax.transAxes,
            fontsize=16, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    output_file = os.path.join(output_dir, f"scouting_correction_map_{jet_type}.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def plot_uncertainty_map(data, output_dir):
    """Create 2D uncertainty map (resolution spread)."""
    jet_type = data['jet_type']
    pt_bins = np.array(data['pt_bins'])
    eta_bins = np.array(data['eta_bins'])
    uncertainties = np.array(data['uncertainties'])
    counts = np.array(data['counts'])
    
    # Check if there's any data
    if counts.max() == 0:
        print(f"Skipping uncertainty map - no data for {jet_type}")
        return
    
    # Mask bins with insufficient statistics
    uncertainties_masked = np.ma.masked_where(counts < 10, uncertainties)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    mesh = ax.pcolormesh(pt_bins, eta_bins, uncertainties_masked.T,
                         cmap='RdYlGn_r', shading='flat', edgecolors='none')
    
    ax.set_xlabel("Jet $p_T$ [GeV]", fontsize=18)
    ax.set_ylabel("Jet $\\eta$", fontsize=18)
    ax.set_xlim(pt_bins[0], pt_bins[-1])
    ax.set_ylim(eta_bins[0], eta_bins[-1])
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    cbar = plt.colorbar(mesh, ax=ax)
    cbar.set_label("Uncertainty (std dev, resolution)", fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    
    hep.cms.label(loc=0, data=False, lumi=None, year=2017, ax=ax)
    
    jet_label = "AK4 Jets" if jet_type == "Jet" else "AK8 Jets"
    ax.text(0.95, 0.95, jet_label, transform=ax.transAxes,
            fontsize=16, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    output_file = os.path.join(output_dir, f"scouting_uncertainty_map_{jet_type}.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def plot_stat_uncertainty_map(data, output_dir):
    """Create 2D statistical uncertainty map (σ/√N)."""
    jet_type = data['jet_type']
    pt_bins = np.array(data['pt_bins'])
    eta_bins = np.array(data['eta_bins'])
    stat_uncertainties = np.array(data['stat_uncertainties'])
    counts = np.array(data['counts'])
    
    # Check if there's any data
    if counts.max() == 0:
        print(f"Skipping statistical uncertainty map - no data for {jet_type}")
        return
    
    # Mask bins with insufficient statistics
    stat_uncertainties_masked = np.ma.masked_where(counts < 10, stat_uncertainties)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    mesh = ax.pcolormesh(pt_bins, eta_bins, stat_uncertainties_masked.T,
                         cmap='plasma', shading='flat', edgecolors='none')
    
    ax.set_xlabel("Jet $p_T$ [GeV]", fontsize=18)
    ax.set_ylabel("Jet $\\eta$", fontsize=18)
    ax.set_xlim(pt_bins[0], pt_bins[-1])
    ax.set_ylim(eta_bins[0], eta_bins[-1])
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    cbar = plt.colorbar(mesh, ax=ax)
    cbar.set_label("Statistical Uncertainty (σ/√N)", fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    
    hep.cms.label(loc=0, data=False, lumi=None, year=2017, ax=ax)
    
    jet_label = "AK4 Jets" if jet_type == "Jet" else "AK8 Jets"
    ax.text(0.95, 0.95, jet_label, transform=ax.transAxes,
            fontsize=16, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    output_file = os.path.join(output_dir, f"scouting_stat_uncertainty_map_{jet_type}.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def plot_statistics_map(data, output_dir):
    """Create statistics map showing number of jets per bin."""
    jet_type = data['jet_type']
    pt_bins = np.array(data['pt_bins'])
    eta_bins = np.array(data['eta_bins'])
    counts = np.array(data['counts'])
    
    # Check if there's any data
    if counts.max() == 0:
        print(f"Skipping statistics map - no data for {jet_type}")
        return
    
    counts_plot = np.ma.masked_where(counts == 0, counts)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    mesh = ax.pcolormesh(pt_bins, eta_bins, counts_plot.T,
                         cmap='YlOrRd', norm=plt.matplotlib.colors.LogNorm(vmin=1, vmax=counts.max()),
                         shading='flat', edgecolors='none')
    
    ax.set_xlabel("Jet $p_T$ [GeV]", fontsize=18)
    ax.set_ylabel("Jet $\\eta$", fontsize=18)
    ax.set_xlim(pt_bins[0], pt_bins[-1])
    ax.set_ylim(eta_bins[0], eta_bins[-1])
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    cbar = plt.colorbar(mesh, ax=ax)
    cbar.set_label("Number of Jets per Bin", fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    
    hep.cms.label(loc=0, data=False, lumi=None, year=2017, ax=ax)
    
    jet_label = "AK4 Jets" if jet_type == "Jet" else "AK8 Jets"
    ax.text(0.95, 0.95, jet_label, transform=ax.transAxes,
            fontsize=16, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    n_good_bins = np.sum(counts > 10)
    ax.text(0.05, 0.05, f"Bins with >10 jets: {n_good_bins}/{counts.size}", 
            transform=ax.transAxes,
            fontsize=14, verticalalignment='bottom', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    output_file = os.path.join(output_dir, f"scouting_statistics_map_{jet_type}.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def plot_correction_vs_pt(data, output_dir):
    """Create 1D profile: correction vs pT."""
    jet_type = data['jet_type']
    pt_bins = np.array(data['pt_bins'])
    eta_bins = np.array(data['eta_bins'])
    corrections = np.array(data['corrections'])
    uncertainties = np.array(data['uncertainties'])
    counts = np.array(data['counts'])
    
    # Check if there's any data
    if counts.max() == 0:
        print(f"Skipping pT profile - no data for {jet_type}")
        return
    
    pt_centers = 0.5 * (pt_bins[:-1] + pt_bins[1:])
    eta_centers = 0.5 * (eta_bins[:-1] + eta_bins[1:])
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot all eta bins within |eta| < 3.5
    n_eta = len(eta_centers)
    for i_eta in range(n_eta):
        # Only plot bins where both edges are within |eta| < 3.5
        if np.abs(eta_bins[i_eta]) >= 3.5 or np.abs(eta_bins[i_eta+1]) > 3.5:
            continue
        eta_label = f"{eta_bins[i_eta]:.1f} < η < {eta_bins[i_eta+1]:.1f}"
        corr_vs_pt = corrections[:, i_eta]
        unc_vs_pt = uncertainties[:, i_eta]
        mask = counts[:, i_eta] > 10
        if np.sum(mask) > 0:
            ax.errorbar(pt_centers[mask], corr_vs_pt[mask], 
                       yerr=unc_vs_pt[mask],
                       marker='o', linestyle='-', label=eta_label, alpha=0.7)
    
    ax.set_xlabel("Scouting Jet $p_T$ [GeV]", fontsize=18)
    ax.set_ylabel("Correction Factor (Offline/Scouting)", fontsize=18)
    ax.axhline(1.0, color='black', linestyle='--', alpha=0.5)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.legend(fontsize=12)
    ax.grid(alpha=0.3)
    
    hep.cms.label(loc=0, data=False, lumi=None, year=2017, ax=ax)
    
    jet_label = "AK4 Jets" if jet_type == "Jet" else "AK8 Jets"
    ax.text(0.95, 0.95, jet_label, transform=ax.transAxes,
            fontsize=16, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    output_file = os.path.join(output_dir, f"scouting_correction_vs_pt_{jet_type}.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def plot_correction_vs_eta(data, output_dir):
    """Create 1D profile: correction vs eta."""
    jet_type = data['jet_type']
    pt_bins = np.array(data['pt_bins'])
    eta_bins = np.array(data['eta_bins'])
    corrections = np.array(data['corrections'])
    uncertainties = np.array(data['uncertainties'])
    counts = np.array(data['counts'])
    
    # Check if there's any data
    if counts.max() == 0:
        print(f"Skipping eta profile - no data for {jet_type}")
        return
    
    eta_centers = 0.5 * (eta_bins[:-1] + eta_bins[1:])
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Average over pT bins
    corr_vs_eta = np.nanmean(corrections, axis=0)
    
    n_bins_with_data = np.sum(counts > 0, axis=0)
    unc_vs_eta = np.sqrt(np.nansum(uncertainties**2, axis=0))
    unc_vs_eta = np.where(n_bins_with_data > 0, unc_vs_eta / n_bins_with_data, 0)
    
    mask = np.sum(counts, axis=0) > 10
    
    ax.errorbar(eta_centers[mask], corr_vs_eta[mask], 
               yerr=unc_vs_eta[mask],
               marker='o', linestyle='-', color='tab:blue')
    
    ax.set_xlabel("Scouting Jet $\\eta$", fontsize=18)
    ax.set_ylabel("Correction Factor (Offline/Scouting)", fontsize=18)
    ax.axhline(1.0, color='black', linestyle='--', alpha=0.5)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(alpha=0.3)
    
    hep.cms.label(loc=0, data=False, lumi=None, year=2017, ax=ax)
    
    jet_label = "AK4 Jets" if jet_type == "Jet" else "AK8 Jets"
    ax.text(0.95, 0.95, jet_label, transform=ax.transAxes,
            fontsize=16, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    output_file = os.path.join(output_dir, f"scouting_correction_vs_eta_{jet_type}.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Plot JEC corrections from C++ output')
    parser.add_argument('json_file', help='JSON file from C++ program')
    parser.add_argument('--output-dir', default='test', help='Output directory for plots')
    args = parser.parse_args()
    
    print("="*80)
    print("JEC Correction Plotter")
    print("="*80)
    print(f"Input: {args.json_file}")
    print(f"Output: {args.output_dir}/")
    
    # Load data
    with open(args.json_file, 'r') as f:
        data = json.load(f)

    jet_type = data['jet_type']
    stats = data['statistics']

    print(f"\nJet type: {jet_type}")
    print(f"Statistics:")
    print(f"  Total matched: {stats['total_matched']}")
    print(f"  Rejected (negative): {stats['rejected_negative_pt']}")
    print(f"  Rejected (outliers): {stats['rejected_outliers']}")
    print(f"  Accepted: {stats['accepted']}")

    # Print mean correction and range summary
    corrections = np.array(data['corrections'])
    counts = np.array(data['counts'])
    valid = (counts >= 10) & (corrections > 0.1) & (corrections < 10.0)
    if np.any(valid):
        mean_corr = np.mean(corrections[valid])
        min_corr = np.min(corrections[valid])
        max_corr = np.max(corrections[valid])
        print(f"  Mean correction: {mean_corr:.4f} (range: {min_corr:.4f} to {max_corr:.4f})")
    else:
        print("  Mean correction: N/A (no valid bins)")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Create plots
    print(f"\nCreating plots...")
    plot_correction_map(data, args.output_dir)
    plot_uncertainty_map(data, args.output_dir)
    plot_stat_uncertainty_map(data, args.output_dir)
    plot_statistics_map(data, args.output_dir)
    plot_correction_vs_pt(data, args.output_dir)
    plot_correction_vs_eta(data, args.output_dir)
    
    print("\n" + "="*80)
    print("Done!")
    print("="*80)


if __name__ == '__main__':
    main()

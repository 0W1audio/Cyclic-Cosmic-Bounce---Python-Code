#!/usr/bin/env python3
"""
Figure 1: Bogoliubov Coefficient Magnitudes
============================================

Two-panel figure showing |α̃_k| and |β̃_k| vs k/λ in long-wavelength regime.

Author: ICB Framework V14
Date: 2025-01-XX
Paper: "Impulsive Cyclic Bounce Framework"
"""

import numpy as np
import matplotlib.pyplot as plt

# ============================================================================
# PHYSICAL CONSTANTS
# ============================================================================
ell_P = 1.616e-35      # Planck length [m]
lambda_SI = 1.2e35     # Impulse strength [m^-1]

# CMB scales
k_CMB_min = 1e-30      # Large scales (ℓ ~ 2)
k_CMB_max = 1e-25      # Small scales (ℓ ~ 3000)
k_pivot = 0.05 / 3.086e22  # 0.05 Mpc^-1 in m^-1

# ============================================================================
# BOGOLIUBOV COEFFICIENTS - UNITARY FORMS
# ============================================================================
def alpha_tilde_k(k, lambda_val=lambda_SI):
    """
    Bogoliubov coefficient α̃_k (unitary form).
    
    From Section 3.3:
    α̃_k = 1/(1 + iλ/(2k))
    
    For λ/k ≫ 1: |α̃_k| ≈ 2k/λ → 0
    
    Parameters:
    -----------
    k : float or array
        Comoving wavenumber [m^-1]
    lambda_val : float
        Impulse strength [m^-1]
        
    Returns:
    --------
    alpha_tilde : complex
        Bogoliubov coefficient (unitary form)
    """
    return 1.0 / (1.0 + 1j * lambda_val / (2.0 * k))


def beta_tilde_k(k, lambda_val=lambda_SI):
    """
    Bogoliubov coefficient β̃_k (unitary form).
    
    From Section 3.3:
    β̃_k = (-iλ/(2k))/(1 + iλ/(2k))
    
    For λ/k ≫ 1: |β̃_k| ≈ 1 (maximal excitation)
    
    Parameters:
    -----------
    k : float or array
        Comoving wavenumber [m^-1]
    lambda_val : float
        Impulse strength [m^-1]
        
    Returns:
    --------
    beta_tilde : complex
        Bogoliubov coefficient (unitary form)
    """
    x = 1j * lambda_val / (2.0 * k)
    return -x / (1.0 + x)


# ============================================================================
# FIGURE GENERATION
# ============================================================================
def generate_figure1(save_path=r'C:\Users\jean\Documents\These\Figures\fig1_bogoliubov_magnitudes.pdf'):
    """
    Generate two-panel figure: |α̃_k| and |β̃_k| vs k/λ.
    
    Parameters:
    -----------
    save_path : str
        Output filename
    """
    print("=" * 70)
    print("Figure 1: Bogoliubov Coefficient Magnitudes")
    print("=" * 70)
    
    # Dimensionless wavenumber k/λ range
    k_over_lambda_min = 1e-70  # Deep CMB
    k_over_lambda_max = 1e-55  # Still CMB
    n_points = 500
    
    k_over_lambda = np.logspace(np.log10(k_over_lambda_min), 
                                 np.log10(k_over_lambda_max), 
                                 n_points)
    
    # Physical wavenumbers
    k_values = k_over_lambda * lambda_SI
    
    print(f"\nk/λ range: {k_over_lambda_min:.2e} to {k_over_lambda_max:.2e}")
    print(f"k range: {k_values.min():.2e} to {k_values.max():.2e} m^-1")
    print(f"Number of points: {n_points}")
    print(f"λ = {lambda_SI:.2e} m^-1")
    
    # Compute coefficients
    alpha_vals = alpha_tilde_k(k_values, lambda_SI)
    beta_vals = beta_tilde_k(k_values, lambda_SI)
    
    # Magnitudes
    alpha_mag = np.abs(alpha_vals)
    beta_mag = np.abs(beta_vals)
    
    # CMB pivot in dimensionless units
    k_pivot_over_lambda = k_pivot / lambda_SI
    alpha_pivot = np.abs(alpha_tilde_k(k_pivot, lambda_SI))
    beta_pivot = np.abs(beta_tilde_k(k_pivot, lambda_SI))
    
    print(f"\nCMB pivot k* = {k_pivot:.2e} m^-1")
    print(f"k*/λ = {k_pivot_over_lambda:.2e}")
    print(f"|α̃(k*)| = {alpha_pivot:.2e}")
    print(f"|β̃(k*)| = {beta_pivot:.6f}")
    
    # ========================================================================
    # CREATE FIGURE - TWO PANELS SIDE BY SIDE
    # ========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Color palette
    color_alpha = '#0173B2'    # Blue
    color_beta = '#DE8F05'     # Orange
    color_pivot = '#CC78BC'    # Purple
    
    # ========================================================================
    # LEFT PANEL: |α̃_k| vs k/λ
    # ========================================================================
    ax1.loglog(k_over_lambda, alpha_mag, '-', color=color_alpha, linewidth=2.5,
              label=r'$|\tilde{\alpha}_k|$')
    
    # Reference line: |α̃| ∝ k/λ
    ref_line = 2 * k_over_lambda  # |α̃| ≈ 2k/λ
    ax1.loglog(k_over_lambda, ref_line, '--', color='black', linewidth=1.5,
              alpha=0.6, label=r'$\propto k/\lambda$ (asymptotic)')
    
    # Mark CMB pivot
    ax1.plot(k_pivot_over_lambda, alpha_pivot, 'o', color=color_pivot,
            markersize=10, markerfacecolor='white', markeredgewidth=2.5,
            label=r'CMB pivot $k_* = 0.05$ Mpc$^{-1}$', zorder=5)
    
    # Labels and formatting
    ax1.set_xlabel(r'Dimensionless wavenumber $k/\lambda$', fontsize=12)
    ax1.set_ylabel(r'Transmission coefficient $|\tilde{\alpha}_k|$', fontsize=12)
    ax1.set_title(r'(a) $|\tilde{\alpha}_k|$ decreases linearly', 
                 fontsize=13, fontweight='bold', pad=15)
    
    ax1.legend(fontsize=10, loc='lower right', framealpha=0.95)
    ax1.grid(True, alpha=0.3, linestyle=':', which='both')
    ax1.tick_params(labelsize=11)
    
    # Set limits
    ax1.set_xlim([k_over_lambda_min, k_over_lambda_max])
    ax1.set_ylim([1e-70, 1e-55])
    
    # ========================================================================
    # RIGHT PANEL: |β̃_k| vs k/λ
    # ========================================================================
    ax2.plot(k_over_lambda, beta_mag, '-', color=color_beta, linewidth=2.5,
            label=r'$|\tilde{\beta}_k|$')
    
    # Reference line: |β̃| = 1 (maximal)
    ax2.axhline(1.0, color='black', linestyle='--', linewidth=1.5,
               alpha=0.6, label=r'$|\tilde{\beta}_k| = 1$ (maximal)')
    
    # Mark CMB pivot
    ax2.plot(k_pivot_over_lambda, beta_pivot, 'o', color=color_pivot,
            markersize=10, markerfacecolor='white', markeredgewidth=2.5,
            label=r'CMB pivot', zorder=5)
    
    # Labels and formatting
    ax2.set_xlabel(r'Dimensionless wavenumber $k/\lambda$', fontsize=12)
    ax2.set_ylabel(r'Particle production coefficient $|\tilde{\beta}_k|$', fontsize=12)
    ax2.set_title(r'(b) $|\tilde{\beta}_k|$ saturates at unity', 
                 fontsize=13, fontweight='bold', pad=15)
    
    ax2.legend(fontsize=10, loc='lower right', framealpha=0.95)
    ax2.grid(True, alpha=0.3, linestyle=':', which='both')
    ax2.tick_params(labelsize=11)
    
    # Set limits and disable automatic offset
    ax2.set_xlim([k_over_lambda_min, k_over_lambda_max])
    ax2.set_ylim([0.9999, 1.0001])
    ax2.ticklabel_format(style='plain', axis='y', useOffset=False)  # KEY FIX!
    
    # ========================================================================
    # FINALIZE AND SAVE
    # ========================================================================
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    print("\n" + "=" * 70)
    print("FIGURE GENERATION COMPLETE")
    print("=" * 70)
    print(f"\n✓ Figure saved: {save_path}")
    print(f"  Format: PDF (vector graphics)")
    print(f"  Resolution: 300 DPI")
    print(f"  Size: {12} × {5} inches")
    print("\n✓ Long-wavelength regime (k/λ ≪ 1):")
    print(f"  |α̃_k| ∝ k/λ → 0 (transmission vanishes)")
    print(f"  |β̃_k| → 1 (maximal particle production)")
    print("=" * 70 + "\n")
    
    return fig, (ax1, ax2)


# ============================================================================
# MAIN EXECUTION
# ============================================================================
if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("ICB FRAMEWORK V14 - FIGURE 1")
    print("Bogoliubov Coefficient Magnitudes")
    print("=" * 70 + "\n")
    
    # Generate figure
    fig, axes = generate_figure1()
    
    print("\n" + "=" * 70)
    print("TO INCLUDE IN LaTeX (Section 3.3)")
    print("=" * 70)
    print(r"""
\begin{figure}[htbp]
  \centering
  \includegraphics[width=\textwidth]{Figures/fig1_bogoliubov_magnitudes.pdf}
  \caption{
    Bogoliubov coefficient magnitudes for the impulsive bounce in the 
    long-wavelength regime ($k/\lambda \ll 1$). (a) Transmission coefficient 
    $|\tilde{\alpha}_k|$ decreases linearly with $k/\lambda$ (solid blue), 
    following the asymptotic form $|\tilde{\alpha}_k| \approx 2k/\lambda$ 
    (dashed black). (b) Particle production coefficient $|\tilde{\beta}_k|$ 
    saturates at unity (maximal excitation) for all cosmological modes 
    $k \ll \lambda$. The CMB pivot scale $k_* = 0.05$ Mpc$^{-1}$ (purple 
    circles) lies deep in this regime where $k_*/\lambda \sim 10^{-61}$, 
    ensuring that the impulsive approximation is extraordinarily accurate.
  }
  \label{fig:bogoliubov_magnitudes}
\end{figure}
    """)
    print("=" * 70 + "\n")
    
    plt.show()
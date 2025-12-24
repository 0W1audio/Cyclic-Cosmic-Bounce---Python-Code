#!/usr/bin/env python3
"""
Figure 4: Particle Production n_k vs Wavenumber k
==================================================

Shows the transition from maximal excitation (n_k ~ 1) for long wavelengths
to negligible excitation for short wavelengths, with cutoff at k ~ λ.

Author: ICB Validation Suite
Date: 2025-01-XX
Paper: "Impulsive Cyclic Bounce Framework" (arXiv:XXXX.XXXXX)
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
# PARTICLE PRODUCTION FORMULA
# ============================================================================
def particle_number(k, lambda_val=lambda_SI):
    """
    Particle number per mode from Bogoliubov coefficients.
    
    From Section 3.3:
    n_k = |β_k|² = (λ/2k)²/(1 + (λ/2k)²)
    
    For λ/k ≫ 1 (CMB): n_k → 1 (maximal)
    For λ/k ≪ 1 (short): n_k → 0
    
    Parameters:
    -----------
    k : float or array
        Comoving wavenumber [m^-1]
    lambda_val : float
        Impulse strength [m^-1]
        
    Returns:
    --------
    n_k : float or array
        Particle number per mode
    """
    x = lambda_val / (2.0 * k)  # λ/(2k)
    return x**2 / (1.0 + x**2)


# ============================================================================
# FIGURE GENERATION
# ============================================================================
def generate_figure4(save_path=r'C:\Users\jean\Documents\fig4_particle_production.pdf'):
    """
    Generate n_k vs k figure showing transition at k ~ λ.
    
    Parameters:
    -----------
    save_path : str
        Output filename
    """
    print("=" * 70)
    print("Figure 4: Particle Production vs Wavenumber")
    print("=" * 70)
    
    # k range: from CMB to trans-Planckian
    k_min = 1e-32       # Beyond CMB
    k_max = 1e40        # Trans-Planckian (beyond cutoff)
    n_k = 10000
    
    k_values = np.logspace(np.log10(k_min), np.log10(k_max), n_k)
    n_k_values = particle_number(k_values, lambda_SI)
    
    print(f"\nk range: {k_min:.2e} to {k_max:.2e} m^-1")
    print(f"Number of points: {n_k}")
    print(f"λ = {lambda_SI:.2e} m^-1")
    print(f"Cutoff k_c ~ λ = {lambda_SI:.2e} m^-1")
    
    # Key scales
    print(f"\nKey scales:")
    print(f"  CMB large scales:  k ~ {k_CMB_min:.2e} m^-1, n_k = {particle_number(k_CMB_min):.6f}")
    print(f"  CMB pivot scale:   k ~ {k_pivot:.2e} m^-1, n_k = {particle_number(k_pivot):.6f}")
    print(f"  CMB small scales:  k ~ {k_CMB_max:.2e} m^-1, n_k = {particle_number(k_CMB_max):.6f}")
    print(f"  Cutoff λ:          k ~ {lambda_SI:.2e} m^-1, n_k = {particle_number(lambda_SI):.6f}")
    print(f"  Trans-Planckian:   k ~ {1e40:.2e} m^-1, n_k = {particle_number(1e40):.6f}")
    
    # ========================================================================
    # CREATE FIGURE
    # ========================================================================
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    # Color palette
    color_main = '#0173B2'      # Blue
    color_CMB = '#029E73'       # Green
    color_cutoff = '#CC78BC'    # Purple
    
    # Plot main curve
    ax.semilogx(k_values, n_k_values, '-', color=color_main, linewidth=3,
               label=r'$n_k = \frac{(\lambda/2k)^2}{1 + (\lambda/2k)^2}$',
               zorder=5)
    
    # Horizontal references
    ax.axhline(1.0, color='black', linestyle='--', linewidth=1.5, alpha=0.7,
              label=r'$n_k = 1$ (maximal)', zorder=2)
    
    ax.axhline(0.5, color='gray', linestyle=':', linewidth=1, alpha=0.5,
              label=r'$n_k = 0.5$ (half-maximal)', zorder=1)
    
    # Vertical references
    # CMB range
    ax.axvspan(k_CMB_min, k_CMB_max, alpha=0.15, color=color_CMB,
              label='CMB scales', zorder=0)
    
    # Cutoff scale
    ax.axvline(lambda_SI, color=color_cutoff, linestyle='--', linewidth=2,
              alpha=0.8, label=r'Cutoff: $k \sim \lambda$', zorder=4)
    
    # Mark CMB pivot
    ax.plot(k_pivot, particle_number(k_pivot), 'o', color='red',
           markersize=10, markerfacecolor='white', markeredgewidth=2,
           label='CMB pivot', zorder=6)
    
    # Annotations
    # CMB region
    ax.annotate('Observable\nCMB modes\n$n_k \\simeq 1$',
               xy=(np.sqrt(k_CMB_min * k_CMB_max), 0.999),
               xytext=(1e-28, 0.85),
               fontsize=11, ha='center',
               arrowprops=dict(arrowstyle='->', lw=1.5, color=color_CMB),
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Cutoff region
    ax.annotate('Cutoff:\nmodes not excited\n$n_k \\to 0$',
               xy=(lambda_SI * 10, 0.05),
               xytext=(1e37, 0.3),
               fontsize=11, ha='center',
               arrowprops=dict(arrowstyle='->', lw=1.5, color=color_cutoff),
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Transition region
    ax.annotate('Transition\n$k \\sim \\lambda$',
               xy=(lambda_SI, 0.5),
               xytext=(lambda_SI / 100, 0.5),
               fontsize=11, ha='right', va='center',
               arrowprops=dict(arrowstyle='->', lw=1.5, color=color_cutoff))
    
    # Labels and formatting
    ax.set_xlabel(r'Comoving wavenumber $k$ [m$^{-1}$]', fontsize=12)
    ax.set_ylabel(r'Particle number per mode $n_k = |\beta_k|^2$', fontsize=12)
    ax.set_title('Particle Production Across the Bounce', 
                fontsize=14, fontweight='bold', pad=15)
    
    ax.set_xlim([k_min, k_max])
    ax.set_ylim([-0.05, 1.1])
    
    ax.legend(fontsize=10, loc='upper right', framealpha=0.95, ncol=2)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.tick_params(labelsize=11)
    
    # Add text box with key insight
    textstr = (r'For $k \ll \lambda$ (all CMB modes):' + '\n'
              r'$n_k \simeq 1$ (maximal excitation)' + '\n'
              r'⟹ Superhorizon amplitude fixed' + '\n'
              r'⟹ Scale-invariant spectrum')
    props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.9, edgecolor='orange', linewidth=2)
    ax.text(0.03, 0.35, textstr, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', bbox=props)
    
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
    print(f"  Size: {10} × {6} inches")
    print("\n✓ All CMB modes maximally excited (n_k ~ 1)")
    print("✓ Natural cutoff at k ~ λ (60 orders above CMB)")
    print("=" * 70 + "\n")
    
    return fig, ax


# ============================================================================
# MAIN EXECUTION
# ============================================================================
if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("ICB VALIDATION SUITE - FIGURE 4")
    print("Particle Production vs Wavenumber")
    print("=" * 70 + "\n")
    
    # Generate figure
    fig, ax = generate_figure4(save_path=r'C:\Users\jean\Documents\fig4_particle_production.pdf')
    
    print("\n" + "=" * 70)
    print("VALIDATION COMPLETE")
    print("=" * 70)
    print("\nTo include in LaTeX document:")
    print("-" * 30)
    print(r"\begin{figure}[htbp]")
    print(r"  \centering")
    print(r"  \includegraphics[width=0.8\textwidth]{fig4_particle_production.pdf}")
    print(r"  \caption{")
    print(r"    Particle number per mode $n_k = |\beta_k|^2$ as a function of")
    print(r"    comoving wavenumber $k$. The impulsive bounce produces maximal")
    print(r"    excitation ($n_k \simeq 1$) for all long-wavelength modes")
    print(r"    $k \ll \lambda$ (green shaded: CMB scales), fixing the superhorizon")
    print(r"    amplitude independently of microphysics. The transition occurs at")
    print(r"    the cutoff scale $k \sim \lambda \sim 10^{35}$ m$^{-1}$ (purple"),
    print(r"    dashed), 60 orders of magnitude above observable scales. Modes with")
    print(r"    $k \gg \lambda$ remain unexcited ($n_k \to 0$), ensuring no")
    print(r"    trans-Planckian problem. The CMB pivot scale (red circle) lies")
    print(r"    deep in the maximal excitation regime.")
    print(r"  }")
    print(r"  \label{fig:particle_production}")
    print(r"\end{figure}")
    print("\n" + "=" * 70 + "\n")
    
    plt.show()

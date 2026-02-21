#!/usr/bin/python3
"""
Simulation of Cross-Correlation Functions (CCF) with an arbitrary number of components.
Each component is simulated by a Gaussian, convolved with an instrumental broadening kernel.
"""

import argparse
import sys
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt

def generate_ccf(rvs, vrads, ampls, vsigs):
    """
    Generate the base CCF signal from multiple Gaussian components.
    """
    ccf = np.zeros_like(rvs)
    for v, a, s in zip(vrads, ampls, vsigs):
        ccf += a * np.sqrt(2 * np.pi) * s * st.norm.pdf(rvs, v, s)
    return ccf

def add_noise(ccf, offset=0.0, amplitude=0.1):
    """
    Add random Gaussian noise to the CCF.
    """
    noise = offset + amplitude * np.random.randn(len(ccf))
    return ccf + noise

def convolve_instrument(rvs, ccf, resolution):
    """
    Convolve the CCF with a Gaussian kernel representing the instrumental broadening.
    """
    velocity_elt = 299792.458 / resolution
    broad_kernel = st.norm.pdf(rvs, 0, velocity_elt)
    
    # Normalize kernel to preserve total flux if needed, 
    # but here we follow original scaling logic
    convolved = np.convolve(ccf, broad_kernel, mode='same')
    
    # Normalize to peak amplitude of 1.0 as in original code
    if np.max(convolved) > 0:
        convolved /= np.max(convolved)
    
    return convolved, velocity_elt

def parse_args():
    parser = argparse.ArgumentParser(
        description="Simulate CCFs with multiple RV components.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-R", "--resolution", type=int, default=60000, help="Spectral resolution (resolving power)")
    parser.add_argument("-v", "--vrads", type=float, nargs='+', default=[-40.0, -10.0, 80.0], help="RV positions (km/s)")
    parser.add_argument("-a", "--ampls", type=float, nargs='+', default=[1.0, 0.7, 0.6], help="Amplitudes (rel. to peak)")
    parser.add_argument("-s", "--vsigs", type=float, nargs='+', default=[10.0, 10.0, 15.0], help="Intrinsic widths (km/s)")
    parser.add_argument("--vmin", type=float, default=-400, help="Min velocity range (km/s)")
    parser.add_argument("--vmax", type=float, default=400, help="Max velocity range (km/s)")
    parser.add_argument("--noise-offset", type=float, default=0.0, help="Noise offset")
    parser.add_argument("--noise-amp", type=float, default=0.1, help="Noise amplitude")
    parser.add_argument("-o", "--output", type=str, default=None, help="Base output filename")
    parser.add_argument("--no-plot", action="store_true", help="Disable plotting")
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Validation
    ncomp = len(args.vrads)
    if len(args.ampls) != ncomp or len(args.vsigs) != ncomp:
        print(f"Error: length of vrads ({ncomp}), ampls ({len(args.ampls)}), and vsigs ({len(args.vsigs)}) must match.")
        sys.exit(1)

    # Setup file names
    base_name = args.output if args.output else f"ccf_{args.resolution}_{ncomp}_comp"
    ofn_dat = f"{base_name}.dat"
    ofn_npy = f"{base_name}.npy"
    ofn_png = f"{base_name}.png"

    # Simulation parameters
    velocity_elt = 299792.458 / args.resolution
    dv = round(velocity_elt / 4, 1)
    nrvs = int(round((args.vmax - args.vmin) / dv))
    rvs = np.linspace(args.vmin, args.vmax, nrvs)

    print(f"--- Simulation Configuration ---")
    print(f"Resolution (R):     {args.resolution}")
    print(f"Velocity element:   {velocity_elt:8.3f} km/s")
    print(f"Sampling step (dv): {dv:8.3f} km/s ({nrvs} points)")
    print(f"Range:              [{args.vmin}, {args.vmax}] km/s")

    # Generate signal
    pure_ccf = generate_ccf(rvs, args.vrads, args.ampls, args.vsigs)
    noisy_ccf = add_noise(pure_ccf, args.noise_offset, args.noise_amp)
    final_ccf, vel_elt = convolve_instrument(rvs, noisy_ccf, args.resolution)

    # Save data
    print(f"--- Saving Outputs ---")
    print(f"Saving ASCII:       {ofn_dat}")
    data = np.array([rvs, final_ccf])
    np.savetxt(ofn_dat, data.T, fmt='%12.4e', header="Velocity(km/s)  CCF")
    
    print(f"Saving NumPy:       {ofn_npy}")
    np.save(ofn_npy, data)

    # Plotting
    if not args.no_plot:
        print(f"Saving Plot:        {ofn_png}")
        
        # Modern styling
        plt.style.use('seaborn-v0_8-muted')
        fig, ax = plt.subplots(figsize=(10, 6), tight_layout=True)
        
        ax.plot(rvs, pure_ccf / (max(pure_ccf) if max(pure_ccf) > 0 else 1), 
                '--', color='gray', alpha=0.5, label='Intrinsic Signal (Normalized)')
        ax.plot(rvs, final_ccf, color='#2c3e50', linewidth=1.5, label='Final CCF (Noisy + Broadened)')
        
        # Color palette for components
        colors = plt.cm.viridis(np.linspace(0, 0.8, ncomp))
        for i in range(ncomp):
            ax.axvline(x=args.vrads[i], color=colors[i], linestyle=':', alpha=0.8,
                      label=fr'v{i+1}: {args.vrads[i]:.1f} km/s ($\sigma$={args.vsigs[i]:.1f})')

        ax.set_xlabel('Velocity [km/s]', fontsize=12)
        ax.set_ylabel('Normalized Intensity', fontsize=12)
        ax.set_title(fr'Simulated CCF ($R$={args.resolution}, $\Delta v$={dv:.2f} km/s)', fontsize=14, pad=15)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.legend(frameon=True, facecolor='white', framealpha=0.9)
        ax.set_xlim(args.vmin, args.vmax)
        
        plt.savefig(ofn_png, dpi=150)
        plt.show()

if __name__ == "__main__":
    main()

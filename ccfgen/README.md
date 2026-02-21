# `ccfgen` — Synthetic CCF Generator

This script simulates **Cross-Correlation Functions (CCF)** with an arbitrary number of radial velocity (RV) components. Each component is modeled as a Gaussian, and the final signal includes instrumental broadening and random noise.

---

## Usage

`ccfgen.py` is a command-line tool. You can run it with default parameters or customize the simulation via arguments.

### Basic execution
```bash
python ccfgen.py
```

### Custom simulation
```bash
python ccfgen.py -R 40000 -v -25.0 10.0 90.0 -a 1.0 0.8 0.4 -s 8 8 12
```

---

## Command-Line Arguments

| Argument | Long Form | Default | Description |
| :--- | :--- | :--- | :--- |
| `-R` | `--resolution` | `60000` | Spectral resolving power ($R = \lambda / \Delta \lambda$). |
| `-v` | `--vrads` | `-40 -10 80` | List of RV positions in km/s. |
| `-a` | `--ampls` | `1.0 0.7 0.6` | List of relative amplitudes for each component. |
| `-s` | `--vsigs` | `10 10 15` | List of intrinsic Gaussian widths ($\sigma$) in km/s. |
| | `--vmin` | `-400` | Minimum velocity for the grid (km/s). |
| | `--vmax` | `400` | Maximum velocity for the grid (km/s). |
| | `--noise-amp` | `0.1` | Amplitude of the random Gaussian noise. |
| `-o` | `--output` | *auto* | Base name for output files. |
| | `--no-plot` | `False` | Run without generating/showing a plot. |

---

## Outputs

The script generates three files per run (based on the resolution and number of components, e.g., `R60k_3comp.dat`):

1.  **`.dat`**: A text file with two columns: `Velocity (km/s)` and `Normalized Intensity`.
2.  **`.npy`**: A NumPy binary file containing the same data for fast loading in Python.
3.  **`.png`**: A high-resolution diagnostic plot showing the intrinsic signal, the noisy/broadened CCF, and the component positions.

---

## Dependencies

Ensure you have the following Python libraries installed:
- `numpy`
- `scipy`
- `matplotlib`

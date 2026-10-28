"""
Spectral and radiometric calibration of pushbroom hyperspectral imagers.

Companion code to the NERC Field Spectroscopy Facility protocol
"Calibration and Characterisation of Hyperspectral Imagers".

Each detector of the unit under test (e.g. a VNIR or SWIR array) is processed
separately: edit the CONFIG at the bottom of this file and run it once per
detector. Datacubes are read as ENVI files (.hdr + raw data) with dimensions
(frames, spatial channels, spectral bands).
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import constants
from scipy.interpolate import splev, splrep
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from spectral import open_image


FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))


# %% 1. Reading datacubes

def load_cube(path):
    return np.asarray(open_image(str(path)).load(), dtype=float)


def nominal_wavelengths(path):
    centers = open_image(str(path)).bands.centers
    if not centers:
        raise ValueError(f"No wavelength information in header: {path}")
    return np.asarray(centers, dtype=float)


def frame_statistics(path):
    cube = load_cube(path)
    return cube.mean(axis=0), cube.std(axis=0, ddof=1), cube.shape[0]


def dark_corrected_frame(light_path, dark_path):
    light_mean, _, _ = frame_statistics(light_path)
    dark_mean, _, _ = frame_statistics(dark_path)
    return light_mean - dark_mean


# %% 2. Spectral calibration

def _gaussian(x, amplitude, centre, sigma):
    return amplitude * np.exp(-0.5 * ((x - centre) / sigma) ** 2)


def find_emission_lines(spectrum, flux_threshold, half_width=3):
    bands = np.arange(spectrum.size)
    peaks, _ = find_peaks(spectrum, height=flux_threshold)
    rows = []
    for peak in peaks:
        lo, hi = max(peak - half_width, 0), min(peak + half_width + 1, spectrum.size)
        if hi - lo < 4:
            continue
        try:
            (amplitude, centre, sigma), _ = curve_fit(
                _gaussian, bands[lo:hi], spectrum[lo:hi], p0=(spectrum[peak], peak, 1.0)
            )
        except RuntimeError:
            continue
        sigma = abs(sigma)
        rows.append((amplitude, centre, sigma, FWHM_PER_SIGMA * sigma))
    return pd.DataFrame(rows, columns=["Amplitude", "Centre Line", "STD", "FWHM"])


def match_emission_lines(detected, lines_nm, wavelengths, tolerance=3.0, max_fwhm=None):
    if max_fwhm is not None:
        detected = detected[detected["FWHM"] < max_fwhm]
    expected_bands = np.interp(lines_nm, wavelengths, np.arange(wavelengths.size))
    matches = []
    for line_nm, expected in zip(lines_nm, expected_bands):
        offset = (detected["Centre Line"] - expected).abs()
        if offset.empty or offset.min() > tolerance:
            print(f"  No peak found within {tolerance} bands of {line_nm} nm (expected band {expected:.1f})")
            continue
        match = detected.loc[offset.idxmin()].copy()
        match["True Wavelength"] = line_nm
        matches.append(match)
    matches = pd.DataFrame(matches)
    shared = matches.duplicated("Centre Line", keep=False)
    for line_nm in matches.loc[shared, "True Wavelength"]:
        print(f"  {line_nm} nm shares a peak with another line (unresolved); excluded")
    return matches[~shared].reset_index(drop=True)


def spectral_calibration(lamps, tolerance=3.0, polynomial_order=3, wavelengths=None):
    matched = []
    for lamp in lamps:
        frame = dark_corrected_frame(lamp["light"], lamp["dark"])
        spectrum = frame.mean(axis=0)
        reference = wavelengths if wavelengths is not None else nominal_wavelengths(lamp["light"])
        detected = find_emission_lines(spectrum, lamp["flux_threshold"])
        lines = match_emission_lines(
            detected, lamp["lines_nm"], reference, tolerance, lamp.get("max_fwhm")
        )
        lines.insert(0, "Lamp", lamp["gas"])
        matched.append(lines)

    lines = pd.concat(matched).sort_values("Centre Line").reset_index(drop=True)
    coefficients, covariance = np.polyfit(
        lines["Centre Line"], lines["True Wavelength"], polynomial_order, cov=True
    )
    lines["Fitted Wavelength"] = np.polyval(coefficients, lines["Centre Line"])
    lines["Residual (nm)"] = lines["True Wavelength"] - lines["Fitted Wavelength"]

    order = np.arange(polynomial_order, -1, -1)
    table = pd.DataFrame({
        "Term": [f"C{n}" for n in order],
        "Coefficient": coefficients,
        "Standard Error": np.sqrt(np.diag(covariance)),
    }).iloc[::-1].reset_index(drop=True)
    return coefficients, table, lines


def band_wavelengths(coefficients, n_bands):
    return np.polyval(coefficients, np.arange(n_bands))


# %% 3. Radiometric calibration

def load_sphere_radiance(csv_path, wavelength_column, radiance_column,
                         uncertainty_column, uncertainty_coverage_factor=2.0,
                         uncertainty_is_percent=True):
    sphere = pd.read_csv(csv_path)
    radiance = sphere[radiance_column].to_numpy(dtype=float)
    uncertainty = sphere[uncertainty_column].to_numpy(dtype=float) / uncertainty_coverage_factor
    if not uncertainty_is_percent:
        uncertainty = 100.0 * uncertainty / radiance
    return pd.DataFrame({
        "Wavelength": sphere[wavelength_column].to_numpy(dtype=float),
        "Radiance": radiance,
        "Relative Standard Uncertainty (%)": uncertainty,
    })


def interpolate_spectrum(x, y, x_new):
    return splev(x_new, splrep(x, y))


def radiometric_calibration(light, dark, integration_time_ms, coefficients, sphere,
                            drift_uncertainty=0.002, uniformity_uncertainty=0.004,
                            coverage_factor=2.0):
    light_mean, light_std, n_light = frame_statistics(light)
    dark_mean, dark_std, n_dark = frame_statistics(dark)

    signal = (light_mean - dark_mean) / integration_time_ms
    wavelengths = band_wavelengths(coefficients, signal.shape[1])
    radiance = interpolate_spectrum(sphere["Wavelength"], sphere["Radiance"], wavelengths)
    u_radiance = interpolate_spectrum(
        sphere["Wavelength"], sphere["Relative Standard Uncertainty (%)"], wavelengths
    ) / 100.0

    gain = radiance / signal

    u_signal = np.sqrt(light_std ** 2 / n_light + dark_std ** 2 / n_dark) / (light_mean - dark_mean)
    u_combined = np.sqrt(
        u_radiance ** 2 + drift_uncertainty ** 2 + uniformity_uncertainty ** 2 + u_signal ** 2
    )
    expanded_percent = 100.0 * coverage_factor * u_combined

    return {
        "wavelengths": wavelengths,
        "radiance": radiance,
        "signal": signal,
        "gain": gain,
        "expanded_uncertainty_percent": expanded_percent,
    }


# %% 4. Photon transfer curve, A* and NESR

def photon_transfer_curve(datasets):
    points = []
    for integration_time, paths in sorted(datasets.items()):
        light_cube = load_cube(paths["light"])
        dark_mean, _, _ = frame_statistics(paths["dark"])
        signal = (light_cube.mean(axis=0) - dark_mean).mean(axis=0)
        variance = light_cube.var(axis=0, ddof=1).mean(axis=0)
        points.append(pd.DataFrame({
            "Band": np.arange(signal.size),
            "Signal": signal,
            "Variance": variance,
            "Integration Time (ms)": integration_time,
        }))
    ptc = pd.concat(points, ignore_index=True)

    gamma, gain, read_noise = np.polyfit(ptc["Signal"], ptc["Variance"], 2)
    fitted = np.polyval((gamma, gain, read_noise), ptc["Signal"])
    residual = np.sum((ptc["Variance"] - fitted) ** 2)
    total = np.sum((ptc["Variance"] - ptc["Variance"].mean()) ** 2)
    fit = {"gamma": gamma, "G": gain, "epsilon": read_noise, "R2": 1.0 - residual / total}
    return ptc, fit


def a_star_nesr(radiometric, ptc_fit):
    wavelengths = radiometric["wavelengths"]
    photon_radiance = (radiometric["radiance"] / 1000.0) / constants.h
    signal = radiometric["signal"].mean(axis=0)
    a_star = ((1.0 / ptc_fit["G"]) * signal / photon_radiance) * 1e12
    nesr = (constants.h * constants.c / wavelengths) * (abs(ptc_fit["epsilon"]) / a_star)
    return pd.DataFrame({
        "Wavelength": wavelengths,
        "Photon Spectral Radiance": photon_radiance,
        "Normalized Signal": signal,
        "A*": a_star,
        "NESR W m-2 sr-1 nm-1": nesr,
    })


# %% 5. Plotting

FSF_COLOURS = [(115, 145, 195), (106, 162, 142), (149, 188, 99), (214, 221, 105),
               (227, 186, 105), (217, 141, 110), (192, 104, 123)]
FSF_CMAP = LinearSegmentedColormap.from_list(
    "fsf", [tuple(channel / 255.0 for channel in colour) for colour in FSF_COLOURS]
)


def set_plot_style():
    plt.rcParams.update({
        "font.family": ["Roboto", "DejaVu Sans"],
        "font.size": 14,
        "legend.fontsize": 12,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    })


def plot_radiometric_coefficients(radiometric, name, output_dir):
    fig, axes = plt.subplots(1, 2, figsize=(14, 10))
    panels = [
        (radiometric["gain"], f"{name} radiometric gain",
         "mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$ / (DN ms$^{-1}$)"),
        (radiometric["expanded_uncertainty_percent"], f"{name} expanded uncertainty (k = 2)", "%"),
    ]
    for ax, (data, title, label) in zip(axes, panels):
        image = ax.imshow(data, aspect="auto", cmap=FSF_CMAP)
        ax.set_title(title, loc="left")
        ax.set_xlabel("Spectral Band")
        ax.set_ylabel("Spatial Channel")
        fig.colorbar(image, ax=ax, shrink=0.75, label=label)
    fig.tight_layout()
    fig.savefig(output_dir / f"{name}_Radiometric_Coefficients.pdf", dpi=300, transparent=True)
    return fig


def plot_photon_transfer_curve(ptc, fit, name, output_dir):
    fig, ax = plt.subplots(figsize=(12, 8))
    times = sorted(ptc["Integration Time (ms)"].unique())
    colours = FSF_CMAP(np.linspace(0, 1, len(times)))
    for time, colour in zip(times, colours):
        subset = ptc[ptc["Integration Time (ms)"] == time]
        ax.scatter(subset["Signal"], subset["Variance"], s=40, color=colour, label=f"{time} ms")
    x = np.linspace(0, ptc["Signal"].max(), 200)
    ax.plot(x, np.polyval((fit["gamma"], fit["G"], fit["epsilon"]), x), c="#5E5E5E", linewidth=4)
    ax.annotate(
        f"$σ^2$ = {fit['gamma']:.3g}$x̄^2$ + {fit['G']:.3g}$x̄$ {fit['epsilon']:+.3g}\n"
        f"$R^2$ = {fit['R2']:.2f}",
        xy=(0.55, 0.08), xycoords="axes fraction", bbox=dict(boxstyle="square", fc="1.0"),
    )
    ax.set_title(f"{name} photon transfer curve", loc="left")
    ax.set_xlabel("Signal ($x̄$), dark subtracted, averaged across spatial channels (DN)")
    ax.set_ylabel("Temporal variance ($σ^2$, DN$^2$)")
    ax.legend(loc="upper left")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(output_dir / f"{name}_Photon_Transfer_Curve.pdf", dpi=300, transparent=True)
    return fig


def plot_a_star_nesr(characterisation, name, output_dir):
    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 6))
    left.plot(characterisation["Wavelength"], characterisation["A*"], c="#87B665", linewidth=4)
    left.set_title(f"{name} A*", loc="left")
    left.set_ylabel("A* (nm$^2$ sr)")
    right.plot(characterisation["Wavelength"], characterisation["NESR W m-2 sr-1 nm-1"],
               c="#5094A7", linewidth=4)
    right.set_title(f"{name} NESR", loc="left")
    right.set_ylabel("NESR (W m$^{-2}$ sr$^{-1}$ nm$^{-1}$)")
    for ax in (left, right):
        ax.set_xlabel("Wavelength (nm)")
        ax.grid(True)
    fig.tight_layout()
    fig.savefig(output_dir / f"{name}_A-star_NESR.pdf", dpi=300, transparent=True)
    return fig


# %% 6. Full workflow

def run(config):
    name = config["name"]
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    set_plot_style()

    print(f"{name}: spectral calibration")
    coefficients, coefficient_table, lines = spectral_calibration(**config["spectral"])
    coefficient_table.to_csv(output_dir / f"{name}_Spectral_Coefficients.csv", index=False)
    lines.to_csv(output_dir / f"{name}_Spectral_Line_Matching.csv", index=False)

    print(f"{name}: radiometric calibration")
    rad = dict(config["radiometric"])
    sphere = load_sphere_radiance(**rad.pop("sphere"))
    radiometric = radiometric_calibration(coefficients=coefficients, sphere=sphere, **rad)
    np.savetxt(output_dir / f"{name}_Radiometric_Gain.csv", radiometric["gain"], delimiter=",")
    np.savetxt(output_dir / f"{name}_Radiometric_Gain_Expanded_Uncertainty_Percent.csv",
               radiometric["expanded_uncertainty_percent"], delimiter=",")

    print(f"{name}: photon transfer curve, A* and NESR")
    ptc, ptc_fit = photon_transfer_curve(config["ptc"]["datasets"])
    ptc.to_csv(output_dir / f"{name}_Photon_Transfer_Curve.csv", index=False)
    pd.DataFrame([ptc_fit]).to_csv(output_dir / f"{name}_Photon_Transfer_Curve_Fit.csv", index=False)
    characterisation = a_star_nesr(radiometric, ptc_fit)
    characterisation.to_csv(output_dir / f"{name}_A-star_NESR.csv", index=False)

    plot_radiometric_coefficients(radiometric, name, output_dir)
    plot_photon_transfer_curve(ptc, ptc_fit, name, output_dir)
    plot_a_star_nesr(characterisation, name, output_dir)
    plt.show()

    return {
        "coefficients": coefficient_table,
        "lines": lines,
        "radiometric": radiometric,
        "ptc": ptc,
        "ptc_fit": ptc_fit,
        "characterisation": characterisation,
    }


# %% 7. Configuration

CONFIG = {
    "name": "UUT_VNIR",
    "output_dir": "output",
    "spectral": {
        "tolerance": 3.0,
        "polynomial_order": 3,
        "lamps": [
            {
                "gas": "HgAr",
                "light": "path/to/hgar_light.hdr",
                "dark": "path/to/hgar_dark.hdr",
                "flux_threshold": 100.0,
                "lines_nm": [404.66, 435.83, 546.07],
            },
            {
                "gas": "Ar",
                "light": "path/to/ar_light.hdr",
                "dark": "path/to/ar_dark.hdr",
                "flux_threshold": 80.0,
                "max_fwhm": None,
                "lines_nm": [696.54, 763.51, 826.45, 965.8],
            },
            {
                "gas": "Kr",
                "light": "path/to/kr_light.hdr",
                "dark": "path/to/kr_dark.hdr",
                "flux_threshold": 60.0,
                "max_fwhm": None,
                "lines_nm": [557.03, 587.1, 877.68, 892.87],
            },
        ],
    },
    "radiometric": {
        "light": "path/to/sphere_light.hdr",
        "dark": "path/to/sphere_dark.hdr",
        "integration_time_ms": 4.0,
        "drift_uncertainty": 0.002,
        "uniformity_uncertainty": 0.004,
        "coverage_factor": 2.0,
        "sphere": {
            "csv_path": "path/to/sphere_radiance.csv",
            "wavelength_column": "Wavelength",
            "radiance_column": "Radiance",
            "uncertainty_column": "Uncertainty",
            "uncertainty_coverage_factor": 2.0,
            "uncertainty_is_percent": True,
        },
    },
    "ptc": {
        "datasets": {
            0.5: {"light": "path/to/ptc_0p5ms_light.hdr", "dark": "path/to/ptc_0p5ms_dark.hdr"},
            1.0: {"light": "path/to/ptc_1p0ms_light.hdr", "dark": "path/to/ptc_1p0ms_dark.hdr"},
            1.5: {"light": "path/to/ptc_1p5ms_light.hdr", "dark": "path/to/ptc_1p5ms_dark.hdr"},
        },
    },
}


if __name__ == "__main__":
    results = run(CONFIG)

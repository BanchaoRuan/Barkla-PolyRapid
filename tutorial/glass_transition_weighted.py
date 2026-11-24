"""Recover breakpoint fit errors from completed PolyRapid cooling and heating runs.

Section S.8 calls the final estimate an inverse weighted arithmetic mean without
specifying the exponent. This implementation uses inverse variance, 1 / SE**2.
The SEs describe the two fits to these trajectories, in kelvin.
"""

import argparse
import csv
import io
import os
import tempfile
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit


def read_ramp(mdp):
    options = {}
    for line in mdp.read_text().splitlines():
        key, sep, value = line.partition(";")[0].partition("=")
        if sep:
            options[key.strip()] = value.strip()
    return float(options["dt"]) * int(options["nsteps"]), float(options["annealing_temp"].split()[-1])


def read_curve(path, end_time):
    data = np.loadtxt(path, comments=("#", "@"), ndmin=2)
    if data.shape[0] < 5 or data.shape[1] != 3 or not np.isfinite(data).all():
        raise ValueError(f"{path}: expected finite time, temperature and density columns")
    time, temperature, density = data.T
    if (not np.all(np.diff(time) > 0) or abs(time[0]) > 0.01
            or abs(time[-1] - end_time) > 0.01):
        raise ValueError(f"{path}: expected a complete ramp from 0 to {end_time:g} ps")
    return temperature, density


def intersecting_lines(x, x0, y0, m1, m2):
    return np.where(x < x0, m1 * (x - x0) + y0, m2 * (x - x0) + y0)


def fit_breakpoint(curve, guess, limit):
    temperature, density = curve
    domain = temperature <= limit
    x, y = temperature[domain], density[domain]
    below = x < guess
    if x.size <= 4 or min(np.unique(x[below]).size, np.unique(x[~below]).size) < 2:
        raise ValueError(f"Insufficient temperatures on each side of the {guess:g} K fit guess")
    m1, b1 = np.polyfit(x[below], y[below], 1)
    m2, b2 = np.polyfit(x[~below], y[~below], 1)
    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)
        fit, covariance = curve_fit(intersecting_lines, x, y,
                                    p0=[guess, (b1 + b2) / 2, m1, m2], maxfev=10000)
    tg, variance = fit[0], covariance[0, 0]
    if not (x.min() < tg < x.max() and np.isfinite(variance) and variance > 0):
        raise ValueError("The fit did not resolve an interior breakpoint with positive finite variance")
    return tg, variance


def weighted_result(material, guess):
    glass = material / "GlassTransition"
    analysis = glass / "Analysis"
    cooling_ps, _ = read_ramp(material / "glasstransition.mdp")
    heating_ps, limit = read_ramp(glass / "Heating" / "heating.mdp")
    cooling = read_curve(analysis / "cooling_denstemp.xvg", cooling_ps)
    heating = read_curve(analysis / "heating_denstemp.xvg", heating_ps)
    rough_tg, _ = fit_breakpoint(cooling, guess, 800)
    if limit != 50 * round((rough_tg + 100) / 50):
        raise ValueError("Cooling fit and recorded heating ramp disagree; check the material and initial guess")
    cool_tg, cool_var = fit_breakpoint(cooling, rough_tg, limit)
    heat_tg, heat_var = fit_breakpoint(heating, rough_tg, limit)
    with (analysis / "Tg_results.csv").open(newline="") as source:
        row, = csv.DictReader(source)
    original = [float(row["Cooling_Tg_linear"]), float(row["Heating_Tg_linear"])]
    if not np.allclose([cool_tg, heat_tg], original, rtol=0, atol=0.001):
        raise ValueError("Refitted breakpoints differ from Tg_results.csv by more than 0.001 K")

    # Scaling both inverse variances by the smaller variance keeps the weights bounded.
    variances = np.array([cool_var, heat_var])
    weighted = np.average([cool_tg, heat_tg], weights=variances.min() / variances)
    values = [cool_tg, np.sqrt(cool_var), heat_tg, np.sqrt(heat_var), weighted]
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(["Cooling_Tg_K", "Cooling_fit_SE_K", "Heating_Tg_K", "Heating_fit_SE_K",
                     "Inverse_variance_Tg_K"])
    writer.writerow([f"{v:.8f}" for v in values])
    output = analysis / "Tg_weighted.csv"
    with tempfile.TemporaryDirectory(dir=analysis, prefix=".tg-fit-") as staging:
        candidate = Path(staging) / output.name
        candidate.write_text(buf.getvalue())
        try:
            os.link(candidate, output)
        except FileExistsError:
            if output.read_text() != buf.getvalue():
                raise ValueError(f"{output} already exists with different results")
    print(f"Inverse-variance Tg = {weighted:.4f} K ({output})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("material", nargs="?", type=Path, default=Path.cwd())
    parser.add_argument("--GuessTg", type=float, default=200,
                        help="initial guess used by the original cooling fit (default: 200 K)")
    args = parser.parse_args()
    try:
        weighted_result(args.material, args.GuessTg)
    except (OSError, ValueError, RuntimeError, KeyError, OptimizeWarning) as exc:
        parser.exit(1, f"Tg analysis: {exc}\n")

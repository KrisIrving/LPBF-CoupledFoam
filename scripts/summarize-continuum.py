#!/usr/bin/env python3
"""Summarize M1 observations, without declaring physical validation passed."""
import argparse
import json
import math
from pathlib import Path

NUMERIC_FIELDS = (
    "time", "dt", "massKg", "outwardMassFluxKgPerS", "kineticEnergyJ",
    "absorbedLaserPowerW", "TminK", "TmaxK", "maxSpeedMS",
)


def summarize(path):
    samples = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("M1_CONTINUUM "):
            continue
        raw = dict(token.split("=", 1) for token in line.split()[1:])
        sample = {name: float(raw[name]) for name in NUMERIC_FIELDS}
        if not all(math.isfinite(value) for value in sample.values()):
            raise ValueError(f"Non-finite diagnostic value in {path.name}")
        sample["solver"] = raw["solver"]
        samples.append(sample)
    if len(samples) < 2:
        raise ValueError(f"Initial and final diagnostic samples required: {path.name}")

    residuals = []
    absorbed_energy = 0.0
    for previous, current in zip(samples, samples[1:]):
        interval = current["time"] - previous["time"]
        if interval <= 0:
            raise ValueError("Diagnostic times must increase; analyze separate restart logs separately")
        residuals.append(
            (current["massKg"] - previous["massKg"]) / interval
            + current["outwardMassFluxKgPerS"]
        )
        absorbed_energy += current["absorbedLaserPowerW"] * interval

    return {
        "solver": samples[-1]["solver"],
        "samples": len(samples),
        "initial": samples[0],
        "final": samples[-1],
        "inventory_change_kg": samples[-1]["massKg"] - samples[0]["massKg"],
        "right_endpoint_mass_residual_max_abs_kg_per_s": max(map(abs, residuals)),
        "right_endpoint_absorbed_energy_j": absorbed_energy,
        "interpretation": (
            "Observations only. The mass residual uses a right-endpoint flux estimate; "
            "it is not a general discrete residual for arbitrary time schemes, porosity, "
            "external sources or mesh changes. No full energy balance is evaluated."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("serial", type=Path)
    parser.add_argument("parallel", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    serial, parallel = summarize(args.serial), summarize(args.parallel)
    result = {
        "schema_version": 1,
        "serial": serial,
        "parallel": parallel,
        "final_inventory_difference_kg": parallel["final"]["massKg"] - serial["final"]["massKg"],
        "physical_validation": "not_evaluated",
        "field_equivalence": "not_evaluated",
    }
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"M1 observations written: {args.output}")


if __name__ == "__main__":
    main()

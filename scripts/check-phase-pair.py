#!/usr/bin/env python3
"""M1.2 isolated algebra audit of actual V3 source expressions, not a CFD test."""
import argparse
import ast
import hashlib
import json
import math
import re
from pathlib import Path


def expression(block, variable):
    # Remove comments before finding the actual assignment.
    block = re.sub(r"//[^\n]*", "", block)
    return re.search(r"\b" + variable + r"\s*=\s*([^;]+);", block).group(1)


def evaluate(expr, current, other, values):
    for token, value in {
        "alpha.thermo().rho()": current[1],
        "alpha2.thermo().rho()": other[1],
        "alpha.thermo().Cv()": current[2],
        "alpha2.thermo().Cv()": other[2],
    }.items():
        expr = expr.replace(token, repr(value))
    tree = ast.parse(expr, mode="eval")
    allowed = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Add, ast.Sub,
               ast.Mult, ast.Div, ast.USub, ast.UAdd, ast.Constant,
               ast.Name, ast.Load, ast.Call)
    for node in ast.walk(tree):
        if not isinstance(node, allowed):
            raise ValueError("Unsupported source expression")
        if isinstance(node, ast.Call) and not (
            isinstance(node.func, ast.Name) and node.func.id == "min"
            and len(node.args) == 2 and not node.keywords
        ):
            raise ValueError("Only scalar min calls are supported")
    return eval(compile(tree, "source-expression", "eval"),
                {"__builtins__": {}, "min": min},
                dict(values, alpha=current[0], alpha2=other[0]))


def audit(source):
    text = source.read_text()
    liquid = text.split('if(alpha2.name()==(alpha.name()+"vapour")){', 1)[1]
    liquid, vapour = liquid.split('else if(alpha.name()==(alpha2.name()+"vapour")){', 1)
    expressions = {name: {v: expression(block, v) for v in ("alphagen", "massgen")}
                   for name, block in (("liquid", liquid), ("vapour", vapour))}
    rows = []
    # Synthetic values are deliberate: no material/reproduction claims.
    for density_ratio in (2.0, 1000.0):
        l, v = (0.5, density_ratio, 500.0), (0.5, 1.0, 1000.0)
        for direction in ("evaporation", "condensation"):
            for cap in (1.0, 0.05):
                values = dict(evaprate=0.1, condrate=0.1, maxrate=cap,
                              evapcoefffield=float(direction == "evaporation"),
                              pair_LHG=2e6)
                sl = evaluate(expressions["liquid"]["alphagen"], l, v, values)
                sv = evaluate(expressions["vapour"]["alphagen"], v, l, values)
                heat = (evaluate(expressions["liquid"]["massgen"], l, v, values)
                        + evaluate(expressions["vapour"]["massgen"], v, l, values))
                mdot = v[1]*sv  # positive liquid -> vapour
                residual = l[1]*sl + v[1]*sv
                assert math.isclose(residual, 0, abs_tol=1e-12)
                assert math.isclose(sl + sv, mdot*(1/v[1]-1/l[1]), abs_tol=1e-12)
                assert (sl + sv > 0) == (direction == "evaporation")
                donor_rate = min(0.1, cap)*(l[0] if direction == "evaporation" else v[0])
                intended_donor_mass = donor_rate*(l[1] if direction == "evaporation" else v[1])
                rows.append(dict(direction=direction, density_ratio=density_ratio,
                                 rate_cap_per_s=cap, liquid_alpha_source_per_s=sl,
                                 vapour_alpha_source_per_s=sv,
                                 mass_residual_kg_per_m3_s=residual,
                                 volume_source_per_s=sl+sv,
                                 transferred_mass_kg_per_m3_s=mdot,
                                 magnitude_over_donor_rate_mass=abs(mdot)/intended_donor_mass,
                                 upstream_temperature_source_kg_K_per_m3_s=heat,
                                 reference_latent_sink_W_per_m3=mdot*values['pair_LHG']))
    return dict(schema_version=1, source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                extracted_expressions=expressions, cases=rows,
                scope="Fixed densities/Cv, identical pair rates; original raw sources only. "
                      "No advection, PCR correction, EOS, pressure or temperature integration. "
                      "Temperature source has different units from latent power; do not compare directly.",
                physical_validation="not_evaluated")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source = root / 'applications/solvers/compressibleLaserbeamFoam/multiphaseMixtureThermo/multiphaseMixtureThermo.C'
    result = audit(source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print('M1.2 raw source pairing: 8 algebra checks passed; CFD validation pending.')
    print('Evaporation transfer / donor-rate mass = 1 / density ratio; interpretation needs review.')
    print(f'Report: {args.output}')


if __name__ == '__main__':
    main()

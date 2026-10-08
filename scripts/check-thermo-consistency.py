#!/usr/bin/env python3
"""Audit pressure-energy identity for the documented synthetic Test1 thermos.

Analytic formulas checked against OpenCFD v2512 source, not a live thermo call.
This is a necessary local thermodynamic identity, not full EOS validation.
"""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path


def number(text, key):
    return float(re.search(r'\b'+key+r'\s+([\deE.+-]+)\s*;', text).group(1))


def audit(root):
    directory = root/'tutorials/compressiblelaserbeamFoam/Test1/constant'
    liquid_path = directory/'thermophysicalProperties.metal1'
    vapour_path = directory/'thermophysicalProperties.metal1vapour'
    liquid, vapour = liquid_path.read_text(), vapour_path.read_text()
    if not re.search(r'thermo\s+hConst\s*;', liquid) or not re.search(r'equationOfState\s+adiabaticPerfectFluid\s*;', liquid):
        raise ValueError('Liquid analytic audit requires hConst + adiabaticPerfectFluid')
    if not re.search(r'thermo\s+eConst\s*;', vapour) or not re.search(r'equationOfState\s+perfectFluid\s*;', vapour):
        raise ValueError('Vapour analytic audit requires eConst + perfectFluid')
    rho0, p0, gamma, bulk = (number(liquid, k) for k in ('rho0', 'p0', 'gamma', 'B'))
    gas_rho0, gas_R = (number(vapour, k) for k in ('rho0', 'R'))
    temperature = 4101.0
    cases = []
    for pressure in (80000.0, 100000.0, 120000.0):
        rho = rho0*((pressure+bulk)/(p0+bulk))**(1/gamma)
        psi = rho/(gamma*(pressure+bulk))
        # hConst Hs = Cp*(T-Tref)+Hsref and EOS H=0; Es=Hs-p/rho.
        native_derivative = -1/rho + pressure*psi/rho**2
        required_derivative = pressure*psi/rho**2  # rho_T = 0
        gap = native_derivative-required_derivative
        assert math.isclose(gap, -1/rho, rel_tol=1e-12)
        # Verify the differentiated native pressure term numerically.
        def native_pressure_term(p):
            density = rho0*((p+bulk)/(p0+bulk))**(1/gamma)
            return -p/density
        finite_difference = (native_pressure_term(pressure+1)-native_pressure_term(pressure-1))/2
        assert math.isclose(finite_difference, native_derivative, rel_tol=1e-8)
        gas_rho = gas_rho0+pressure/(gas_R*temperature)
        gas_psi = 1/(gas_R*temperature)
        gas_rho_T = -pressure/(gas_R*temperature**2)
        gas_required = (temperature*gas_rho_T+pressure*gas_psi)/gas_rho**2
        assert abs(gas_required) < 1e-12
        cases.append(dict(pressure_Pa=pressure, temperature_K=temperature,
            liquid_density_kg_per_m3=rho,
            liquid_native_de_dp_J_per_kg_Pa=native_derivative,
            liquid_required_de_dp_J_per_kg_Pa=required_derivative,
            liquid_identity_gap_J_per_kg_Pa=gap,
            vapour_native_de_dp_J_per_kg_Pa=0.0,
            vapour_required_de_dp_J_per_kg_Pa=gas_required,
            liquid_pressure_identity_passed=False,
            vapour_pressure_identity_passed=True))
    return dict(schema_version=1,
        source_dictionary_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in (liquid_path, vapour_path)},
        identity='(de/dp)_T = [T*(d rho/dT)_p + p*(d rho/dp)_T]/rho^2',
        analytic_formula_version='Checked against OpenCFD v2512 EOS and caloric source',
        cases=cases, physical_validation='not_evaluated',
        interpretation='The liquid native caloric/EOS pairing fails this necessary pressure-energy '
                       'identity by -1/rho. A constant phase reference offset cannot change that derivative. '
                       'This does not establish the entire solver energy error or validate the vapour '
                       'model beyond this one identity. No production properties or equations were changed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit(Path(__file__).resolve().parents[1])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('Analytic checks completed: liquid pressure-energy identity fails; vapour passes this identity only.')
    print(f'Report: {args.output}')

"""Material data and magnetic co-energy densities.

B-H curves are tabulated as (B [T], H [A/m]). They are extended beyond the last point
with a vacuum slope so the Newton iteration stays well defined in deep saturation.
"""
import math
from ngsolve import BSpline

mu0 = 4e-7*math.pi

BH = {
    # Co-Fe (Hiperco 50 class), annealed
    'hiperco50': ([0, 0.5, 1.0, 1.5, 1.8, 2.0, 2.1, 2.2, 2.25, 2.3, 2.35, 2.45, 2.6],
                  [0, 40, 80, 150, 250, 500, 800, 1500, 2500, 5000, 10000, 40000, 150000]),
    # low-carbon steel (1010 class)
    'steel1010': ([0, 0.5, 1.0, 1.3, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.3],
                  [0, 150, 300, 600, 1100, 1800, 3000, 5500, 10000, 18000, 35000, 120000]),
    # quenched and tempered Cr-Mo steel (40Cr / 4140 / 4340 class), estimated
    'qt40cr':    ([0, 0.5, 1.0, 1.3, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.3],
                  [0, 700, 1500, 3000, 5500, 8000, 12000, 19000, 32000, 55000, 100000, 240000]),
    # non-oriented electrical steel M330-35A
    'm330_35a':  ([0, 0.5, 1.0, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1],
                  [0, 55, 95, 130, 170, 250, 500, 1300, 3600, 8000, 16000, 30000, 60000]),
}


def add_bh(name, B, H):
    """Register a user B-H curve (lists of B in T and H in A/m, both starting at 0)."""
    BH[name] = (list(B), list(H))


def energy_cf(name, b, sf=1.0, hscale=1.0):
    """Magnetic energy density w(|B|) = int_0^|B| H(b') db' as an NGSolve CoefficientFunction.

    sf      lamination stacking factor: B_eff = sf*B_mat(H) + (1-sf)*mu0*H
    hscale  multiplies H at every B (e.g. 1.5 to model cutting damage)
    """
    B, H = BH[name]
    B = list(B); H = [h*hscale for h in H]
    if sf != 1.0:
        B = [sf*bb + (1-sf)*mu0*hh for bb, hh in zip(B, H)]
    for _ in range(7):                       # extend with the vacuum slope
        H.append(H[-1] + 0.5/mu0); B.append(B[-1] + 0.5)
    sp = BSpline(2, [0] + B, H)              # piecewise-linear H(B)
    return sp.Integrate()(b)

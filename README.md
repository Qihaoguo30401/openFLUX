# openFLUX

openFLUX is a 3D nonlinear magnetostatic finite-element solver for axial-flux permanent-magnet (AFPM) machines, built on [NGSolve/Netgen](https://ngsolve.org).

It builds a parametric 3D model of one periodic sector of the machine, solves the magnetic field with nonlinear iron, and returns torque and coil flux linkage for any rotor position and coil currents.

## Method

- **Formulation:** magnetic vector potential **A** on Nédélec (H(curl)) edge elements, solved over a periodic sector.
- **Nonlinear iron:** the iron enters through its magnetic energy density. The energy functional is minimised with Newton's method and a backtracking line search, and each linear system is solved with BDDC-preconditioned conjugate gradients.
- **Torque:** the Arkkio air-gap volume integral.
- **Geometry:**
  - Q slots and NP poles, with concentrated tooth coils.
  - Rounded tooth corners.
  - Optional Halbach (tangentially magnetised) pieces between the poles.
  - Glue joints between magnet pieces.
  - An optional rotor back-iron ring.
  - An optional symmetry plane at the magnet mid-plane, for double-stator machines with an ironless rotor.
- **Materials:** Hiperco 50, 1010 steel, quenched-and-tempered 40Cr/4140 and M330-35A are included. You can add your own B-H curves with `openflux.add_bh`. Laminated cores take a stacking factor.

Validation is in [`examples/benchmark_energies2020`](examples/benchmark_energies2020). It models a published 45-slot / 10-pole double-stator AFPM prototype (Energies 2020, 13(8), 2122). Its torque constant is within 5.6–6.7 % of the measurement at 100–175 A, and above it at lower current. The comparison and the full method are in the accompanying report, [`paper/main.pdf`](paper/main.pdf), with LaTeX source in `paper/`.

## Install

```bash
pip install -e .            # or: pip install ngsolve numpy
pip install matplotlib      # for the example post-processing
```

Python 3.10 or newer is required.

## Quick start

```python
import math
from openflux import AxialFluxModel

p = dict(Q=12, NP=8, nsec=4, Ri=40e-3, Ro=70e-3, ws=6e-3, hs=15e-3, ys=6e-3, g=1.0e-3,
         hm=4e-3, alpha=0.8, yr=5e-3, bi_ri=38e-3, bi_ro=72e-3, sym_top=False)
m = AxialFluxModel(theta_m=0.0, pgeo=p, Br=1.25, steel='m330_35a',
                   mesh_kw=dict(h_gap=0.8e-3, h_mag=2e-3, h_iron=3e-3, h_coil=3e-3))
NI = [300.0, -150.0, -150.0] # ampere-turns of each tooth coil in the sector (Q/nsec = 3 coils)
res = m.solve(NI)
print(res['T'], res['lam'])  # torque [N m], flux linkage per turn of each coil [Wb]
```

The sector must contain a whole number of slots and an even number of poles, because the boundary condition is periodic rather than anti-periodic. The full parameter list is documented at the top of [`openflux/geometry.py`](openflux/geometry.py).

To reproduce the benchmark:

```bash
cd examples/benchmark_energies2020
python run_benchmark.py phase 0.0                        # current angle of maximum torque
python run_benchmark.py 0,68,139,209 0.0 -1.7196         # load points at rotor position 0 deg
python analyse.py                                        # compare with the measurements
```

Each 3D solve uses roughly 4–6 GB of memory and takes 5–15 minutes on one core, depending on the mesh. Run solves one at a time on a 16 GB machine.

## License

openFLUX is released under the [Apache License 2.0](LICENSE).

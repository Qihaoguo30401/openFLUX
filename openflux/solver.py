"""3D nonlinear magnetostatic solver for axial-flux PM machines (NGSolve).

Magnetic vector potential A on Nedelec (HCurl) elements over a periodic sector.
The nonlinear iron enters through its magnetic energy density, and the total energy
functional is minimised with Newton's method and a backtracking line search; the
linearised systems are solved with BDDC-preconditioned CG. Torque is computed with
the Arkkio air-gap volume integral.

Coil excitations are given directly as ampere-turns per tooth coil.
"""
import math, time
from ngsolve import *
from . import geometry
from .materials import energy_cf, mu0


class AxialFluxModel:
    """Mesh one rotor position and solve it for any set of coil currents.

    theta_m   rotor mechanical angle [rad]
    pgeo      geometry parameters, overriding geometry.P
    mesh_kw   local mesh sizes: h_gap, h_mag, h_iron, h_coil [m]
    Br, mur   magnet remanence [T] and recoil permeability
    steel     stator B-H curve name in materials.BH, sf its stacking factor
    ring_steel  B-H curve of the rotor back-iron ring (if yr > 0)
    order     HCurl element order (1 is adequate for torque with a fine gap mesh)
    """
    def __init__(s, theta_m, pgeo=None, mesh_kw=None, Br=1.2, mur=1.05, steel='hiperco50', sf=1.0,
                 ring_steel='steel1010', order=1):
        geo, s.info = geometry.build(theta_m, p=pgeo, **(mesh_kw or {}))
        s.mesh = Mesh(geo.GenerateMesh(maxh=10e-3)); s.p = s.info['p']
        s.Br, s.mur, s.steel, s.sf, s.ring_steel = Br, mur, steel, sf, ring_steel
        s.fes = Periodic(HCurl(s.mesh, order=order, dirichlet='outer', nograds=True))
        print(f'mesh {s.mesh.ne} tets, {s.fes.ndof} dofs', flush=True)
        s.S = s.info['tc']*s.info['hs']
        s.ncoil = s.p['Q']//s.p['nsec']
        s.dirs = {}
        for name, (k, d) in s.info['coilmat'].items():
            if d[0] == 'bar': s.dirs[name] = (k, CF((d[1][0], d[1][1], 0)))
            else:
                B = d[1]; rx = x-B[0]; ry = y-B[1]; rr = sqrt(rx*rx+ry*ry)+1e-9
                s.dirs[name] = (k, CF((-ry/rr, rx/rr, 0)))
        s.gfu = GridFunction(s.fes)
        h1 = Periodic(H1(s.mesh, order=1, dirichlet='outer')); s.h1 = h1
        ph_, q_ = h1.TnT(); s.ah1 = BilinearForm(grad(ph_)*grad(q_)*dx).Assemble()
        s.ah1inv = s.ah1.mat.Inverse(h1.FreeDofs(), inverse='sparsecholesky')

    def coil_dir_cf(s, w):
        return s.mesh.MaterialCF({n: w[k]*cf for n, (k, cf) in s.dirs.items()}, default=CF((0,0,0)))

    def solve(s, NI, tol=1e-4, maxit=40):
        """NI: ampere-turns of each tooth coil in the sector (length Q/nsec). Returns dict(T, lam)."""
        mesh, fes = s.mesh, s.fes; u = fes.TrialFunction()
        J = s.coil_dir_cf([v/s.S for v in NI])
        if max(abs(v) for v in NI) > 0:
            q_ = s.h1.TestFunction(); f1 = LinearForm(J*grad(q_)*dx('coil.*')).Assemble()
            phi = GridFunction(s.h1); phi.vec.data = s.ah1inv*f1.vec; J = J - grad(phi)
        nu0 = 1/mu0; num = 1/(mu0*s.mur); B = curl(u); b = sqrt(B*B + 1e-12)
        magnames = '|'.join(s.info['magmat'].keys())
        a = BilinearForm(fes, symmetric=True)
        a += Variation(energy_cf(s.steel, b, s.sf)*dx('stator'))
        if 'ring' in mesh.GetMaterials(): a += Variation(energy_cf(s.ring_steel, b)*dx('ring'))
        a += Variation((0.5*nu0*B*B)*dx('air|gap|coil.*'))
        Mv = mesh.MaterialCF({n: CF(d) for n, d in s.info['magdir'].items()}, default=CF((0,0,0)))
        a += Variation((0.5*num*B*B - num*s.Br*(Mv*B))*dx(magnames))
        a += Variation((-J*u)*dx)
        a += Variation(1e-6*nu0*0.5*u*u*dx)
        pre = Preconditioner(a, 'bddc'); gfu = s.gfu
        res = gfu.vec.CreateVector(); du = gfu.vec.CreateVector(); t = time.time()
        for it in range(maxit):
            a.AssembleLinearization(gfu.vec); a.Apply(gfu.vec, res); pre.Update()
            inv = CGSolver(a.mat, pre.mat, tol=1e-8, maxiter=2000, printrates=False)
            du.data = inv*res; e0 = a.Energy(gfu.vec); w = 1.0
            for _ in range(8):
                gfu.vec.data -= w*du
                if a.Energy(gfu.vec) <= e0 + 1e-12*abs(e0): break
                gfu.vec.data += w*du; w *= 0.5
            err = math.sqrt(abs(InnerProduct(du, res)))
            if it == 0: err0 = err
            if err < tol*err0 or err < 1e-12: break
        print(f'  newton {it} its {time.time()-t:.0f}s', flush=True)
        return s.post()

    def post(s):
        mesh = s.mesh; Bv = curl(s.gfu); g = s.p['g']; ns = s.p['nsec']
        r = sqrt(x*x+y*y)+1e-12; Bt = (-y*Bv[0]+x*Bv[1])/r; Bz = Bv[2]
        # Arkkio torque over the whole machine [N m] (one stator side if sym_top)
        T = ns*Integrate(r*Bt*Bz/mu0*dx('gap'), mesh)/g
        # flux linkage per turn of each tooth coil in the sector [Wb]
        lam = [Integrate(s.coil_dir_cf([1.0 if j == k else 0 for j in range(s.ncoil)])*s.gfu*dx('coil.*'), mesh)/s.S
               for k in range(s.ncoil)]
        return dict(T=T, lam=lam)


# backwards-compatible name
ModelAF = AxialFluxModel

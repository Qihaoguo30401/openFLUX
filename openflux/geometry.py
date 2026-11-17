"""Parametric 3D geometry of a single-stator / single-rotor axial-flux PM machine sector (Netgen OCC).

Q slots / NP poles with concentrated tooth coils, built over one of `nsec` periodic sectors.
Options: Halbach (tangential) pieces between poles, a rotor back-iron ring, a symmetry plane
at the magnet mid-plane (double-stator machines with an ironless rotor), rounded tooth corners
and glue joints between magnet pieces.

Units: metres. Sector edges lie on slot centre lines so no tooth or coil is cut.
Q/nsec and NP/nsec must be integers and the sector must be electrically periodic
(an even number of poles per sector), since the boundary condition is periodic, not anti-periodic.
"""
import math, numpy as np
from netgen.occ import *

# Default parameters. Every key can be overridden through `p` in build().
#   Q, NP     slots, poles            nsec    periodic sectors in the full machine
#   Ri, Ro    active inner/outer radius of teeth and magnets
#   ws, hs    slot (tooth-coil) width, tooth height      ys  stator yoke thickness
#   g         air gap                 hm      magnet thickness (axial)
#   alpha     pole-arc ratio of the axially magnetised pieces
#   hal       None, or +1/-1: add tangentially magnetised Halbach pieces between the poles
#   yr        rotor back-iron ring thickness (0 = none), between radii bi_ri and bi_ro
#   sym_top   True: symmetry plane at the magnet mid-plane (model half of a double-stator machine)
#   glue      glue joint width between magnet pieces      rc_out, rc_in  tooth corner radii
#   Rbox      (r_in, r_out) of the surrounding air region
P = dict(Q=45, NP=10, nsec=5, Ri=70e-3, Ro=122.5e-3, ws=7e-3, hs=20e-3, ys=8.5e-3, g=1.125e-3,
         hm=5e-3, yr=0.0, alpha=2/3, bi_ri=68e-3, bi_ro=124.5e-3, sym_top=True)

def build(theta_m=0.0, p=None, h_gap=0.35e-3, h_mag=0.7e-3, h_iron=1.2e-3, h_coil=1.0e-3):
    """theta_m: rotor mechanical angle (rad). Returns (geo, info)."""
    p = dict(P, **(p or {}))
    Q, NP, Ri, Ro, ws, hs, ys, g, hm, yr = [p[k] for k in ('Q','NP','Ri','Ro','ws','hs','ys','g','hm','yr')]
    tc = ws/2                       # coil side fills half the slot (liner merged into the coil)
    nsec = p.get('nsec', 3); nt = Q//nsec; npl = NP//nsec
    a0 = -math.pi/Q; span = 2*math.pi/nsec
    z_t1 = ys+hs; z_m0 = z_t1+g; z_m1 = z_m0+hm; z_r1 = z_m1+yr
    Rin, Rout = p.get('Rbox', (Ri - 12e-3, Ro + 13e-3)); zlo = -6e-3
    zhi = z_m1 if p.get('sym_top') else z_r1+6e-3

    def wedge(r0, r1, z0, z1, th0=a0, th1=a0+span):
        c = Cylinder((0,0,z0), Z, r=r1, h=z1-z0) - Cylinder((0,0,z0), Z, r=r0, h=z1-z0)
        # sector cut with two half-spaces through the z axis
        big = 1.0
        hs0 = Box((-big,0,-big),(big,big,big)).Rotate(Axis((0,0,0),Z), math.degrees(th0))
        hs1 = Box((-big,-big,-big),(big,0,big)).Rotate(Axis((0,0,0),Z), math.degrees(th1))
        return c * hs0 * hs1

    # ---- stator: yoke + teeth
    yoke_ri = p.get('yoke_ri', Ri - tc - 0.5e-3); yoke_ro = p.get('yoke_ro', Ro + tc + 0.5e-3)
    stator = wedge(yoke_ri, yoke_ro, 0, ys)
    teeth = []; coils = []
    for k in range(nt):
        th = 2*math.pi*k/Q; ds = math.pi/Q
        def edge(sign, r):
            a = th + sign*ds
            n = np.array([-math.sin(a), math.cos(a)])*(-sign)
            u = np.array([math.cos(a), math.sin(a)])
            t = math.sqrt(r*r - (ws/2)**2)
            return t*u + n*ws/2
        pts = [edge(-1,Ri), edge(-1,Ro), edge(+1,Ro), edge(+1,Ri)]
        # make CCW
        ar = sum(pts[i][0]*pts[(i+1)%4][1]-pts[(i+1)%4][0]*pts[i][1] for i in range(4))
        if ar < 0: pts = pts[::-1]
        # optional rounded tooth corners: radius rc_out at the Ro corners, rc_in at the Ri corners
        rc = [p.get('rc_out', 0.0) if np.linalg.norm(q) > (Ri+Ro)/2 else p.get('rc_in', 0.0) for q in pts]
        tan_in, tan_out, cen, turn = [], [], [], []
        for i in range(4):
            Pv, A_, C_ = pts[i], pts[i-1], pts[(i+1)%4]
            ta = (Pv-A_)/np.linalg.norm(Pv-A_); tb = (C_-Pv)/np.linalg.norm(C_-Pv)
            ext = math.atan2(ta[0]*tb[1]-ta[1]*tb[0], ta@tb)          # exterior (left) turn angle
            d = rc[i]*math.tan(ext/2); nin = np.array([-ta[1], ta[0]])  # left normal = inward for CCW
            tan_in.append(Pv-ta*d); tan_out.append(Pv+tb*d); cen.append(Pv-ta*d+nin*rc[i]); turn.append(math.degrees(ext))
        wp = WorkPlane(Axes((0,0,ys), n=Z, h=X))
        t0 = (pts[1]-pts[0])/np.linalg.norm(pts[1]-pts[0])
        wp.MoveTo(*tan_out[0]).Direction(*t0)
        for i in range(1,5):
            j = i % 4
            wp.LineTo(*tan_in[j])
            if rc[j] > 0: wp.Arc(rc[j], turn[j])
            else: wp.Rotate(turn[j])
        tooth = wp.Face().Extrude(hs)
        teeth.append(tooth)
        # coil pieces: 4 straight bars + 4 corner sectors, each with a known current direction
        for i in range(4):
            j = (i+1) % 4
            A, B = tan_out[i], tan_in[j]
            t = (B-A)/np.linalg.norm(B-A); n = np.array([t[1], -t[0]])      # outward normal (CCW)
            wp = WorkPlane(Axes((0,0,ys), n=Z, h=X))
            q = [A, A+tc*n, B+tc*n, B]
            wp.MoveTo(*q[0]);  [wp.LineTo(*x) for x in q[1:]]
            bar = wp.Close().Face().Extrude(hs)
            bar.mat(f'coil{k}_b{i}'); coils.append((bar, k, ('bar', t)))
            ang = turn[j]; R = rc[j]
            wp = WorkPlane(Axes((0,0,ys), n=Z, h=X))
            if R > 0:
                wp.MoveTo(*B).Direction(*n).Line(tc).Rotate(90).Arc(R+tc, ang).Rotate(90).Line(tc).Rotate(90).Arc(R, -ang)
                cor = wp.Face().Extrude(hs); Bc = cen[j]
            else:
                wp.MoveTo(*B).Direction(*n).Line(tc).Rotate(90).Arc(tc, ang)
                cor = wp.LineTo(*B).Face().Extrude(hs); Bc = B
            cor.mat(f'coil{k}_c{i}'); coils.append((cor, k, ('corner', Bc)))
    stator = stator + sum(teeth[1:], teeth[0])
    stator.mat('stator')

    # ---- magnets (rotor rotated by theta_m)
    pole = 2*math.pi/NP; half = p['alpha']*pole/2; dg = p.get('glue', 0.0)/2/((Ri+Ro)/2)
    mags = []; mdir = {}
    for k in range(-2, npl+2):
        c = k*pole + pole/2 + theta_m
        m = wedge(Ri, Ro, z_m0, z_m1, c-half+dg, c+half-dg)
        # clip to sector
        m = m * wedge(Rin*0.9, Rout*1.1, zlo-1, zhi+1)
        if m.solids and m.mass > 1e-12:
            sgn = 1 if (k % 2 == 0) else -1
            m.mat(f'mag{"N" if sgn>0 else "S"}{k}'); mags.append((m, sgn))
            mdir[m.solids[0].name] = (0.0, 0.0, float(sgn))
    if p.get('hal'):   # tangential (parallel-magnetised) pieces between poles; boundary k-1|k has M_theta = hal*sgn_k
        for k in range(-2, npl+3):
            b = k*pole + theta_m; sgk = 1 if (k % 2 == 0) else -1
            m = wedge(Ri, Ro, z_m0, z_m1, b-(pole/2-half)+dg, b+(pole/2-half)-dg)
            m = m * wedge(Rin*0.9, Rout*1.1, zlo-1, zhi+1)
            if m.solids and m.mass > 1e-12:
                s_ = p['hal']*sgk
                m.mat(f'magT{k}'); mags.append((m, 0))
                mdir[m.solids[0].name] = (-s_*math.sin(b), s_*math.cos(b), 0.0)
    gap = wedge(Ri-4e-3, Ro+4e-3, z_t1, z_m0); gap.mat('gap')
    air = wedge(Rin, Rout, zlo, zhi); air.mat('air')
    solids = [stator] + [c[0] for c in coils] + [m[0] for m in mags] + [gap]
    if p.get('yr', 0) > 0 and not p.get('sym_top'):
        ring = wedge(p['bi_ri'], p['bi_ro'], z_m1, z_r1); ring.mat('ring'); solids.append(ring)
    # air = box minus everything else
    rest = air - sum(solids[1:], solids[0])
    rest.mat('air')
    shape = Glue(solids + [rest])
    # mesh sizes
    for so in shape.solids:
        nm = so.name
        if nm == 'gap': so.maxh = h_gap
        elif nm.startswith('mag') or nm == 'ring': so.maxh = h_mag
        elif nm == 'stator': so.maxh = h_iron
        elif nm.startswith('coil'): so.maxh = h_coil
        else: so.maxh = 5e-3
    # periodic faces: plane at a0 and a0+span
    def onplane(f, th):
        c = f.center; return abs(-math.sin(th)*c.x + math.cos(th)*c.y) < 1e-7 and (math.cos(th)*c.x+math.sin(th)*c.y) > 0
    fa = [f for f in shape.faces if onplane(f, a0)]
    fb = [f for f in shape.faces if onplane(f, a0+span)]
    rot = Rotation(Axis((0,0,0), Z), math.degrees(span))
    npair = 0
    for f in fa:
        c = f.center; cr = np.array([c.x*math.cos(span)-c.y*math.sin(span), c.x*math.sin(span)+c.y*math.cos(span), c.z])
        best = min(fb, key=lambda q: np.linalg.norm(np.array([q.center.x,q.center.y,q.center.z])-cr))
        d = np.linalg.norm(np.array([best.center.x,best.center.y,best.center.z])-cr)
        assert d < 1e-7, ('unmatched periodic face', d)
        f.Identify(best, 'per', IdentificationType.PERIODIC, rot); npair += 1
    for f in shape.faces:
        c = f.center; r = math.hypot(c.x, c.y)
        if abs(r-Rin) < 1e-6 or abs(r-Rout) < 1e-6 or abs(c.z-zlo) < 1e-9: f.bc('outer')
        elif abs(c.z-zhi) < 1e-9: f.bc('sym' if p.get('sym_top') else 'outer')
    for f in fa: f.bc('perA')
    for f in fb: f.bc('perB')
    info = dict(mags=mags,
                z=dict(t1=z_t1, m0=z_m0, m1=z_m1, r1=z_r1), tc=tc, hs=hs, npair=npair, p=p)
    # coil info keyed by material name
    info['coilmat'] = {}
    for (s, k, d) in coils:
        info['coilmat'][s.solids[0].name] = (k, d)
    info['magmat'] = {m[0].solids[0].name: m[1] for m in mags}
    info['magdir'] = mdir
    return OCCGeometry(shape), info

if __name__ == '__main__':
    import time
    t = time.time(); geo, info = build(0.0)
    print('geo', time.time()-t, 'periodic pairs', info['npair'])
    t = time.time(); m = geo.GenerateMesh(maxh=4e-3)
    print('mesh', time.time()-t, m.ne, 'elements')

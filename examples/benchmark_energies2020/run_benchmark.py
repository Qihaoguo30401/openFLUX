# Benchmark A: Energies 2020, 13(8), 2122 (DSSR axial-flux PMSM, 45 slots / 10 poles, strip-wound M330-35A,
# ironless rotor with 10 mm N33EH magnets, 1.125 mm gaps, OD/ID 245/140, slots 7 x 20 mm, yoke 8.5 mm).
# Half model (one stator + half the magnet thickness, symmetry plane at the rotor mid-plane), 1/5 sector.
# Usage (from this folder):
#   python run_benchmark.py phase <pos_deg_mech>                       find the max-torque current angle g0
#   python run_benchmark.py <NI list, e.g. 0,68,139> <pos_deg_mech> <g0> [Br]   load points
# Results are appended to bench2122.jsonl; analyse.py compares them with the published measurements.
import sys, os, math, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from openflux import AxialFluxModel
pos = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
Br = float(sys.argv[4]) if len(sys.argv) > 4 else 1.15
p = dict(Q=45, NP=10, nsec=5, Ri=70e-3, Ro=122.5e-3, ws=7e-3, hs=20e-3, ys=8.5e-3, g=1.125e-3,
         hm=5e-3, yr=0.0, alpha=120/180, sym_top=True, Rbox=(58e-3, 135e-3))
pp = 5; Qs = 9
def coil_mmf(X, ge):
    """X = peak ampere-turns of one coil branch; ge = electrical angle of the current vector. Returns tooth-coil NI."""
    sect = {0:('A',1),1:('C',-1),2:('B',1),3:('A',-1),4:('C',1),5:('B',-1)}
    ph0 = {'A':0, 'B':-2*math.pi/3, 'C':2*math.pi/3}
    i = []
    for c in range(Qs):
        a = (c*40) % 360; ph, sg = sect[int(((a+30) % 360)//60)]
        i.append(sg*X*math.cos(ge + ph0[ph]))
    S = [i[j] - i[(j-4) % Qs] for j in range(Qs)]          # outward slot current, slot j between tooth j and j+1
    C = [0.0]
    for k in range(Qs-1): C.append(C[-1] + S[k])           # C_{k+1} - C_k = S_k
    m = sum(C)/Qs; return [c-m for c in C]
t = time.time()
m = AxialFluxModel(math.radians(pos), p, dict(h_gap=0.9e-3, h_mag=3e-3, h_iron=4e-3, h_coil=5e-3), Br=Br, steel='m330_35a', sf=0.95)
print('setup', time.time()-t, flush=True)
mode = sys.argv[1]
out = open('bench2122.jsonl', 'a')
if mode == 'phase':        # find the current angle of maximum torque at low current
    X = 100.0; T = []
    for ge in (0.0, math.pi/2):
        T.append(m.solve(coil_mmf(X, ge + pp*math.radians(pos)))['T'])
    g0 = math.atan2(T[1], T[0]); res = dict(mode='phase', pos=pos, T=T, g0=g0, Tmax=math.hypot(*T))
    print(json.dumps(res), flush=True); out.write(json.dumps(res)+'\n')
else:
    g0 = float(sys.argv[3]); Xs = [float(v) for v in mode.split(',')]
    for X in Xs:
        r = m.solve(coil_mmf(X, g0 + pp*math.radians(pos)))
        res = dict(mode='load', pos=pos, X=X, T=r['T'], Br=Br, sec=time.time()-t)
        print(json.dumps(res), flush=True); out.write(json.dumps(res)+'\n'); out.flush()

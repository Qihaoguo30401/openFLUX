# Post-process benchmark A (Energies 2020, 13(8), 2122): torque constant vs current, compared with the
# published measurements and the authors' own 3D FEA. Reads bench2122.jsonl, writes benchmark_results.json and benchmarkA_kT.png.
import json, math, sys, collections
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d = [json.loads(l) for l in open('bench2122.jsonl')]
ph = [r for r in d if r['mode'] == 'phase'][0]
cog = {r['pos']: r['T'] for r in d if r['mode'] == 'load' and r['X'] == 0}   # cogging (no-load) torque per position
load = [r for r in d if r['mode'] == 'load' and r['X'] > 0]
Ke_LL = 0.62                       # measured line-line voltage constant, V·s/rad (RMS, mechanical)
kT_emf = math.sqrt(3)*Ke_LL        # torque per A_rms implied by the measured back-EMF (sinusoidal, Id = 0)
k_lin = 2*math.hypot(ph['T'][0]-cog[0.0], ph['T'][1]-cog[0.0])/100           # N·m per peak coil-branch ampere-turn, both stators, low current
n = kT_emf/(k_lin*math.sqrt(2)/2)  # coil turns / parallel paths implied by the measured back-EMF
print(f'kT from back-EMF {kT_emf:.3f} N·m/A, implied turns/paths n = {n:.2f}')
nn = round(n)
I_meas = [24, 49, 74, 100, 124, 150, 175]
kT_meas = [0.926, 0.964, 0.977, 0.975, 0.969, 0.960, 0.947]
kT_pFEA = [1.053, 1.044, 1.037, 0.977, 0.951, 0.932, 0.905]
byX = collections.defaultdict(list)
for r in load: byX[r['X']].append(2*(r['T'] - cog[r['pos']]))   # both stators, ripple-corrected
rows = []
for X, I, km, kp in zip(sorted(byX), I_meas, kT_meas, kT_pFEA):
    Ia = X/(nn*math.sqrt(2)/2); T = sum(byX[X])/len(byX[X]); k = T/Ia
    rows.append(dict(I=I, X=X, npos=len(byX[X]), T_ours=T, kT_ours=k, kT_meas=km, kT_paperFEA=kp,
                     err_ours=100*(k/km-1), err_paper=100*(kp/km-1)))
    print(f'{I:4d} A  ours {k:.3f}  meas {km:.3f}  paper-FEA {kp:.3f}  ours-vs-meas {100*(k/km-1):+5.1f} %  paper {100*(kp/km-1):+5.1f} %  ({len(byX[X])} pos)')
json.dump(dict(Ke_LL=Ke_LL, kT_emf=kT_emf, n_implied=n, rows=rows), open('benchmark_results.json', 'w'), indent=1)
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(I_meas, kT_meas, 'o-', color='#222', label='Measured (Energies 2020, 13, 2122)')
ax.plot(I_meas, kT_pFEA, 's--', color='#c0392b', label="Authors' 3D FEA")
ax.plot(I_meas, [r['kT_ours'] for r in rows], 'd-', color='#2471a3', label='openFLUX 3D (Ke-anchored)')
ax.axhline(kT_emf, color='#888', lw=0.8, ls=':'); ax.text(178, kT_emf-0.003, '√3·Ke (measured back-EMF)', ha='right', va='top', fontsize=8, color='#666')
ax.set_xlabel('Phase current (A rms)'); ax.set_ylabel('Torque constant (N·m/A rms)')
ax.set_title('Benchmark A: 45-slot/10-pole DSSR axial-flux motor'); ax.grid(alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig('benchmarkA_kT.png', dpi=150)

"""Step 4: quantify what the plots suggest, and test the anomalies instead of assuming.
 A. I/Q noise anomaly across captures
 B. Velocity step: size, time, ratio
 C. Onset: rise-time fit
 D. Second track search
 E. Does SNR fall like a receding point target (R^-n)? Independent check of the 'away' sign
"""
import numpy as np
from scipy import signal, optimize
from analyze import load, C
from zoom import ridge

EV = {1: (0.2, 1.5), 2: (0.45, 1.56), 3: (1.7, 2.58), 4: (1.2, 1.96)}
D = {i: load(i) for i in range(1, 5)}
IQ = {i: D[i][0] for i in D}
FS = 50000
LAM = C / 24e9


def sec(t):
    print("\n" + "=" * 8, t)


# ---------------------------------------------------------------- A
sec("A. Noise-only region (first 10000 samples = 0.2 s, before any event)")
N = 10000
for i in IQ:
    x = IQ[i][:N]
    print(f" cap{i}: std I={x.real.std():.2f} Q={x.imag.std():.2f}  kurt I={((x.real-x.real.mean())**4).mean()/x.real.var()**2:.2f} "
          f"Q={((x.imag-x.imag.mean())**4).mean()/x.imag.var()**2:.2f}")
print(" regress I_k on I_1 (slope, residual std):")
for k in (2, 3, 4):
    a, b = IQ[1][:N].real, IQ[k][:N].real
    s = (a @ b) / (a @ a)
    print(f"   cap{k}: slope={s:.3f} resid std={np.std(b - s*a):.2f}   (I std {b.std():.2f})")
print(" same for Q:")
for k in (2, 3, 4):
    a, b = IQ[1][:N].imag, IQ[k][:N].imag
    s = (a @ b) / (a @ a)
    print(f"   cap{k}: slope={s:.3f} resid std={np.std(b - s*a):.2f}   (Q std {b.std():.2f})")
a, b = IQ[1][:N].real, IQ[2][:N].real
xc = signal.correlate(b - b.mean(), a - a.mean(), mode="full") / (N * a.std() * b.std())
lag = xc.argmax() - (N - 1)
print(f" cap1 vs cap2 I cross-corr peak {xc.max():.3f} at lag {lag}")
print(" windowed I-corr cap1 vs cap2 (2000-sample windows) over common length:")
L = min(IQ[1].size, IQ[2].size)
row = []
for s0 in range(0, L - 2000, 2000):
    c = np.corrcoef(IQ[1].real[s0:s0 + 2000], IQ[2].real[s0:s0 + 2000])[0, 1]
    row.append(f"{s0/FS:.2f}s:{c:+.2f}")
print("  ", " ".join(row))
print(" Gaussianity: first-10000 I values quantised to int8; unique values / range:",
      np.unique(IQ[1][:N].real).size)
# integer structure of the I sequence: are cap1 and cap2 samples integer-equal often?
print(" fraction of samples where I1==I2 exactly:", np.mean(IQ[1][:N].real == IQ[2][:N].real).round(3))
print(" spectral flatness of noise-only I vs Q (cap1, 8 x averaged 1024-pt Welch; band std in dB):")
for name, v in (("I", IQ[1][:N].real), ("Q", IQ[1][:N].imag)):
    f, p = signal.welch(v, fs=FS, nperseg=1024)
    bands = np.array_split(10 * np.log10(p), 8)
    print(f"   {name}: " + " ".join(f"{b.mean():.1f}" for b in bands))

# ---------------------------------------------------------------- B, C
sec("B/C. Tracks, step, rise")
TR = {}
for i in IQ:
    iq = IQ[i]
    r = ridge(iq, FS, LAM, 2048, 256, *EV[i], min_snr=10)
    # keep main connected track: drop isolated jumps >15 m/s
    first = int(np.argmax(r[:, 2] > 20))   # seed the track at the first strong frame, not at noise
    r = r[first:]
    keep = [0]
    for k in range(1, len(r)):
        if abs(r[k, 1] - r[keep[-1], 1]) < 8:
            keep.append(k)
    r = r[keep]
    TR[i] = r
    dv = np.diff(r[:, 1]); dt = np.diff(r[:, 0])
    j = np.argmin(dv)  # most negative jump
    pre, post = r[j, 1], r[j + 1, 1]
    # robust pre/post using a few points either side
    pre_m = np.median(r[max(0, j - 3):j + 1, 1]); post_m = np.median(r[j + 1:j + 5, 1])
    print(f" cap{i}: onset~{r[0,0]:.3f}s ... end {r[-1,0]:.3f}s | step at t={r[j,0]:.3f}-{r[j+1,0]:.3f}s "
          f"v {pre_m:.2f} -> {post_m:.2f} m/s  drop={pre_m-post_m:.2f} m/s  ratio post/pre={post_m/pre_m:.3f}  "
          f"(jump spans {dt[j]*1e3:.0f} ms between ridge frames)")
    # decay slope before step (after plateau, last 0.3 s before step) and after
    m1 = (r[:, 0] > r[j, 0] - 0.30) & (r[:, 0] <= r[j, 0])
    m2 = (r[:, 0] > r[j + 1, 0]) & (r[:, 0] < r[j + 1, 0] + 0.30)
    if m1.sum() > 4 and m2.sum() > 4:
        s1 = np.polyfit(r[m1, 0], r[m1, 1], 1)[0]; s2 = np.polyfit(r[m2, 0], r[m2, 1], 1)[0]
        print(f"        slope before step {s1:+.2f} m/s^2 (v~{pre_m:.1f}), after step {s2:+.2f} m/s^2 (v~{post_m:.1f})")
    # peak speed
    print(f"        peak v_r={r[:,1].max():.2f} m/s = {r[:,1].max()*3.6:.1f} km/h at t={r[r[:,1].argmax(),0]:.3f}; "
          f"onset v_r={r[1,1]:.2f}")
    # rise fit on onset..peak+: v = vinf - A exp(-(t-t0)/tau) on first 0.25 s
    t0 = r[0, 0]
    m = r[:, 0] < t0 + 0.30
    # use frames with decent snr
    m &= r[:, 2] > 25
    tt = r[m, 0] - t0; vv = r[m, 1]
    f = lambda t, vinf, A, tau, s: vinf + s * t - A * np.exp(-t / tau)
    try:
        p, _ = optimize.curve_fit(f, tt, vv, p0=[vv.max(), 3, 0.05, -0.5], maxfev=20000)
        print(f"        rise fit: v_inf={p[0]:.2f}, amplitude={p[1]:.2f} m/s, tau={p[2]*1e3:.0f} ms, slope={p[3]:+.2f} m/s^2")
    except Exception as e:
        print("        rise fit failed", e)
    # integrated radial distance onset->step
    tstep = r[j, 0]
    mm = r[:, 0] <= tstep
    dist = np.trapezoid(r[mm, 1], r[mm, 0])
    print(f"        radial path from onset to step = {dist:.1f} m over {tstep-t0:.3f}s")

# ---------------------------------------------------------------- D
sec("D. Second-track search (mask +/-1.5 m/s around main ridge, look for others > 12 dB over noise)")
for i in IQ:
    f, t, Z = signal.stft(IQ[i], fs=FS, window="hann", nperseg=1024, noverlap=1024 - 128,
                          return_onesided=False, boundary=None, padded=False)
    o = np.argsort(f); f = f[o]; P = np.abs(Z[o]) ** 2; v = f * LAM / 2
    floor = np.median(P)
    r = TR[i]
    cnt = 0; tracks = []
    for k in range(P.shape[1]):
        if not (r[0, 0] - 0.01 <= t[k] <= r[-1, 0]):
            continue
        vm = np.interp(t[k], r[:, 0], r[:, 1])
        col = P[:, k].copy()
        col[(abs(v - vm) < 1.5) | (v < 0) | (v > 70)] = 0
        j = col.argmax()
        snr = 10 * np.log10(col[j] / floor)
        if snr > 12:
            cnt += 1; tracks.append((t[k], v[j], v[j] - vm, snr))
    print(f" cap{i}: frames with 2nd peak: {cnt}")
    for tr_ in tracks[::max(1, len(tracks) // 8)]:
        print(f"     t={tr_[0]:.3f} v2={tr_[1]:.2f} (main {tr_[1]-tr_[2]:.2f}, delta {tr_[2]:+.2f})  snr={tr_[3]:.1f}")

# ---------------------------------------------------------------- E
sec("E. SNR(t) vs receding point target: SNR = K - 10 n log10(R0 + integral v dt)")
for i in IQ:
    r = TR[i]
    t0 = r[0, 0]
    m = r[:, 0] > t0 + 0.06   # skip onset transient
    tt, vv, ss = r[m, 0], r[m, 1], r[m, 2]
    # exclude step discontinuity effect: keep all, integral uses measured v
    cum = np.concatenate([[0], np.cumsum(0.5 * (vv[1:] + vv[:-1]) * np.diff(tt))])

    def res(p, n=None):
        K, R0 = p[0], p[1]
        nn = 4.0 if n is None else p[2]
        R = max(R0, 1e-3) + cum
        return ss - (K - 10 * nn * np.log10(R))

    p4 = optimize.least_squares(res, [60, 3.0], bounds=([-100, 0.1], [300, 200]))
    pn = optimize.least_squares(lambda p: res(p, True), [60, 3.0, 4.0], bounds=([-100, 0.1, 0.5], [300, 200, 12]))
    print(f" cap{i}: SNR fall {ss[0]:.1f}->{ss[-1]:.1f} dB over radial path {cum[-1]:.1f} m | "
          f"n=4: R0={p4.x[1]:.2f} m rms={np.sqrt(np.mean(p4.fun**2)):.2f} dB | "
          f"free n={pn.x[2]:.2f}, R0={pn.x[1]:.2f} m rms={np.sqrt(np.mean(pn.fun**2)):.2f} dB")

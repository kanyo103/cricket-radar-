"""Step 2: ridge-track the Doppler line, find onset, speed, deceleration, step events, secondary tracks."""
import numpy as np
from scipy import signal
from analyze import load, C
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NPERSEG, HOP = 512, 64  # 10.2 ms window, 1.28 ms hop


def stft(iq, fs):
    f, t, Z = signal.stft(iq, fs=fs, window="hann", nperseg=NPERSEG,
                          noverlap=NPERSEG - HOP, return_onesided=False, boundary=None, padded=False)
    order = np.argsort(f)
    return f[order], t, np.abs(Z[order]) ** 2


def main():
    for i in range(1, 5):
        iq, m = load(i)
        fs, fc = m["sample_rate_hz"], m["fc_hz"]
        lam = C / fc
        f, t, P = stft(iq, fs)
        PdB = 10 * np.log10(P + 1e-12)
        floor = np.median(PdB)
        # robust noise sigma per bin (dB)
        mad = np.median(np.abs(PdB - floor))
        pk_idx = P.argmax(axis=0)
        pk_db = PdB[pk_idx, np.arange(P.shape[1])]
        snr = pk_db - floor
        thr = 15  # dB above median
        det = snr > thr
        print(f"\n=== capture{i}  (bin = {fs/NPERSEG:.1f} Hz = {fs/NPERSEG*lam/2:.2f} m/s, hop {HOP/fs*1e3:.2f} ms)")
        print(f"  noise median {floor:.1f} dB, MAD {mad:.2f} dB; detections(>{thr}dB): {det.sum()} / {det.size}")
        idx = np.where(det)[0]
        if idx.size == 0:
            continue
        # contiguous runs
        runs = np.split(idx, np.where(np.diff(idx) > 20)[0] + 1)
        print("  runs:", [(round(t[r[0]], 3), round(t[r[-1]], 3), len(r)) for r in runs])
        r = max(runs, key=len)
        tt = t[r]; fd = f[pk_idx[r]]; v = fd * lam / 2
        print(f"  main track: t={tt[0]:.3f}..{tt[-1]:.3f}s ({tt[-1]-tt[0]:.3f}s)")
        print(f"  v_start={v[0]:.2f} m/s  v_max={v.max():.2f} @ t={tt[v.argmax()]:.3f}  v_end={v[-1]:.2f} m/s  ({v.max()*3.6:.1f} km/h max)")
        # print coarse track every 50 ms
        for tq in np.arange(tt[0], tt[-1], 0.05):
            k = np.argmin(abs(tt - tq))
            print(f"    t={tt[k]:.3f}  v={v[k]:6.2f} m/s  snr={snr[r][k]:5.1f} dB")
        np.save(f"track{i}.npy", np.c_[tt, v, snr[r]])


if __name__ == "__main__":
    main()

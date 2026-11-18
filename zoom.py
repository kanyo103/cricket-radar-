"""Step 3: finer analysis. Zoomed spectrograms, sub-bin ridge tracking with a long window, tail structure."""
import numpy as np
from scipy import signal
from analyze import load, C
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def ridge(iq, fs, lam, nper, hop, t0, t1, vmin=0, vmax=60, min_snr=8):
    """Sub-bin ridge track within [t0,t1] and velocity band [vmin,vmax]. Returns t, v, snr(dB)."""
    f, t, Z = signal.stft(iq, fs=fs, window="hann", nperseg=nper, noverlap=nper - hop,
                          return_onesided=False, boundary=None, padded=False)
    o = np.argsort(f); f = f[o]; P = np.abs(Z[o]) ** 2
    v = f * lam / 2
    band = (v >= vmin) & (v <= vmax)
    floor = np.median(P)
    out = []
    for k in range(P.shape[1]):
        if not (t0 <= t[k] <= t1):
            continue
        col = np.where(band, P[:, k], 0)
        j = col.argmax()
        snr = 10 * np.log10(P[j, k] / floor)
        if snr < min_snr or j in (0, len(f) - 1):
            continue
        a, b, c = np.log(P[j - 1, k]), np.log(P[j, k]), np.log(P[j + 1, k])
        d = 0.5 * (a - c) / (a - 2 * b + c)
        fv = f[j] + d * (f[1] - f[0])
        out.append((t[k], fv * lam / 2, snr))
    return np.array(out)


def main():
    ev = {1: (0.2, 1.5), 2: (0.45, 1.56), 3: (1.7, 2.58), 4: (1.2, 1.96)}
    fig, axs = plt.subplots(4, 1, figsize=(14, 16))
    for ax, i in zip(axs, range(1, 5)):
        iq, m = load(i)
        fs, fc = m["sample_rate_hz"], m["fc_hz"]; lam = C / fc
        t0, t1 = ev[i]
        s0, s1 = int(t0 * fs), int(t1 * fs)
        f, t, S = signal.spectrogram(iq[s0:s1], fs=fs, nperseg=1024, noverlap=960, window="hann",
                                     return_onesided=False, mode="psd")
        f = np.fft.fftshift(f); S = np.fft.fftshift(S, axes=0)
        v = f * lam / 2
        sel = (v > 0) & (v < 60)
        SdB = 10 * np.log10(S[sel] + 1e-12)
        ax.pcolormesh(t + t0, v[sel], SdB, shading="auto", vmin=np.median(SdB), vmax=SdB.max())
        r = ridge(iq, fs, lam, 2048, 256, t0, t1)
        if len(r):
            ax.plot(r[:, 0], r[:, 1], "r.", ms=2, alpha=.6)
        ax.set_title(f"capture{i} zoom (red = ridge, 2048-pt window)"); ax.set_ylabel("v [m/s]")
    axs[-1].set_xlabel("time [s]")
    fig.savefig("zoom.png", dpi=80, bbox_inches="tight")

    for i in range(1, 5):
        iq, m = load(i); fs, fc = m["sample_rate_hz"], m["fc_hz"]; lam = C / fc
        r = ridge(iq, fs, lam, 2048, 512, *ev[i])
        print(f"\ncapture{i}: ridge points={len(r)}  span {r[0,0]:.3f}..{r[-1,0]:.3f}s")
        for row in r[::max(1, len(r) // 22)]:
            print(f"   t={row[0]:.3f}  v={row[1]:6.2f}  snr={row[2]:5.1f}")


if __name__ == "__main__":
    main()

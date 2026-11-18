"""Step 1: load the captures, validate against metadata, and look at basic signal quality + spectrograms."""
import json, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal

C = 299_792_458.0


def load(i):
    d = json.load(open(f"capture{i}.json"))
    meta = d["metadata"]["/"]["attrs"]
    n = d["metadata"]["/iq"]["sample_count"]
    raw = np.array(d["data"]["/iq"], dtype=np.float64)
    assert raw.size == 2 * n, (raw.size, n)
    iq = raw[0::2] + 1j * raw[1::2]
    return iq, meta


if __name__ == "__main__":
    for i in range(1, 5):
        iq, m = load(i)
        fs, fc = m["sample_rate_hz"], m["fc_hz"]
        lam = C / fc
        print(f"\n=== capture{i}: N={iq.size}  dur={iq.size/fs:.3f}s  fs={fs}  fc={fc/1e9}GHz  lambda={lam*1000:.2f}mm")
        print(f"  I mean={iq.real.mean():.2f} std={iq.real.std():.2f} min/max={iq.real.min():.0f}/{iq.real.max():.0f}")
        print(f"  Q mean={iq.imag.mean():.2f} std={iq.imag.std():.2f} min/max={iq.imag.min():.0f}/{iq.imag.max():.0f}")
        print(f"  I/Q corr={np.corrcoef(iq.real, iq.imag)[0,1]:.3f}  |iq| max={np.abs(iq).max():.1f}")

        f, t, S = signal.spectrogram(iq, fs=fs, nperseg=1024, noverlap=768,
                                     window="hann", return_onesided=False, mode="psd")
        f = np.fft.fftshift(f); S = np.fft.fftshift(S, axes=0)
        SdB = 10 * np.log10(S + 1e-12)
        floor = np.median(SdB)
        print(f"  spectrogram: median={floor:.1f} dB, max={SdB.max():.1f} dB, peak above floor={SdB.max()-floor:.1f} dB")
        fig, ax = plt.subplots(figsize=(14, 5))
        v = f * lam / 2
        im = ax.pcolormesh(t, v, SdB, shading="auto", vmin=floor, vmax=SdB.max())
        ax.set_xlabel("time [s]"); ax.set_ylabel("radial velocity [m/s] (+ = away)")
        ax.set_title(f"capture{i}"); fig.colorbar(im, label="dB")
        fig.savefig(f"spec{i}.png", dpi=90, bbox_inches="tight"); plt.close(fig)

"""Check whether the noise-only regions are independent between captures (real hardware) or reused (simulation)."""
import numpy as np
from analyze import load

iqs = [load(i)[0] for i in range(1, 5)]
N = 10000  # first 0.2 s: before any event in all captures
for a in range(4):
    for b in range(a + 1, 4):
        x, y = iqs[a][:N], iqs[b][:N]
        cI = np.corrcoef(x.real, y.real)[0, 1]
        cQ = np.corrcoef(x.imag, y.imag)[0, 1]
        print(f"cap{a+1} vs cap{b+1}  corr I={cI:+.3f}  corr Q={cQ:+.3f}")

# How many leading samples of I match exactly / closely?
print()
for a in range(1, 4):
    d = np.abs(iqs[0].real[:2000] - iqs[a].real[:2000])
    print(f"cap1 vs cap{a+1}: mean |dI| over first 2000 = {d.mean():.2f}; frac |dI|<=3: {(d<=3).mean():.2f}")
# Autocorrelation / spectrum flatness of the noise
x = iqs[0][:20000]
S = np.abs(np.fft.fft(x * np.hanning(x.size))) ** 2
S = np.fft.fftshift(S)
print("\nnoise PSD flatness (std of 20 band means, dB):",
      np.std(10 * np.log10([b.mean() for b in np.array_split(S, 20)])).round(2))
print("kurtosis I:", ((x.real - x.real.mean()) ** 4).mean() / x.real.var() ** 2)

# Cricket radar prototype: what the four captures show

## Data
Four JSON files, each one continuous int8 I/Q stream. The metadata gives carrier 24 GHz, sample rate 50 kHz, and `velocity_sign = positive_away`.

| File | Samples | Duration |
|---|---|---|
| capture1 | 75,854 | 1.52 s |
| capture2 | 78,281 | 1.57 s |
| capture3 | 128,944 | 2.58 s |
| capture4 | 98,001 | 1.96 s |

Doppler to radial velocity: v = f·λ/2, with λ = c/fc = 12.49 mm. The unambiguous range is ±156 m/s. A 2048-point window gives about 0.15 m/s resolution.

Sanity checks passed: I and Q are balanced (std about equal, no DC offset, correlation about 0) and the sample counts match the metadata.

## What the data contains
Each capture is flat noise, then a single narrow Doppler line appears abruptly. Method: spectrogram, ridge tracking with parabolic sub-bin interpolation, onset taken at the first frame above 20 dB SNR.

| Capture | Onset | Peak radial speed | Sudden drop | Time after onset | Radial path to drop |
|---|---|---|---|---|---|
| 1 | 0.26 s | 19.2 m/s (69 km/h) | 18.8 to 16.4 m/s (x0.87) | 0.80 s | 15.1 m |
| 2 | 0.51 s | 23.2 m/s (83 km/h) | 22.8 to 19.9 m/s (x0.87) | 0.58 s | 13.0 m |
| 3 | 1.74 s | 28.8 m/s (104 km/h) | 27.6 to 24.2 m/s (x0.88) | 0.60 s | 16.7 m |
| 4 | 1.22 s | 34.7 m/s (125 km/h) | 32.5 to 29.8 m/s (x0.92) | 0.31 s | 10.1 m |

Observations:
1. Each event starts abruptly and then rises quickly (4 to 12 m/s over about 0.1 s) to a peak.
2. After the peak the speed declines slowly.
3. Each event has one sudden drop of 2.4 to 3.4 m/s (8 to 13%), completed within the 5 ms frame spacing, after which the line continues at lower speed.
4. Capture 4 has a second parallel track about 3 to 4 m/s above the main one. It peaks at 38.4 m/s and ends around t = 1.43 s. I have no explanation for it.
5. SNR falls by about 30 dB over about 20 m of integrated radial path. A receding point target fits SNR = K - 10·n·log10(R) with fitted n of 3.7 to 4.8 (about 1 dB RMS), consistent with the radar equation (R^-4). This supports "away" independently of the metadata sign flag.

## Interpretation (hypotheses, not measurements)
- **Four cricket deliveries**, moving away from the radar.
- **The sudden drop is the bounce.** A 12% speed loss in under 5 ms is typical of a ball hitting the pitch. Release-to-bounce path of 10 to 17 m is plausible.
- **Speeds are lower bounds.** The radar measures only the radial component, v·cos(theta). The early rise is too fast to be real acceleration (4 to 12 m/s in about 0.1 s), so it is more likely geometry: the angle between the ball's path and the radar line of sight changing.
- Capture 4's second track could be a second object or a second scattering component. Unverified.

## What the radar could do
Measure radial speed continuously, detect delivery and bounce, and measure speed lost at bounce.

It could not give range, position, lateral movement or spin from these data. Range change can only be estimated by integrating velocity, which needs a starting range and accumulates error. I did not check for chirp or frame structure. If the prototype is FMCW, range information may exist in a form I did not use.

## Concern about the data
In the first 0.2 s of every file, before any event, the I channel is the same noise sequence in all four captures (correlation 0.999, with only a gain change between files). Q is independent. Real receiver noise should not repeat. The noise is also very Gaussian (kurtosis about 3.0) and spectrally flat. Either the data is synthetic with a fixed seed, or there is a fixed-pattern fault on I. The files cannot tell me which. I would ask the team about this before trusting the prototype's noise performance. Also, "It works!" is not backed by any ground truth: no reference speeds were supplied, so accuracy is unverified.

## Limits of this analysis
- Rise-time and slope numbers come from few points. Treat them as rough. The capture 4 slope before the drop is distorted by the second track.
- Capture-to-capture speed differences (69 to 125 km/h) may reflect geometry or different bowlers, not necessarily different true speeds.
- Output of the full run is in `findings_output.txt`.

## Reproduce
```
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python analyze.py       # loads data, validates, writes spec*.png
.venv/bin/python noise_check.py   # I/Q noise correlation between captures
.venv/bin/python track.py         # coarse ridge tracking
.venv/bin/python zoom.py          # fine ridge tracking, writes zoom.png
.venv/bin/python findings.py      # drop, rise, second track, SNR vs range
```
Files: `analyze.py` (loader), `noise_check.py`, `track.py`, `zoom.py`, `findings.py`, figures `spec1-4.png` and `zoom.png`.

## AI use
I used an AI coding assistant (Kiro) as a collaborator. I framed the problem and told it not to assume anything about the data, so it inspected the
files before analysing them and tested its conclusions instead of taking them as given.I chose how analyze.py loads the captures. I read the metadata (interleaved_iq,
order real then imag), split the flat array into I and Q, and checked the sample counts. It wrote the Python scripts and first drafts of the findings, and I ran
them and directed the steps. I reviewed the results and decided which conclusions to keep. Interpretations that go beyond the data, such as the bounce
explanation, are marked as hypotheses.

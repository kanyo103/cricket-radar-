# Cricket radar capture analysis

Analysis of four I/Q captures from a 24 GHz cricket radar prototype. The full write-up is in **[REPORT.md](REPORT.md)**.

## Data note
The capture files (`capture1.json` to `capture4.json`) and the task instructions (`Instructions.txt`) were supplied by the company as part of the application. They are **not included in this repository**.

To rerun the analysis, place the four capture files in the repository root, next to the scripts. The scripts expect this layout:

```
capture1.json  capture2.json  capture3.json  capture4.json
analyze.py  noise_check.py  track.py  zoom.py  findings.py
```

Each capture file has this structure (taken from its metadata):
- `/iq`: int8 samples, interleaved real then imaginary, shape `[N, 2]`
- `fc_hz` = 24,000,000,000, `sample_rate_hz` = 50,000, `velocity_sign` = `positive_away`

## Run
```
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python analyze.py       # loads and validates data, writes spec*.png
.venv/bin/python noise_check.py   # I/Q noise correlation between captures
.venv/bin/python track.py         # coarse ridge tracking
.venv/bin/python zoom.py          # fine ridge tracking, writes zoom.png
.venv/bin/python findings.py      # sudden speed drop, rise, second track, SNR vs range
```
`findings_output.txt` is the saved output of `findings.py`. `spec1.png` to `spec4.png` and `zoom.png` are the figures from the report.

## Files
| File | Purpose |
|---|---|
| `REPORT.md` | Findings, interpretation, limitations, AI use |
| `analyze.py` | Loader (`load()`), basic signal checks, spectrograms |
| `noise_check.py` | Compares the noise-only start of each capture |
| `track.py` | Coarse Doppler ridge tracking |
| `zoom.py` | Fine ridge tracking (`ridge()`), zoomed spectrograms |
| `findings.py` | Quantifies the speed drop, rise, second track, SNR fall |

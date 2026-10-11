# Shooting protocol for a sunset (or sunrise) sequence

The retrieval uses the **change in sun colour with air mass** within one sequence. Anything that changes the colour for a reason other than the atmosphere (white balance, tone curve, clipping, clouds) is an error. This protocol removes those causes.

## Before you go

- [ ] Camera can shoot **RAW** (any interchangeable-lens camera; many phones support DNG).
- [ ] Tripod, plus an ND filter or solar film if the sun is still high (start of sequence).
- [ ] Phone with automatic network/GPS time, to check the camera clock.
- [ ] Clear view toward the sunset azimuth down to the horizon. Check with any sun-position app.
- [ ] Note site latitude, longitude (decimal degrees) and altitude (m).

## Camera settings (keep them fixed for the whole sequence)

| Setting | Value | Why |
|---|---|---|
| File format | RAW (or RAW + JPEG) | Linear data. JPEG biased slopes by 2–3 % in testing |
| White balance | Daylight (fixed, never Auto) | The retrieval relies on constant channel gains |
| ISO | Lowest native (100–200) | Noise |
| Focus | Manual, on the sun/horizon | Autofocus hunting changes the disk size |
| Lens | 100–300 mm if available | Sun disk ≥ 15 px in radius gives better photometry |
| Image stabilisation | Off on a tripod | Avoids drift |
| Exposure | Manual. Adjust between frames if needed, but **never let the sun disk clip** | Only channel ratios matter; clipping destroys them |

Check the histogram: the sun disk must sit below the right edge in **all three** channels. As the sun reddens, blue gets very dark. That is expected, but if blue drops into the noise, the tool rejects the frame (SNR < 30).

## During the sequence

1. Photograph the phone clock showing seconds. You will compute `--clock-offset-s` = phone time − camera time.
2. Start at about 15° sun elevation (roughly 1 h 15 min before sunset at mid-latitudes; less in summer at high latitude, about 1 h in Lebanon).
3. Take one frame every **30–60 s** until the sun touches the horizon or an obstacle.
4. Keep the sun inside the frame and away from the edges. Reframe if needed (the pipeline finds the sun in each frame).
5. Write down conditions:
   - cloud at the sun's position (yes/no, when)
   - visible haze, smoke smell or dust
   - wind direction
   - anything unusual

For **sunrise**, reverse the order: start at the horizon and continue to about 15°.

## After

```bash
python -m skycolour process <folder> --out results/<site>_<date> \
    --lat <lat> --lon <lon> --tz <IANA zone> --clock-offset-s <seconds> [--altitude-m <m>]
```

Then:

- Open `colour_index_vs_airmass.png`. The points should lie on a straight line. Curvature means the aerosol changed during the sequence, or clouds crossed the sun.
- Read `qc_reason` in `frames.csv` for every rejected frame.
- If mountains or buildings hide the last degrees, set `--airmass-max` to the air mass of the last clear frame.

## Site notes for the planned comparison

| Site | Practical notes |
|---|---|
| Montréal | Winter: cold, dry, very clean (small aerosol slope expected). Summer: humidity and wildfire smoke episodes (large slope, high α). |
| Paris | Frequent cloud: pick clear evenings. Spring nitrate episodes are a good contrast case. |
| Lebanon, coast | Sunsets over the sea give an open horizon. Spring dust (low α) is the case where colour underestimates AOD most. |
| Lebanon, mountains (~1800 m) | Use `--altitude-m`. Often above the coastal haze layer, giving a clean reference. |

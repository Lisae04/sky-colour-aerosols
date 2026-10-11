# data/

Put each photo sequence in its own folder, e.g. `data/raw/montreal_2026-10-09_sunset/`.
`data/raw/` is git-ignored so large RAW files stay off GitHub.
Keep a one-line note per sequence (site, time, conditions) in `data/sessions.csv`:

```
folder,site,lat,lon,altitude_m,tz,clock_offset_s,notes
montreal_2026-10-09_sunset,Montreal,45.5017,-73.5673,40,America/Toronto,0,clear; light haze west
```

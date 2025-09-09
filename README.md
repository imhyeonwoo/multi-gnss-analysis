# multi-gnss-analysis

MATLAB-based project for analyzing multi-GNSS (GPS, Galileo, GLONASS, BeiDou) data.  
Developed as part of experimental studies in guidance and control at Konkuk University.  

## Features
- Least Squares (LS) positioning
- Coordinate transformations: ECEF ↔ ENU/LLH
- Satellite selection based on elevation angle and C/N0
- Skyplot generation
- NovAtel vs. Android GNSS performance comparison
- Export to Google Earth (KML trajectory)

## Directory Structure
```bash
multi-gnss-analysis/
├─ matlab/ # Core MATLAB functions (LS.m, xyz2enu.m, etc.)
├─ data/ # RINEX, NovAtel, Android GNSS/IMU logs
├─ output/ # Figures, KML trajectories
└─ README.md
```

## Getting Started
1. Place raw GNSS data into the `data/` directory.
2. Run `matlab/main.m`.
3. Results will be saved in `output/` (plots and KML files).

## License
MIT

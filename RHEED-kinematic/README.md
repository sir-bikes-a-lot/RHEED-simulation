Simulate reflection high-energy electron diffraction (RHEED) patterns for arbitrary atomic arrangements using the kinematic diffraction approximation.
Based on: Ayahiko Ichimaya and Philip I Cohen, "Reflection High Energy Electron Diffraction", Chapters 4, 10, and 11, Cambridge University Press (2004).

Function "CalcRHEED()" returns the screen coordinates and diffraction intensity for the VASP POSCAR or molecular dynamics dump file specified in "filename".
Additionally, the crystal parameters and Kikuchi line locations are returned.

Function "PlotRHEED()" plots the RHEED pattern in greyscale, and returns the broadened and/or clipped intensity pattern.

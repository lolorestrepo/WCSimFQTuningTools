`compute_STable.py` Computes the Scattering Table from simulation files. 

Run it by: `python compute_STable.py /path/to/files [--nbins nzs nRs nR_PMT nz_PMT nphi nzd ntheta] [--zedge zedge in cm] [--redge redge in cm] [--vaxis 0|1|2 (default 2)] [--wcsimlib /path/to/wcsimlib/] [--fitqun /path/to/fitqun/source/code/] [-v]`

The `/path/to/files` argument point to the directory of the simulation files. The script ignores files ending with `_flat.root`.

The `nzs nRs nR_PMT nz_PMT nphi nzd ntheta` arguments are the number of histogram bins for each variable (see defaults in script):
```
zs: source vertical position
Rs: source R position
R_PMT: PMT radial position for PMTs at the bottom and top caps
z_PMT: PMT vertical position for PMTs at the side
phi: azimuth angle of source w.r.t. PMT position
zd: source vertical propagation direction 
theta: source azimuth propagation direction
```

The `vaxis` variable is used to indicate the vector component of the vertical axis of the detector 0, 1, 2 (default). For WCTE, the vertical component is the Y axis, ie vector component number 1.

The `zedge` and `redge` variables indicate the histogram bounds in cm. The vertical dimensions `(zs, z_PMT)` histogram will be produced between `[-zedge, zedge]` and the horizontal dimensions `(Rs, R_PMT)` between `[0, redge]`.

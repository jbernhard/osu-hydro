# OSU hydro

_The Ohio State University viscous hydrodynamics code for relativistic heavy-ion collisions_

## Attribution

I am not the author of this software:
it was originally written by Dr. Huichao Song and further modified and enhanced by Zhi Qiu, Chun Shen, and Jia Liu.
Alternate versions may be found at https://github.com/chunshen1987/VISHNew and https://github.com/JiaLiu11/iEBE/tree/hydro/EBE-Node/VISHNew.

Some references:

- https://inspirehep.net/record/771325
- https://inspirehep.net/record/785597

## Usage

This is a stripped down and cleaned up version of the Ohio State University viscous 2+1D hydrodynamics code (also known elsewhere as "VISHNew").
I started from Jia Liu's version including shear and bulk viscosity at https://github.com/JiaLiu11/iEBE/tree/hydro/EBE-Node/VISHNew, however this repository does not share git history with Jia's because I used `git filter-branch` to remove a number of large files from old commits.

I made this version with one task in mind: running hydro as part of an event-by-event hybrid model, i.e. to be preceded by a Monte Carlo initial condition model and followed by particlization and hadronic transport models.

To compile and install, follow the standard CMake sequence:

    mkdir build && cd build
    cmake .. [-DCMAKE_INSTALL_PREFIX=<prefix>]
    make install

This will place the compiled binary `osu-hydro` in `<prefix>/bin/` and the configuration file `osu-hydro.conf` and equation of state data `eos.dat` in `<prefix>/share/osu-hydro/`.

At runtime, `osu-hydro` searches for `osu-hydro.conf` and `eos.dat` in the following locations:

1. The current directory.
2. User data directory `$XDG_DATA_HOME/osu-hydro/`, with `$XDG_DATA_HOME` defaulting to `~/.local/share/` if unset.
3. System directory `/usr/local/share/osu-hydro/`.
4. System directory `/usr/share/osu-hydro/`.

If either of the files cannot be found, the program will exit with an error message.

Thus, `osu-hydro` may be executed from anywhere provided it can locate the data files.
However, it always expects the initial condition data files (described below) to be in the current working directory, and it will also write output files in the current directory.

### Equation of state

The EOS is generated at build time by the script `eos/eos.py`, which connects a hadron resonance gas EOS at low temperature to a lattice EOS (HotQCD, http://inspirehep.net/record/1307761) at high temperature.
See the help output (`eos.py --help`) for more options.
This depends on [frzout](https://github.com/jbernhard/frzout) to compute the hadron resonance gas part.

### Configuration

The config file `osu-hydro.conf` has brief comments for each parameter.
Note that the code expects the parameters in a specific order, so __do not change the order of parameters in the config file__.
Parameters may also be passed on the command line with `key=value` syntax (case-insensitive), e.g.

    osu-hydro edec=0.30

### Temperature-dependent viscosities

The shear viscosity is parametrized as

    (eta/s)(T) = min + slope*(T - Tc) * (T/Tc)^curvature

for T > Tc = 0.154 GeV, where `min`, `slope`, and `curvature` are input parameters.
The `min` and `slope` must be non-negative; `curvature` may be negative but probably not below -1, since this would cause decreasing (eta/s)(T).

For T < Tc (HRG phase), eta/s is constant, controlled by a single input parameter.

These parameters may be set in the config file or on the command line with keys `etas_hrg`, `etas_min`, `etas_slope`, `etas_curv`.

The bulk viscosity is parametrized as a Cauchy distribution with tunable peak location, max value, and width:

    (zeta/s)(T) = max / (1 + ((T - T0)/width)^2)

The location (temperature T0), max, and width may be set in the config file or on the command line with keys `zetas_t0`, `zetas_max`, `zetas_width`.

### Grid settings

The computational grid is always square and centered at the origin.
The size may be set in the config file or on the command line by the `nls` parameter.
`LS` ("lattice size") means the number of cells from the origin to the grid edge, not counting the cell precisely at the origin; the full grid is (2\*LS + 1) × (2\*LS + 1).

The grid cell size is set in the config file or on the command line by the `dxy` parameter.
Ensure the cell size is small enough to resolve the system geometry.
For Monte Carlo initial conditions with Gaussian nucleons of width `w`, `dxy = 0.15*w` is a good compromise between computation time and precision.

__Example__: With the default `DXY = 0.1 fm` and `LS = 200`, the grid is 401 × 401 with cells centered at -20.0, -19.9, ..., 0.0, 0.1, ..., 20.0 fm.

__WARNING__:
_The grid size is defined by the **center** of the outermost grid cell, which is inconsistent with programs that use the cell **edge**, such as [trento](https://github.com/Duke-QCD/trento) and [freestream](https://github.com/Duke-QCD/freestream).
The two definitions differ by half a cell width.
For example, the default 401 × 401 grid has its outermost cells centered at ±20.0 fm, and the cells are 0.1 fm wide, so the physical extent of the grid is a half-cell larger: ±20.05 fm.
To generate a matching grid with `trento`, use options `--grid-max 20.05 --grid-step 0.1`._

The timestep is set in the config file or on the command line by the `dt` parameter.
The SHASTA algorithm used in this code requires the timestep to be less than one-half the spatial step.
However, a smaller timestep around one-quarter to one-fifth the spatial step is a good compromise between computation time and precision.

### Running initial conditions

By default (with option `InitialURead = 1`), the initial energy density, flow, and shear viscous tensor must be provided.
Create the following __binary__ data files:

- `ed.dat` - energy density
- `u{1,2}.dat` - flow velocity in x and y directions
- `pi{11,12,22}.dat` - shear tensor components xx, xy, yy (the remaining components are computed internally)
- The bulk pressure is computed internally as the deviation from ideal pressure, Π = e/3 - p(e).

If `InitialURead = 0`, only the energy density is required.
Alternatively, set `IEin = 1` and provide the entropy density `sd.dat`.

All files must contain an array of doubles in binary format.
They can be easily created in Python using [numpy.ndarray.tofile](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.tofile.html), e.g., if `ed` is a numpy array containing the energy density, then `ed.tofile('ed.dat')` writes it as a binary file.

__WARNING__:
_The input data files must have precisely the expected grid size, as described above in [grid settings](#grid-settings).
This is not checked by the code._

During time evolution, the code prints out a line for each timestep

     timestep_number  tau  max_energy_density  max_temperature

Evolution stops after the max energy density drops below the configured decoupling energy density.

The only output file is the hypersurface data `surface.dat`, a __binary__ file containing an N × 16 array of doubles, where N is the number of freeze-out volume elements and the 16 columns are

Index | Quantity                        | Index | Viscous pressures [GeV/fm<sup>3</sup>]
------|---------------------------------|-------|---------------------------------------
0     | τ [fm]                          | 8     | π<sup>tt</sup>
1     | x [fm]                          | 9     | π<sup>tx</sup>
2     | y [fm]                          | 10    | π<sup>ty</sup>
3     | dσ<sub>t</sub> [fm<sup>2</sup>] | 11    | π<sup>xx</sup>
4     | dσ<sub>x</sub> [fm<sup>2</sup>] | 12    | π<sup>xy</sup>
5     | dσ<sub>y</sub> [fm<sup>2</sup>] | 13    | π<sup>yy</sup>
6     | v<sub>x</sub>                   | 14    | π<sup>zz</sup>
7     | v<sub>y</sub>                   | 15    | Π

The file can be read in Python with

```python
surface = np.fromfile('surface.dat', dtype='float64').reshape(-1, 16)
```

__Note__:
The surface normal vectors dσ<sub>μ</sub> are output in covariant form, which is standard.
From [84fc226](https://github.com/jbernhard/osu-hydro/commit/84fc2261ead02a690b12424e71a44584b91c01e1) (July 2016) until [9995ee1](https://github.com/jbernhard/osu-hydro/commit/9995ee157c142c3a9fb76eba4a2fe4913d4770a7) (April 2017), they were contravariant (dσ<sup>μ</sup>).

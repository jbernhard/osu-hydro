# vishnew

_viscous 2+1D Israel-Stewart hydrodynamics for relativistic heavy-ion collisions_

## Attribution

I am not the author of this software:
it was originally written by Dr. Huichao Song and further modified and enhanced by Zhi Qiu, Chun Shen, and Jia Liu.
Alternate versions may be found at https://github.com/chunshen1987/VISHNew and https://github.com/JiaLiu11/iEBE/tree/hydro/EBE-Node/VISHNew.

Some references:

- http://inspirehep.net/record/771325
- http://inspirehep.net/record/785597

## Usage

This is a stripped down and cleaned up version of the Ohio State University viscous 2+1D hydrodynamics code.
I started from Jia Liu's version including shear and bulk viscosity at https://github.com/JiaLiu11/iEBE/tree/hydro/EBE-Node/VISHNew, however this repository does not share git history with Jia's because I used `git filter-branch` to remove a number of large files from old commits.

I made this version with one task in mind: running hydro as part of an event-by-event hybrid model, i.e. to be preceded by a Monte Carlo initial condition model and followed by particlization and hadronic transport models.

To compile and install, follow the standard CMake sequence:

    mkdir build && cd build
    cmake ..
    make install

This will place several files in `<prefix>/vishnew`: the compiled binary `vishnew`, configuration file `vishnew.conf`, and equation of state table `eos.dat`.

### Equation of state

The EOS table is generated at build time by the script `eos/generate-eos-table.py`, which connects a hadron resonance gas EOS at low temperature to a lattice EOS (HotQCD, http://inspirehep.net/record/1307761) at high temperature.
Run `./generate-eos-table.py --plot` to visualize several important thermodynamic quantities.

The EOS table has columns

    energy_density  pressure  entropy_density  temperature

in natural units of GeV and fm.
It has 100000 rows with evenly-spaced energy density steps, as expected by the `vishnew` code.

### Configuration

The config file `vishnew.conf` has brief comments for each parameter.
Note that the code expects the parameters in a specific order, so __do not change the order of parameters in the config file__.
Parameters may also be passed on the command line with a `key=value` syntax, e.g.

    vishnew Edec=0.30

### Temperature-dependent viscosities

The shear viscosity is parametrized as

    (eta/s)(T) = min + slope*(T - Tc) * (T/Tc)^curvature

for T > Tc = 0.154 GeV, where `min`, `slope`, and `curvature` are input parameters.
The `min` and `slope` must be non-negative; `curvature` may be negative but probably not below -1, since this would cause decreasing (eta/s)(T).

For T < Tc (HRG phase), eta/s is constant, controlled by a single input parameter.

These parameters may be set in the config file or on the command line with keys `etas_hrg`, `etas_min`, `etas_slope`, `etas_curv`.

The bulk viscosity is parametrized as a Cauchy distribution with its peak at Tc and tunable max and width:

    (zeta/s)(T) = max / (1 + ((T - Tc)/width)^2)

The max and width may be set in the config file or on the command line with keys `zetas_max`, `zetas_width`.

### Running initial conditions

By default (with option `InitialURead = 1`), the initial energy density, flow, and shear viscous tensor must be provided.
Place the following files in the `vishnew` folder:

- `ed.dat` - energy density
- `u{1,2}.dat` - flow velocity in x and y directions
- `pi{11,12,22}.dat` - components of the shear tensor (the remaining components are computed internally)
- The bulk pressure is computed internally as the deviation from ideal pressure, Π = e/3 - p(e).

If `InitialURead = 0`, only the energy density is required.
Alternatively, set `IEin = 1` and provide the entropy density `sd.dat`.

All files must contain square block-style grids and nothing else (comments will not work).
The grid size depends on the `LS` (lattice size) parameter: grid size = 2\*LS + 1.
The grid step is hard-coded to 0.1 fm.
So e.g. with the default `LS = 130`, the grid is 261x261 from -13 to +13 fm.

During time evolution, the code prints out a line for each timestep

     timestep_number  tau  max_energy_density  max_temperature

Evolution stops after the max energy density drops below the configured decoupling energy density.

The only output file is the hypersurface data `surface.dat`, beginning with a commented header containing the thermodynamic freeze-out quantities

    # e = <energy density>
    # p = <pressure>
    # s = <entropy density>
    # T = <temperature>

and followed by data columns

    tau x y dsigma^t dsigma^x dsigma^y v_x v_y pi00 pi01 pi02 pi11 pi12 pi22 pi33 PI

where each row represents a volume element.
All units are GeV and fm.

**Note**:
The normal vector `dsigma` is output in contravariant form (index raised).
Previous versions output the covariant vector (index lowered).

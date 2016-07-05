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

The EOS table is generated at build time by the script `eos/generate-eos-table.py`, which blends a hadron resonance gas EOS at low temperature into the HotQCD EOS (http://inspirehep.net/record/1307761) at high temperature.
Run `./generate-eos-table.py --plot` to visualize several important thermodynamic quantities.

The EOS table has columns

    energy_density  pressure  entropy_density  temperature

in natural units of GeV and fm.
It has exactly 155500 rows with evenly-spaced energy density steps, as expected by the `vishnew` code.

### Configuration

The config file `vishnew.conf` has brief comments for each parameter.
Note that the code expects the parameters in a specific order, so __do not change the order of parameters in the config file__.
Parameters may also be passed on the command line with a `key=value` syntax, e.g.

    vishnew Edec=0.30

### Running initial conditions

Place a file `initial.dat` in the `vishnew` folder and run the executable.

The initial file must be a square block-style grid of energy or entropy density, with no other content in the file (comments will not work).
The grid size depends on the `LS` (lattice size) parameter: grid size = 2\*LS + 1.
The grid step is hard-coded to 0.1 fm.
So e.g. with the default `LS = 130`, the initial condition must be a 261x261 grid from -13 to +13 fm.

If the `IEin` parameter is set to 1, the initial condition will be interpreted as entropy density;
if `IEin = 0`, it will be interpreted as energy density.

During time evolution, the code prints out a line for each timestep

     timestep_number  tau  max_energy_density  max_temperature

Evolution stops after the max energy density drops below the configured decoupling energy density.

The only output file is the hypersurface data `surface.dat`, beginning with a commented header containing the thermodynamic freeze-out quantities

    # e = <energy density>
    # p = <pressure>
    # s = <entropy density>
    # T = <temperature>

and followed by data columns

    tau x y dsigma_t dsigma_x dsigma_y v_x v_y pi00 pi01 pi02 pi11 pi12 pi22 pi33 PI

where each row represents a volume element.
All units are GeV and fm.

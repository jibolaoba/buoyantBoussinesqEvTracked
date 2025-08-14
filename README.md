# buoyantBoussinesqEvTracked

## Overview

This repository contains a custom CFD (Computational Fluid Dynamics) solver based on the OpenFOAM framework, specifically designed for simulating buoyancy-driven flows using the Boussinesq approximation with evaporation tracking. The solver, named `buoyantBoussinesqEvTracked`, extends standard OpenFOAM solvers to incorporate phase change phenomena, such as evaporation at fluid interfaces, while accounting for buoyancy effects in incompressible flows. It is particularly useful for applications involving natural convection, heat transfer, and multiphase flows with evaporation, such as in environmental engineering, HVAC systems, or industrial processes like boiling and condensation.

The solver integrates the Boussinesq approximation for buoyancy (treating density variations as linear with temperature in the momentum equation) and implements a simplified evaporation modeling, that bypass chemical composition effects at the liquid–vapor interface and during evaporation-driven mass transfer.
This project is built on OpenFOAM (Open Field Operation and Manipulation), an open-source CFD toolbox.

## Features

- **Buoyancy Modeling**: Implements the Boussinesq approximation for efficient simulation of buoyancy-driven flows without solving full compressible equations.
- **Evaporation Tracking**: Includes models for phase change at interfaces, tracking evaporation rates based on temperature, concentration, or saturation conditions.
- **Multiphase Support**: Handles two-phase flows (e.g., liquid-gas) with interface capturing or tracking techniques.
- **Heat Transfer**: Coupled energy equation solving for temperature fields influencing buoyancy and evaporation.
- **Customizable**: Easily extendable for specific boundary conditions, turbulence models (e.g., k-epsilon, LES), or additional physics.
- **Parallel Computing**: Compatible with OpenFOAM's parallel processing for large-scale simulations.
- **Post-Processing**: Integrated with ParaView or OpenFOAM utilities for visualization and analysis.

## Dependencies

- **OpenFOAM**: Version 8 or later (tested on v2112 and v2206). Install from [official sources](https://openfoam.org/download/).
- **Compiler**: GCC 7+ or compatible (included in OpenFOAM installation).
- **Optional Libraries**:
  - Swak4Foam or funkySetFields for advanced field manipulations.
  - MPI for parallel runs (e.g., OpenMPI).
- No additional external libraries are required beyond standard OpenFOAM.

Ensure your environment is set up with OpenFOAM sourced (e.g., `source /opt/openfoam8/etc/bashrc`).

## Installation

1. **Clone the Repository**:
   ```
   git clone https://github.com/jibolaoba/buoyantBoussinesqEvTracked.git
   cd buoyantBoussinesqEvTracked
   ```

2. **Compile the Solver**:
   - Navigate to the solver directory (e.g., `applications/solver`).
   - Run the compilation script:
     ```
     wmake
     ```
   - This builds the executable `buoyantBoussinesqEvTracked` in your OpenFOAM user bin directory.

3. **Verify Installation**:
   - Check if the solver is available:
     ```
     buoyantBoussinesqEvTracked -help
     ```
   - It should display usage options.

If compilation fails, ensure your OpenFOAM environment is correctly set up and check for missing dependencies.

## Repository Structure

- **`Allwmake`**: Script to compile all components (solvers, libraries).
- **`applications/`**: Contains the main solver source code.
  - `solver/buoyantBoussinesqEvTracked/`: Core solver files.
    - `buoyantBoussinesqEvTracked.C`: Main executable entry point.
    - `createFields.H`: Field initialization (velocity, pressure, temperature, phase fraction).
    - `UEqn.H`: Momentum equation with Boussinesq term.
    - `pEqn.H`: Pressure correction.
    - `TEqn.H`: Energy/temperature equation.
    - `evaporationModel.H`: Custom evaporation source terms.
- **`cases/`**: Example simulation cases.
  - `tutorialCase/`: A sample case for a buoyant plume with evaporation.
    - `0/`: Initial conditions (U, p, T, alpha.phase).
    - `constant/`: Physical properties, transport models.
    - `system/`: Control dictionaries (controlDict, fvSchemes, fvSolution).
- **`libs/`**: Custom libraries for models (e.g., evaporation models, buoyancy terms).
- **`doc/`**: Documentation files (if any additional PDFs or notes).
- **`README.md`**: This file.
- **`.gitignore`**: Standard ignores for temporary files.

## Usage

### Running a Simulation

1. **Prepare a Case Directory**:
   - Copy an example case: `cp -r cases/tutorialCase myCase`.
   - Navigate to `myCase`.

2. **Set Up the Case**:
   - Edit `system/controlDict` for run time, write intervals, etc.
   - Define boundary conditions in `0/` directory.
   - Specify models in `constant/transportProperties` (e.g., Prandtl number, beta for Boussinesq).
   - For evaporation based on relative humidity, set parameters like Specific heat capacity, beta, Thermal expansion coefficient, Evaporation coefficient and mass transfer coefficients.

3. **Mesh the Domain**:
   - Run `blockMesh` or use snappyHexMesh for complex geometries.

4. **Execute the Solver**:
   ```
   buoyantBoussinesqEvTracked
   ```
   - For parallel run: `decomposePar`, then `mpirun -np * buoyantBoussinesqEvTracked -parallel` (where * represents the number of multiple sub-domains for parallel execution across processors), followed by `reconstructPar`.

5. **Post-Process**:
   - Use `paraFoam` to visualize in ParaView.
   - Compute additional fields with `postProcess` utilities.

### Example: Buoyant Plume with Evaporation

- **Setup**: A hot liquid pool evaporating into air, with buoyancy driving the flow.
- **Key Files**:
  - `constant/transportProperties`: Define phases, viscosity, thermal expansion coefficient (beta), reference temperature.
  - `constant/evaporationProperties`: Evaporation model parameters (e.g., heat of vaporization).
- **Run**: Execute the solver and monitor residuals in the log file.

## Troubleshooting

- **Compilation Errors**: Check OpenFOAM version compatibility. Clean with `wclean` and retry.
- **Runtime Issues**: Ensure mesh quality (`checkMesh`). Adjust solver tolerances in `system/fvSolution`.
- **Evaporation Not Triggering**: Verify parameters such as relative humidity, set parameters like Specific heat capacity, beta, Thermal expansion coefficient, Evaporation coefficient and mass transfer coefficients matches experimentally proven results.

## Contributing

Contributions are welcome! Please fork the repository, make changes, and submit a pull request. Ensure code follows OpenFOAM coding standards. Add tests or example cases for new features.

## License

This project is licensed under the GNU General Public License v3 (GPLv3). See the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Built upon the OpenFOAM foundation.
- Inspired by standard solvers like `buoyantBoussinesqPimpleFoam` and multiphase extensions like `interFoam`.

For questions, contact the repository owner or open an issue.

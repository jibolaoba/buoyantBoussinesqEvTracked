# buoyantBoussinesqEvTracked

## Overview

This repository contains a custom CFD solver based on the OpenFOAM framework, specifically designed for simulating buoyancy-driven flows using the Boussinesq approximation with evaporation tracking. The solver, named `buoyantBoussinesqEvTracked`, extends standard OpenFOAM solvers to incorporate phase change phenomena, such as evaporation at fluid interfaces, while accounting for buoyancy effects in incompressible flows. It is particularly useful for applications involving natural convection, heat transfer, and multiphase flows with evaporation.

The solver integrates the Boussinesq approximation for buoyancy (treating density variations as linear with temperature in the momentum equation) and implements a simplified evaporation modeling that bypasses chemical composition effects at the liquid–vapor interface and during evaporation-driven mass transfer.

## Features

- **Buoyancy Modeling**: Implements the Boussinesq approximation for efficient simulation of buoyancy-driven flows without solving full compressible equations.
- **Evaporation Tracking**: Includes models for phase change at interfaces, tracking evaporation rates based on temperature, concentration, or saturation conditions.
- **Multiphase Support**: Handles two-phase flows (e.g., liquid-gas) with interface capturing or tracking techniques.
- **Heat Transfer**: Coupled energy equation solving for temperature fields influencing buoyancy and evaporation.
- **Customizable**: Easily extendable for specific boundary conditions, turbulence models (e.g., k-epsilon, LES), or additional physics.
- **Parallel Computing**: Compatible with OpenFOAM's parallel processing for large-scale simulations.
- **Post-Processing**: Integrated with ParaView or OpenFOAM utilities for visualization and analysis.

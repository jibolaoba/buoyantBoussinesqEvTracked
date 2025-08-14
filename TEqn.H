#include "readTransportProperties.H"
{
    alphat = turbulence->nut()/Pr;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff
    (
        IOobject
        (
            "alphaEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        turbulence->nu()/Pr + alphat
    );

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
      + fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve(mesh.solver("T"));

    radiation->correct();

    fvOptions.correct(T);

    // Recalculating kinematic density
    rhok = 1.0 - beta*(T - TRef);
}

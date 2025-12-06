/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "dynamicRefineFvMesh.H"
#include "cellQuality.H"

///
// --- Timing support (v2412 compatible)
#include "cpuTime.H"
#include <sstream>

inline std::string fmt8(double x)
{
    std::ostringstream os;
    os.setf(std::ios::fixed, std::ios::floatfield);
    os.precision(8);
    os << x;
    return os.str();
}

struct ScopedCpu
{
    cpuTime t; double& acc;
    explicit ScopedCpu(double& a) : acc(a) {}
    ~ScopedCpu() { acc += t.elapsedCpuTime(); }
};
///

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using isoAdvector phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing.\n"
        "The solver is derived from interFoam"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    #include "porousCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
///
/*
    	#include "readDyMControls.H"
        #include "porousCourantNo.H"
        #include "porousAlphaCourantNo.H"
        #include "setDeltaT.H"
*/
        // Per-step accumulators and total
        cpuTime tStep;
        double accCtrl=0.0, accReconst=0.0, accMeshUpd=0.0, accRemap=0.0, accMRF=0.0;
        double accCorrPhi=0.0, accAlpha=0.0, accMix=0.0, accUEqn=0.0, accPEqn=0.0;
        double accTurb=0.0, accWrite=0.0;

        { ScopedCpu s(accCtrl);
            #include "readDyMControls.H"
            #include "porousCourantNo.H"
            #include "porousAlphaCourantNo.H"
            #include "setDeltaT.H"
        }
///
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
///
        cpuTime tPIMPLE; // reference (not summed)
///
    	while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                if (isA<dynamicRefineFvMesh>(mesh))
                {
///
//                    advector.surf().reconstruct();
                    {
                    	ScopedCpu s(accReconst);
                        advector.surf().reconstruct();
                    }                	
///                	
                }

///            	
//                mesh.update();
                {
                	ScopedCpu s(accMeshUpd);
                    mesh.update();
                }
///
            	
                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    if (isA<dynamicRefineFvMesh>(mesh))
                    {
///
/*
                    	advector.surf().mapAlphaField();
                        alpha2 = 1.0 - alpha1;
                        alpha2.correctBoundaryConditions();
                        rho == alpha1*rho1 + alpha2*rho2;
                        rho.correctBoundaryConditions();
                        rho.oldTime() = rho;
                        alpha2.oldTime() = alpha2;
*/
                        {
                        	ScopedCpu s(accRemap);
                            advector.surf().mapAlphaField();
                            alpha2 = 1.0 - alpha1;
                            alpha2.correctBoundaryConditions();
                            rho = alpha1*rho1 + alpha2*rho2;   // ← 代入に修正
                            rho.correctBoundaryConditions();
                            rho.oldTime() = rho;
                            alpha2.oldTime() = alpha2;
                        }
///
                    }

///                	
//                    MRF.update();
                    {
                        ScopedCpu s(accMRF);
                        MRF.update();
                    }                	
///
                	
                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

///                    	
//                        #include "correctPhi.H"
                        {
                            ScopedCpu s(accCorrPhi);
                            #include "correctPhi.H"                    	
///
                    	
                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
///
                        }
///
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "alphaControls.H"

///
//        	#include "alphaEqnSubCycle.H"
        	{
        		ScopedCpu s(accAlpha);
                #include "alphaEqnSubCycle.H"
            }
///
        	
///        	
//            mixture.correct();
            {
            	ScopedCpu s(accMix);
                mixture.correct();
            }        	
///
        	
        	
            if (pimple.frozenFlow())
            {
                continue;
            }

///        	
//            #include "UEqn.H"
            cpuTime tUEqn;
            #include "UEqn.H"
            accUEqn += tUEqn.elapsedCpuTime();    	
///
        	
            // --- Pressure corrector loop
            while (pimple.correct())
            {
///
//            	#include "pEqn.H"
                {
                	cpuTime tPEqn;
                    #include "pEqn.H"
                    accPEqn += tPEqn.elapsedCpuTime();
                }            	
///
            }

            if (pimple.turbCorr())
            {
///
//            	turbulence->correct();
                {
                	ScopedCpu s(accTurb);
                    turbulence->correct();
                }            	
///
            }
        }

///    	
//        runTime.write();
        {
        	ScopedCpu s(accWrite);
            runTime.write();
        }    	
///
    	
        runTime.printExecutionTime(Info);

///
        // Report
        Info<< "PIMPLE loop: " << fmt8(tPIMPLE.elapsedCpuTime()) << " s" << nl;
        double accTotal =
            accCtrl + accReconst + accMeshUpd + accRemap + accMRF +
            accCorrPhi + accAlpha + accMix + accUEqn + accPEqn +
            accTurb + accWrite;
    	
        double tTotal = tStep.elapsedCpuTime();
        // % helper (fixed, 2 decimals)
        auto pct = [&](double v)->std::string {
            std::ostringstream os;
            os.setf(std::ios::fixed, std::ios::floatfield);
            os.precision(2);
            os << (tTotal > 0.0 ? (v * 100.0 / tTotal) : 0.0);
            return os.str();
        };
        // Breakdown with %
        Info<< "=== Time-step breakdown ===" << nl
            << "  controls(Courant/deltaT) : " << fmt8(accCtrl)     << " [s]  (" << pct(accCtrl)     << " %)" << nl
            << "  reconstruct(surf)        : " << fmt8(accReconst)  << " [s]  (" << pct(accReconst)  << " %)" << nl
            << "  mesh.update()            : " << fmt8(accMeshUpd)  << " [s]  (" << pct(accMeshUpd)  << " %)" << nl
            << "  remap(mapAlpha etc.)     : " << fmt8(accRemap)    << " [s]  (" << pct(accRemap)    << " %)" << nl
            << "  MRF.update()             : " << fmt8(accMRF)      << " [s]  (" << pct(accMRF)      << " %)" << nl
            << "  correctPhi(+relative)    : " << fmt8(accCorrPhi)  << " [s]  (" << pct(accCorrPhi)  << " %)" << nl
            << "  alphaEqnSubCycle         : " << fmt8(accAlpha)    << " [s]  (" << pct(accAlpha)    << " %)" << nl
            << "  mixture.correct()        : " << fmt8(accMix)      << " [s]  (" << pct(accMix)      << " %)" << nl
            << "  UEqn(assemble+solve)     : " << fmt8(accUEqn)     << " [s]  (" << pct(accUEqn)     << " %)" << nl
            << "  pEqn(correctors total)   : " << fmt8(accPEqn)     << " [s]  (" << pct(accPEqn)     << " %)" << nl
            << "  turbulence.correct()     : " << fmt8(accTurb)     << " [s]  (" << pct(accTurb)     << " %)" << nl
            << "  write()                  : " << fmt8(accWrite)    << " [s]  (" << pct(accWrite)    << " %)" << nl;

    	Info << "Sections total: "  << fmt8(accTotal) << " [s]" << nl;
        Info << "Time-step total: " << fmt8(tTotal)   << " [s]" << nl;
        Info << "Residual: "       << fmt8(tTotal - accTotal) << " [s]" << nl;
///
    
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

// pEqn.H for interIsoFoam_time with GNN output functionality
// PIMPLE algorithm for two-phase incompressible flow

{
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi, Uf))
    );
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        // ===== GNN Output Section =====
        const bool doGnnWrite =
            writePressureSystem
         && pimple.finalNonOrthogonalIter()
         && (runTime.timeIndex() % gnnWriteInterval == 0);

        // Matrix size logging (global and local)
        const label nCellsLocal = p_rghEqn.diag().size();
        const label nUpperLocal = p_rghEqn.upper().size();
        const label nnzLocal    = nCellsLocal + 2*nUpperLocal;

        const label nCellsGlobal = returnReduce(nCellsLocal, sumOp<label>());
        const label nUpperGlobal = returnReduce(nUpperLocal, sumOp<label>());
        const label nnzGlobal    = returnReduce(nnzLocal,   sumOp<label>());

        if (Pstream::master())
        {
            Info<< "p_rgh-matrix global: nCells=" << nCellsGlobal
                << "  nUpper=" << nUpperGlobal
                << "  nnz="    << nnzGlobal
                << nl;
        }

        Pout<< "p_rgh-matrix(proc " << Pstream::myProcNo() << "): "
            << "nCellsMat=" << nCellsLocal
            << "  nUpper="  << nUpperLocal
            << "  nnz="     << nnzLocal
            << nl;

        // ===== GNN data output (A, b, coordinates, etc.) =====
        if (doGnnWrite)
        {
            // 1. Get A (lduMatrix structure)
            const lduMatrix& A = p_rghEqn;

            const scalarField& diag  = A.diag();
            const scalarField& lower = A.lower();
            const scalarField& upper = A.upper();

            const lduAddressing& addr     = A.lduAddr();
            const labelUList&    lowerAddr = addr.lowerAddr();
            const labelUList&    upperAddr = addr.upperAddr();

            // 2. Right-hand side vector b
            const scalarField& b = p_rghEqn.source();

            // 3. Cell center coordinates
            const vectorField& C = mesh.C();

            // 4. Output directory and filename
            fileName outDir(runTime.path()/"gnn");
            mkDir(outDir);
            fileName outName
            (
                outDir
              / (
                    "pEqn_"
                  + runTime.timeName()
                  + "_rank" + Foam::name(Pstream::myProcNo())
                  + ".dat"
                )
            );

            // Mesh size
            const label nCells = mesh.nCells();
            const label nFaces = lowerAddr.size();

            // Boundary mesh (for CoCell and WALL_FACES output)
            const fvBoundaryMesh&  fvb = mesh.boundary();
            const polyBoundaryMesh& pbm = mesh.boundaryMesh();

            // ---------------------------------------------------------
            // 1) Cell quality: skewness / non-orthogonality
            // ---------------------------------------------------------
            cellQuality cq(mesh);

            tmp<scalarField> tCellNonOrtho = cq.nonOrthogonality();
            const scalarField& cellNonOrtho = tCellNonOrtho();

            tmp<scalarField> tCellSkew = cq.skewness();
            const scalarField& cellSkew = tCellSkew();

            // ---------------------------------------------------------
            // 2) Cell aspect ratio (bounding box ratio)
            // ---------------------------------------------------------
            const cellList&  cells = mesh.cells();
            const faceList&  faces = mesh.faces();
            const pointField& pts  = mesh.points();

            scalarField cellAspect(nCells, 1.0);

            for (label celli=0; celli<nCells; ++celli)
            {
                const labelList& cFaces = cells[celli];

                scalar minX = GREAT, maxX = -GREAT;
                scalar minY = GREAT, maxY = -GREAT;
                scalar minZ = GREAT, maxZ = -GREAT;

                forAll(cFaces, fi)
                {
                    const label facei = cFaces[fi];
                    const face& f = faces[facei];

                    forAll(f, pi)
                    {
                        const point& p = pts[f[pi]];
                        minX = Foam::min(minX, p.x());
                        maxX = Foam::max(maxX, p.x());
                        minY = Foam::min(minY, p.y());
                        maxY = Foam::max(maxY, p.y());
                        minZ = Foam::min(minZ, p.z());
                        maxZ = Foam::max(maxZ, p.z());
                    }
                }

                const scalar dx = maxX - minX;
                const scalar dy = maxY - minY;
                const scalar dz = maxZ - minZ;

                const scalar minDim = Foam::max(SMALL, Foam::min(dx, Foam::min(dy, dz)));
                const scalar maxDim = Foam::max(dx, Foam::max(dy, dz));

                cellAspect[celli] = maxDim/minDim;
            }

            // ---------------------------------------------------------
            // 3) diagContrast: max/min ratio of diagonal components among neighbors
            // ---------------------------------------------------------
            scalarField diagMax(nCells, 0.0);
            scalarField diagMin(nCells, GREAT);

            for (label i = 0; i < nCells; ++i)
            {
                const scalar d = Foam::mag(diag[i]);
                diagMax[i] = d;
                diagMin[i] = d;
            }

            for (label f = 0; f < nFaces; ++f)
            {
                const label own = lowerAddr[f];
                const label nei = upperAddr[f];

                const scalar dOwn = Foam::mag(diag[own]);
                const scalar dNei = Foam::mag(diag[nei]);

                diagMax[own] = Foam::max(diagMax[own], dNei);
                diagMin[own] = Foam::min(diagMin[own], dNei);

                diagMax[nei] = Foam::max(diagMax[nei], dOwn);
                diagMin[nei] = Foam::min(diagMin[nei], dOwn);
            }

            scalarField diagContrast(nCells, 1.0);
            for (label i = 0; i < nCells; ++i)
            {
                diagContrast[i] = diagMax[i]/(diagMin[i] + SMALL);
            }

            // ---------------------------------------------------------
            // 4) Cell size and size jump
            // ---------------------------------------------------------
            const scalarField& V = mesh.V();
            scalarField cellSize(nCells, 0.0);

            for (label celli = 0; celli < nCells; ++celli)
            {
                cellSize[celli] = Foam::pow(V[celli], 1.0/3.0);
            }

            scalarField sizeJump(nCells, 1.0);

            for (label f = 0; f < nFaces; ++f)
            {
                const label own = lowerAddr[f];
                const label nei = upperAddr[f];

                const scalar hi = cellSize[own];
                const scalar hj = cellSize[nei];

                scalar ratio = hi/(hj + VSMALL);
                if (ratio < 1.0) ratio = 1.0/ratio;

                sizeJump[own] = Foam::max(sizeJump[own], ratio);
                sizeJump[nei] = Foam::max(sizeJump[nei], ratio);
            }

            // ---------------------------------------------------------
            // 4') Cell Courant number
            // ---------------------------------------------------------
            scalarField CoCell(nCells, 0.0);

            const scalarField& phiIn = phiHbyA.internalField();
            const auto&        phiBf = phiHbyA.boundaryField();

            const label  nInternalFaces = mesh.nInternalFaces();
            const scalar deltaT         = runTime.deltaTValue();

            for (label celli = 0; celli < nCells; ++celli)
            {
                const labelList& cFaces = cells[celli];
                scalar sumMagPhi = 0.0;

                forAll(cFaces, fi)
                {
                    const label facei = cFaces[fi];

                    scalar phiVal = 0.0;
                    if (facei < nInternalFaces)
                    {
                        phiVal = phiIn[facei];
                    }
                    else
                    {
                        const label patchi     = pbm.whichPatch(facei);
                        const label localFacei = facei - pbm[patchi].start();
                        phiVal = phiBf[patchi][localFacei];
                    }

                    sumMagPhi += Foam::mag(phiVal);
                }

                CoCell[celli] = 0.5*sumMagPhi*deltaT/(V[celli] + VSMALL);
            }

            // ---------------------------------------------------------
            // 5) File output
            // ---------------------------------------------------------
            OFstream os(outName);

            os  << "nCells " << nCells << '\n'
                << "nFaces " << nFaces << '\n';

            // Node (cell) information
            os << "CELLS\n";
            for (label i = 0; i < nCells; ++i)
            {
                const vector& ci = C[i];
                os  << i << ' '
                    << ci.x() << ' ' << ci.y() << ' ' << ci.z() << ' '
                    << diag[i]        << ' '
                    << b[i]           << ' '
                    << cellSkew[i]    << ' '
                    << cellNonOrtho[i]<< ' '
                    << cellAspect[i]  << ' '
                    << diagContrast[i]<< ' '
                    << V[i]           << ' '
                    << cellSize[i]    << ' '
                    << sizeJump[i]    << ' '
                    << CoCell[i]
                    << '\n';
            }

            // Edge information
            os << "EDGES\n";
            for (label f = 0; f < nFaces; ++f)
            {
                os  << f << ' '
                    << lowerAddr[f] << ' ' << upperAddr[f] << ' '
                    << lower[f] << ' ' << upper[f] << '\n';
            }

            // ---------------------------------------------------------
            // Wall boundary face information
            // ---------------------------------------------------------
            os << "WALL_FACES\n";

            const vectorField& Cf = mesh.Cf();
            const vectorField& Sf = mesh.Sf();

            forAll(fvb, patchi)
            {
                const fvPatch& fvp = fvb[patchi];

                if (fvp.type() == "wall")
                {
                    const polyPatch& pp        = pbm[patchi];
                    const label      start     = pp.start();
                    const label      size      = pp.size();
                    const labelList& faceCells = pp.faceCells();
                    const word&      patchName = fvp.name();

                    for (label i = 0; i < size; ++i)
                    {
                        const label faceI  = start + i;
                        const label cellI  = faceCells[i];
                        const vector& cf   = Cf[faceI];
                        const vector& sf   = Sf[faceI];

                        os << patchi      << ' '
                           << patchName   << ' '
                           << i           << ' '
                           << faceI       << ' '
                           << cellI       << ' '
                           << cf.x() << ' ' << cf.y() << ' ' << cf.z() << ' '
                           << sf.x() << ' ' << sf.y() << ' ' << sf.z()
                           << '\n';
                    }
                }
            }

            // ---------------------------------------------------------
            // 6) CSR format output for matrix A
            // ---------------------------------------------------------
            {
                const label nRows = nCells;
                const label nnz   = nCells + 2*nFaces;
                fileName csrName(outDir/("A_csr_" + runTime.timeName() + "_rank" + Foam::name(Pstream::myProcNo()) + ".dat"));
                OFstream osCsr(csrName);

                osCsr << "nRows " << nRows << '\n'
                      << "nCols " << nRows << '\n'
                      << "nnz "   << nnz   << '\n';

                List<label> rowNNZ(nRows, 1);
                for (label f = 0; f < nFaces; ++f)
                {
                    const label own = lowerAddr[f];
                    const label nei = upperAddr[f];
                    ++rowNNZ[own];
                    ++rowNNZ[nei];
                }

                List<label> rowPtr(nRows + 1, 0);
                for (label i = 0; i < nRows; ++i)
                {
                    rowPtr[i+1] = rowPtr[i] + rowNNZ[i];
                }

                osCsr << "ROW_PTR\n";
                for (label i = 0; i <= nRows; ++i)
                {
                    osCsr << rowPtr[i] << (i < nRows ? ' ' : '\n');
                }

                List<label>  colInd(nnz);
                List<scalar> val(nnz);

                List<label> rowPos(rowPtr);

                for (label i = 0; i < nRows; ++i)
                {
                    const label idx = rowPos[i]++;
                    colInd[idx] = i;
                    val[idx]    = diag[i];
                }

                for (label f = 0; f < nFaces; ++f)
                {
                    const label own = lowerAddr[f];
                    const label nei = upperAddr[f];

                    {
                        const label idx = rowPos[own]++;
                        colInd[idx] = nei;
                        val[idx]    = upper[f];
                    }

                    {
                        const label idx = rowPos[nei]++;
                        colInd[idx] = own;
                        val[idx]    = lower[f];
                    }
                }

                osCsr << "COL_IND\n";
                for (label k = 0; k < nnz; ++k)
                {
                    osCsr << colInd[k] << (k+1 < nnz ? ' ' : '\n');
                }

                osCsr << "VALUES\n";
                for (label k = 0; k < nnz; ++k)
                {
                    osCsr << val[k] << (k+1 < nnz ? ' ' : '\n');
                }
            }
        } // if (doGnnWrite)

        // Solve
        p_rghEqn.solve();

        // Output solution vector x and true residual r = Ax - b
        if (doGnnWrite)
        {
            const scalarField& x = p_rgh.internalField();
            const scalarField& b = p_rghEqn.source();
            const label nCells = mesh.nCells();

            // Output x
            {
                fileName outDir(runTime.path()/"gnn");
                mkDir(outDir);
                fileName outName
                (
                    outDir
                  / (
                        "x_"
                      + runTime.timeName()
                      + "_rank" + Foam::name(Pstream::myProcNo())
                      + ".dat"
                    )
                );

                OFstream osX(outName);

                for (label i=0; i<nCells; ++i)
                {
                    osX << i << ' ' << x[i] << '\n';
                }
            }

            // Get true residual using fvScalarMatrix::residual()
            tmp<scalarField> tR = p_rghEqn.residual();
            const scalarField& r = tR();

            scalar sumR2  = 0.0;
            scalar sumB2  = 0.0;
            scalar sumAx2 = 0.0;
            scalar maxAbsR = 0.0;

            for (label i=0; i<nCells; ++i)
            {
                const scalar ri = r[i];
                const scalar bi = b[i];
                const scalar Axi = ri + bi;

                sumR2  += ri*ri;
                sumB2  += bi*bi;
                sumAx2 += Axi*Axi;
                maxAbsR = Foam::max(maxAbsR, Foam::mag(ri));
            }

            const scalar normR  = Foam::sqrt(sumR2);
            const scalar normB  = Foam::sqrt(sumB2) + SMALL;
            const scalar normAx = Foam::sqrt(sumAx2) + SMALL;

            Info<< "pEqn true(Ax-b) via residual(): "
                << "||r||_2="        << normR
                << ", max|r_i|="     << maxAbsR
                << ", ||r||/||b||="  << normR/normB
                << ", ||r||/||Ax||=" << normR/normAx
                << nl;

            // Output per-cell residual
            {
                fileName outDir(runTime.path()/"gnn");
                mkDir(outDir);
                fileName rName
                (
                    outDir
                  / (
                        "rTrue_"
                      + runTime.timeName()
                      + "_rank" + Foam::name(Pstream::myProcNo())
                      + ".dat"
                    )
                );

                OFstream osR(rName);

                for (label i=0; i<nCells; ++i)
                {
                    osR << i << ' ' << r[i] << '\n';
                }
            }
        }

        // Final non-orthogonal iteration: correct flux
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            // Output flux divergence
            if (doGnnWrite)
            {
                tmp<volScalarField> tDivPhiHbyA = fvc::div(phiHbyA);
                const scalarField& divPhiHbyA =
                    tDivPhiHbyA().internalField();

                tmp<volScalarField> tDivPhi = fvc::div(phi);
                const scalarField& divPhi =
                    tDivPhi().internalField();

                const scalarField& V = mesh.V();
                const label nCells = mesh.nCells();

                fileName outDir(runTime.path()/"gnn");
                mkDir(outDir);

                fileName divName
                (
                    outDir
                  / (
                        "divPhi_"
                      + runTime.timeName()
                      + "_rank" + Foam::name(Pstream::myProcNo())
                      + ".dat"
                    )
                );

                OFstream osDiv(divName);

                for (label i=0; i<nCells; ++i)
                {
                    osDiv
                        << i << ' '
                        << divPhiHbyA[i] << ' '
                        << divPhi[i]     << ' '
                        << V[i]
                        << '\n';
                }
            }

            p == p_rgh + rho*gh;

            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"
}

/****************************************************************************
 * Copyright (c) 2025 University of Utah
 * Copyright (c) 2025 University of Minnesota Duluth
 *
 * Copyright (c) 2025 Behnam Bozorgmehr
 * Copyright (c) 2025 Jeremy A. Gibbs
 * Copyright (c) 2025 Fabien Margairaz
 * Copyright (c) 2025 Eric R. Pardyjak
 * Copyright (c) 2025 Zachary Patterson
 * Copyright (c) 2025 Rob Stoll
 * Copyright (c) 2025 Lucas Ulmer
 * Copyright (c) 2025 Pete Willemsen
 *
 * This file is part of QES-Winds
 *
 * GPL-3.0 License
 *
 * QES-Winds is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * QES-Winds is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QES-Winds. If not, see <https://www.gnu.org/licenses/>.
 ****************************************************************************/

/**
 * @file Solver.cpp
 * @brief Declares and defines variables required for both solvers.
 *
 * An abstract class that is the basis for the windfield convergence algorithm and
 * contains information needed to run the simulation as well as functions widely
 * used by different solver methods. There are several special member variables that
 * should be accesible to all solvers. They are declared in this class.
 *
 * @sa CPUSolver
 * @sa DynamicParallelism
 * @sa GlobalMemory
 * @sa SharedMemory
 */

#include "Solver.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

// duplication of this macro
#define CELL(i, j, k, w) ((i) + (j) * (nx + (w)) + (k) * (nx + (w)) * (ny + (w)))

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


/**
 * Shows progress of the solving process by printing the percentage.
 */
void Solver::printProgress(float percentage)
{
  int val = (int)(percentage * 100);
  int lpad = (int)(percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
}


/**
 * Assigns values read by WINDSInputData to variables
 * used in the solvers.
 *
 * @note this is only meant to work with QES-Winds!
 *
 * @param WID :document this:
 * @param WGD :document this:
 */
Solver::Solver(qes::Domain domain_in, const float &tolerance)
  : domain(std::move(domain_in)),
    alpha1(1),
    alpha2(1),
    eta(pow((alpha1 / alpha2), 2.0)),
    A(pow((domain.dx() / domain.dy()), 2.0)),
    B(eta * pow((domain.dx() / domain.dz()), 2.0))

{
  tol = tolerance;

  lambda.resize(domain.numCellCentered(), 0.0);
  lambda_old.resize(domain.numCellCentered(), 0.0);
  R.resize(domain.numCellCentered(), 0.0);
}

void Solver::buildRHS(const WINDSGeneralData *w)
{
  if (!std::isfinite(tol) || tol <= 0)
    throw std::runtime_error("Solver tolerance must be positive and finite");
  std::fill(R.begin(), R.end(), 0.0);
  // Only rows on which the Poisson equation is solved. In particular k+1
  // must be inside dz_array; the prescribed top layer has no RHS equation.
  for (int k=1; k<domain.nz()-2; ++k)
    for (int j=1; j<domain.ny()-2; ++j)
      for (int i=1; i<domain.nx()-2; ++i) {
        const auto c=domain.cell(i,j,k), f=domain.face(i,j,k);
        R[c] = -2.0*alpha1*alpha1 * (
          (double(w->e[c])*w->u0[f+1]-double(w->f[c])*w->u0[f])*domain.dx()
          + (double(w->g[c])*w->v0[f+domain.nx()]-double(w->h[c])*w->v0[f])*domain.dy()
          + double(w->m[c])*0.5*(domain.dz_array[k]+domain.dz_array[k+1])
              * w->w0[f+domain.nx()*domain.ny()]
          - double(w->n[c])*0.5*(domain.dz_array[k]+domain.dz_array[k-1])*w->w0[f]);
        if (!std::isfinite(R[c])) throw std::runtime_error("Non-finite solver input");
      }
}

void Solver::maskSolidFaces(WINDSGeneralData *w) const
{
  const auto solid = [&](int i,int j,int k) {
    if (i<0 || j<0 || k<1 || i>=domain.nx()-1 || j>=domain.ny()-1 || k>=domain.nz()-1)
      return false;
    const int flag=w->icellflag[domain.cell(i,j,k)];
    return flag==0 || flag==2;
  };
  // Each face has one writer. Cell-owned writes race with neighbouring
  // velocity corrections (and with each other under OpenMP/CUDA).
  for (int k=1;k<domain.nz();++k)
    for (int j=0;j<domain.ny();++j)
      for (int i=0;i<domain.nx();++i) {
        const auto f=domain.face(i,j,k);
        const bool here=solid(i,j,k);
        if (here || solid(i-1,j,k)) w->u[f]=0;
        if (here || solid(i,j-1,k)) w->v[f]=0;
        if (here || solid(i,j,k-1)) w->w[f]=0;
      }
}

void Solver::recoverVelocity(WINDSGeneralData *w) const
{
  w->u=w->u0; w->v=w->v0; w->w=w->w0;
  // Include the face entering the prescribed top layer: omitting it leaves
  // the last solved layer with an uncorrected outgoing vertical flux.
  for (int k=1;k<domain.nz()-1;++k)
    for (int j=1;j<domain.ny()-1;++j)
      for (int i=1;i<domain.nx()-1;++i) {
        const auto c=domain.cell(i,j,k), f=domain.face(i,j,k);
        w->u[f] += (0.5/(alpha1*alpha1))*w->f[c]*domain.dx()*(lambda[c]-lambda[c-1]);
        w->v[f] += (0.5/(alpha1*alpha1))*w->h[c]*domain.dy()*(lambda[c]-lambda[c-(domain.nx()-1)]);
        // The prescribed top row is not scaled by Wall::solverCoefficients.
        // Its lower face must use the scaled upper-face coefficient below it.
        const double vertical = k==domain.nz()-2
          ? double(w->m[c-(domain.nx()-1)*(domain.ny()-1)])*domain.dz_array[k-1]
          : double(w->n[c])*domain.dz_array[k];
        w->w[f] += (0.5/(alpha2*alpha2))*vertical
                    *(lambda[c]-lambda[c-(domain.nx()-1)*(domain.ny()-1)]);
      }
  maskSolidFaces(w);
}

void Solver::reportConvergence(const WINDSGeneralData *w, int iterations, double maxUpdate, int limit)
{
  diagnostics = {};
  diagnostics.iterations=iterations;
  diagnostics.maxUpdate=maxUpdate;
  double squares=0;
  size_t count=0;
  for (int k=1;k<domain.nz()-2;++k)
    for (int j=1;j<domain.ny()-2;++j)
      for (int i=1;i<domain.nx()-2;++i) {
        const auto c=domain.cell(i,j,k), f=domain.face(i,j,k);
        if (w->icellflag[c]==0 || w->icellflag[c]==2) continue;
        // Difference form avoids cancellation of large nearly equal potentials.
        const double equation = w->e[c]*(lambda[c+1]-lambda[c])
          + w->f[c]*(lambda[c-1]-lambda[c])
          + w->g[c]*(lambda[c+domain.nx()-1]-lambda[c])
          + w->h[c]*(lambda[c-domain.nx()+1]-lambda[c])
          + w->m[c]*(lambda[c+(domain.nx()-1)*(domain.ny()-1)]-lambda[c])
          + w->n[c]*(lambda[c-(domain.nx()-1)*(domain.ny()-1)]-lambda[c])-R[c];
        if (!std::isfinite(equation) || !std::isfinite(maxUpdate))
          throw std::runtime_error("Non-finite wind solver result");
        diagnostics.maxEquationResidual=std::max(diagnostics.maxEquationResidual,std::abs(equation));
        // Direct staggered flux check on regular/stair-step air cells. Cut-cell
        // fluxes require area/volume weighting and are intentionally not claimed here.
        if (w->icellflag[c]==7 || w->icellflag[c]==8) continue;
        const double div=(double(w->u[f+1])-w->u[f])/domain.dx()
          +(double(w->v[f+domain.nx()])-w->v[f])/domain.dy()
          +(double(w->w[f+domain.nx()*domain.ny()])-w->w[f])/domain.dz_array[k];
        if (!std::isfinite(div)) throw std::runtime_error("Non-finite final velocity divergence");
        diagnostics.maxDivergence=std::max(diagnostics.maxDivergence,std::abs(div));
        squares+=div*div; ++count;
        diagnostics.netFlux+=div*domain.dx()*domain.dy()*domain.dz_array[k];
      }
  diagnostics.rmsDivergence=count ? std::sqrt(squares/count) : 0;
  diagnostics.converged=std::isfinite(maxUpdate) && maxUpdate<=tol && diagnostics.maxEquationResidual<=tol;
  printf("[Solver] %s after %d iterations (limit %d): max update=%.9g, equation residual=%.9g, tolerance=%.9g\n",
    diagnostics.converged ? "CONVERGED" : "NOT CONVERGED", iterations,limit,maxUpdate,diagnostics.maxEquationResidual,double(tol));
  printf("[Solver] Solved regular/stair air cells: %zu; divergence max=%.9g s^-1, RMS=%.9g s^-1, net flux=%.9g m^3/s (excludes prescribed outer ring/top and cut cells)\n",
    count,diagnostics.maxDivergence,diagnostics.rmsDivergence,diagnostics.netFlux);
  if (!diagnostics.converged)
    std::cerr << "[WARNING] Wind solve did not meet tolerance; output is NOT a converged solution.\n";
}

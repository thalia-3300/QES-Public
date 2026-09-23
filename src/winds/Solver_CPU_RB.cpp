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

/** @file Solver_CPU_RB.cpp */

#include "Solver_CPU_RB.h"

/** :document this:
 * Start by writing a one sentence description here
 *
 * Complete by continuting to write implementation details here.
 * (remove :document this: tag when done)
 */
void Solver_CPU_RB::solve(WINDSGeneralData *WGD, const int &itermax)
{
  auto startOfSolveMethod = std::chrono::high_resolution_clock::now();// Start recording execution time

  /***************************************************************
   *********   Divergence of the initial velocity field   ********
   ***************************************************************/
  long icell_face;// cell-face index
  long icell_cent;// cell-centered index

  // R.resize(WGD->numcell_cent, 0.0);
  // lambda.resize(WGD->numcell_cent, 0.0);
  // lambda_old.resize(WGD->numcell_cent, 0.0);

  auto startSolveSection = std::chrono::high_resolution_clock::now();
  buildRHS(WGD);
  // INSERT CANOPY CODE

  /***************************************************************
   **********************   SOR Solver   *************************
   ***************************************************************/

  int iter = 0;
  double error;
  double max_error = std::numeric_limits<double>::infinity();
  // int i_max, j_max, k_max;

  std::cout << "[Solver]\t Running Red/Black CPU Solver ..." << std::endl;

  while (iter < itermax && max_error > tol) {
    max_error = 0.0; // Reset once, outside the OpenMP team.

// Save previous iteration values for error calculation
#pragma omp parallel private(icell_cent, icell_face, error) default(none) shared(WGD, lambda, lambda_old, max_error, R)
    {
#pragma omp for
      for (auto k = 0u; k < lambda.size(); ++k) {
        lambda_old[k] = lambda[k];
      }
      // end of omp for (with implicit barrier)


      // main SOR formulation loop

      // Red nodes pass
#pragma omp for
      for (int k = 1; k < domain.nz() - 2; ++k) {
        for (int j = 1; j < domain.ny() - 2; ++j) {
          for (int i = 1; i < domain.nx() - 2; ++i) {

            if (((i + j + k) % 2) == 0) {
              // icell_cent = i + j * (WGD->nx - 1) + k * (WGD->nx - 1) * (WGD->ny - 1);
              icell_cent = domain.cell(i, j, k);
              lambda[icell_cent] = (omega / (double(WGD->e[icell_cent]) + WGD->f[icell_cent] + WGD->g[icell_cent] + WGD->h[icell_cent] + WGD->m[icell_cent] + WGD->n[icell_cent]))
                                     * (WGD->e[icell_cent] * lambda[domain.cellAdd(icell_cent, 1, 0, 0)]
                                        + WGD->f[icell_cent] * lambda[domain.cellAdd(icell_cent, -1, 0, 0)]
                                        + WGD->g[icell_cent] * lambda[domain.cellAdd(icell_cent, 0, 1, 0)]
                                        + WGD->h[icell_cent] * lambda[domain.cellAdd(icell_cent, 0, -1, 0)]
                                        + WGD->m[icell_cent] * lambda[domain.cellAdd(icell_cent, 0, 0, 1)]
                                        + WGD->n[icell_cent] * lambda[domain.cellAdd(icell_cent, 0, 0, -1)] - R[icell_cent])
                                   + (1.0f - omega) * lambda[icell_cent];// SOR formulation
            }
          }
        }
      }
      // end of omp for (with implicit barrier)

      // Black nodes pass
#pragma omp for
      for (int k = 1; k < domain.nz() - 2; ++k) {
        for (int j = 1; j < domain.ny() - 2; ++j) {
          for (int i = 1; i < domain.nx() - 2; ++i) {
            if (((i + j + k) % 2) == 1) {
              // icell_cent = i + j * (WGD->nx - 1) + k * (WGD->nx - 1) * (WGD->ny - 1);
              icell_cent = domain.cell(i, j, k);
              lambda[icell_cent] = (omega / (double(WGD->e[icell_cent]) + WGD->f[icell_cent] + WGD->g[icell_cent] + WGD->h[icell_cent] + WGD->m[icell_cent] + WGD->n[icell_cent]))
                                     * (WGD->e[icell_cent] * lambda[domain.cellAdd(icell_cent, 1, 0, 0)]
                                        + WGD->f[icell_cent] * lambda[domain.cellAdd(icell_cent, -1, 0, 0)]
                                        + WGD->g[icell_cent] * lambda[domain.cellAdd(icell_cent, 0, 1, 0)]
                                        + WGD->h[icell_cent] * lambda[domain.cellAdd(icell_cent, 0, -1, 0)]
                                        + WGD->m[icell_cent] * lambda[domain.cellAdd(icell_cent, 0, 0, 1)]
                                        + WGD->n[icell_cent] * lambda[domain.cellAdd(icell_cent, 0, 0, -1)] - R[icell_cent])
                                   + (1.0f - omega) * lambda[icell_cent];// SOR formulation
            }
          }
        }
      }
      // end of omp for (with implicit barrier)

      // Mirror boundary condition (lambda (@k=0) = lambda (@k=1))
#pragma omp for
      for (int j = 0; j < domain.ny() - 1; ++j) {
        for (int i = 0; i < domain.nx() - 1; ++i) {
          // icell_cent = i + j * (WGD->nx - 1);// Lineralized index for cell centered values
          icell_cent = domain.cell(i, j, 0);
          lambda[icell_cent] = lambda[domain.cellAdd(icell_cent, 0, 0, 1)];
        }
      }
      // end of omp for (with implicit barrier)

      // Error calculation
#pragma omp for reduction(max \
                          : max_error)
      for (int k = 1; k < domain.nz() - 1; ++k) {
        for (int j = 0; j < domain.ny() - 1; ++j) {
          for (int i = 0; i < domain.nx() - 1; ++i) {
            // icell_cent = i + j * (WGD->nx - 1) + k * (WGD->nx - 1) * (WGD->ny - 1);// Lineralized index for cell centered values
            icell_cent = domain.cell(i, j, k);
            error = fabs(lambda[icell_cent] - lambda_old[icell_cent]);
            if (error > max_error) {
              max_error = error;
            }
          }
        }
      }
      // end of omp for (with implicit barrier)
    }
    // end of omp parallel workshare
    iter += 1;
  }

  // std::cout << "Solved!\n";

  // std::cout << "Number of iterations:" << iter << "\n";// Print the number of iterations
  // std::cout << "Error:" << max_error << "\n";
  // std::cout << "tol:" << tol << "\n";


  recoverVelocity(WGD);
  reportConvergence(WGD, iter, max_error, itermax);

  auto finish = std::chrono::high_resolution_clock::now();// Finish recording execution time
  std::chrono::duration<float> elapsedTotal = finish - startOfSolveMethod;
  std::chrono::duration<float> elapsedSolve = finish - startSolveSection;
  std::cout << "\t\t elapsed time: " << elapsedTotal.count() << " s\n";// Print out elapsed execution time
  // std::cout << "Elapsed solve time: " << elapsedSolve.count() << " s\n";// Print out elapsed execution time
}

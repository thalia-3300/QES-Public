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
 * @file Solver_GPU_SharedMemory.cu
 * @brief Child class of the Solver that runs the convergence
 * algorithm using DynamicParallelism on a single GPU.
 *
 * @sa Solver
 * @sa Solver_GPU_SharedMemory
 */

#include "Solver_GPU_SharedMemory.h"

using namespace std::chrono;
using namespace std;
using std::ofstream;
using std::ifstream;
using std::istringstream;
using std::string;
using std::cerr;
using std::endl;
using std::vector;
using std::cout;
using std::to_string;

#define BLOCKSIZE 1024
#define cudaCheck(x) _cudaCheck(x, #x, __FILE__, __LINE__)


template<typename T>
void Solver_GPU_SharedMemory::_cudaCheck(T e, const char *func, const char *call, const int line)
{
  if (e != cudaSuccess) {
    printf("\"%s\" at %d in %s\n\treturned %d\n-> %s\n", func, line, call, (int)e, cudaGetErrorString(e));
    exit(EXIT_FAILURE);
  }
}

#include "Solver_GPU_Kernels.cuh"

void Solver_GPU_SharedMemory::solve(WINDSGeneralData *WGD, const int &itermax)
{
  const auto start=std::chrono::high_resolution_clock::now();
  std::cout << "[Solver] Running shared-reduction double-precision CUDA SOR\n";
  buildRHS(WGD);
  const auto result=qesgpu::solve(domain,WGD,R,lambda,omega,tol,itermax,true);
  // A single shared reconstruction avoids CPU/GPU drift and concurrent writes
  // from solid-cell masks. The iterative Poisson solve remains on the GPU.
  recoverVelocity(WGD);
  reportConvergence(WGD,result.first,result.second,itermax);
  std::cout << "\t elapsed time: "
            << std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-start).count() << " s\n";
}

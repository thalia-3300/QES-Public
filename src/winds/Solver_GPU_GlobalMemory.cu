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
 ***************************************************************************/

/**
 * @file Solver_GPU_GlobalMemory.cu
 * @brief Child class of the Solver that runs the convergence
 * algorithm using Global Memory on a single GPU.
 *
 * @sa Solver
 * @sa Solver_GPU_GlobalMemory
 */

#include "Solver_GPU_GlobalMemory.h"

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
void Solver_GPU_GlobalMemory::_cudaCheck(T e, const char *func, const char *call, const int line)
{
  if (e != cudaSuccess) {
    printf("\"%s\" at %d in %s\n\treturned %d\n-> %s\n", func, line, call, (int)e, cudaGetErrorString(e));
    exit(EXIT_FAILURE);
  }
}

#include "Solver_GPU_Kernels.cuh"

Solver_GPU_GlobalMemory::Solver_GPU_GlobalMemory(qes::Domain domain_in, const float &tolerance)
  : Solver(std::move(domain_in), tolerance)
{
  std::cout << "-------------------------------------------------------------------" << std::endl;
  std::cout << "[Solver]\t Initializing Global Memory Solver (GPU) ..." << std::endl;
  int deviceCount = 0;
  cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

  if (error_id != cudaSuccess) {
    std::cerr << "\n===================================================================" << std::endl;
    std::cerr << "[ERROR]\t cudaGetDeviceCount returned "
              << static_cast<int>(error_id) << "\n\t-> "
              << cudaGetErrorString(error_id) << std::endl;
    std::cerr << "===================================================================" << std::endl;
    exit(EXIT_FAILURE);
  }

  // This function call returnsq 0 if there are no CUDA capable devices.
  if (deviceCount == 0) {
    std::cerr << "\n===================================================================" << std::endl;
    std::cerr << "[ERROR]\t There are no available device(s) that support CUDA\n";
    std::cerr << "===================================================================" << std::endl;
    exit(EXIT_FAILURE);
  }

  cudaSetDevice(0);
}

void Solver_GPU_GlobalMemory::solve(WINDSGeneralData *WGD, const int &itermax)
{
  const auto start=std::chrono::high_resolution_clock::now();
  std::cout << "[Solver] Running global-reduction double-precision CUDA SOR\n";
  buildRHS(WGD);
  const auto result=qesgpu::solve(domain,WGD,R,lambda,omega,tol,itermax,false);
  // A single shared reconstruction avoids CPU/GPU drift and concurrent writes
  // from solid-cell masks. The iterative Poisson solve remains on the GPU.
  recoverVelocity(WGD);
  reportConvergence(WGD,result.first,result.second,itermax);
  std::cout << "\t elapsed time: "
            << std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-start).count() << " s\n";
}

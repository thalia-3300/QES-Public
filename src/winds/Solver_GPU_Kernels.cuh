// Shared, checked CUDA implementation for the global/shared reduction backends.
#pragma once
#include <cuda_runtime.h>
#include <math_constants.h>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qesgpu {
inline void check(cudaError_t result)
{
  if (result != cudaSuccess)
    throw std::runtime_error(std::string("CUDA wind solver: ")+cudaGetErrorString(result));
}

template<class T> class Buffer {
public:
  T *data=nullptr;
  explicit Buffer(size_t size) { check(cudaMalloc(reinterpret_cast<void **>(&data),size*sizeof(T))); }
  ~Buffer() { if (data) cudaFree(data); }
  Buffer(const Buffer &)=delete;
  Buffer &operator=(const Buffer &)=delete;
  void upload(const std::vector<T> &v) { check(cudaMemcpy(data,v.data(),v.size()*sizeof(T),cudaMemcpyHostToDevice)); }
};

struct Coefficients { const float *e,*f,*g,*h,*m,*n; };
constexpr int blockSize=256;

static __global__ void relax(double *potential, const double *rhs, Coefficients a,
                      int nx,int ny,int nz,double omega,int color)
{
  const size_t id=size_t(blockIdx.x)*blockDim.x+threadIdx.x;
  const size_t plane=size_t(nx-1)*(ny-1);
  const size_t k=id/plane, j=(id%plane)/(nx-1), i=id%(nx-1);
  if (i>0 && i<size_t(nx-2) && j>0 && j<size_t(ny-2) && k>0 && k<size_t(nz-2)
      && (i+j+k)%2==color) {
    const double diagonal=double(a.e[id])+a.f[id]+a.g[id]+a.h[id]+a.m[id]+a.n[id];
    if (diagonal>0) {
      const double sum=a.e[id]*potential[id+1]+a.f[id]*potential[id-1]
        +a.g[id]*potential[id+nx-1]+a.h[id]*potential[id-nx+1]
        +a.m[id]*potential[id+plane]+a.n[id]*potential[id-plane]-rhs[id];
      potential[id]=(omega/diagonal)*sum+(1.0-omega)*potential[id];
    }
  }
}

static __global__ void mirrorGround(double *potential,size_t plane)
{
  const size_t id=size_t(blockIdx.x)*blockDim.x+threadIdx.x;
  // Cell-centred plane, NOT nx*ny (which overwrites the first physical layer).
  if (id<plane) potential[id]=potential[id+plane];
}

__device__ inline void atomicMaxPositive(double *address,double value)
{
  auto *bits=reinterpret_cast<unsigned long long *>(address);
  unsigned long long old=atomicCAS(bits,0ULL,0ULL);
  while (value>__longlong_as_double(old)) {
    const auto assumed=old;
    old=atomicCAS(bits,assumed,__double_as_longlong(value));
    if (old==assumed) break;
  }
}

template<bool sharedReduction>
__global__ void maxUpdate(const double *potential,const double *previous,size_t count,double *result)
{
  const size_t id=size_t(blockIdx.x)*blockDim.x+threadIdx.x;
  double value=id<count ? fabs(potential[id]-previous[id]) : 0;
  if (!isfinite(value)) value=CUDART_INF;
  if constexpr (sharedReduction) {
    __shared__ double values[blockSize];
    values[threadIdx.x]=value;
    __syncthreads();
    for (int stride=blockSize/2;stride;stride/=2) {
      if (threadIdx.x<stride) values[threadIdx.x]=fmax(values[threadIdx.x],values[threadIdx.x+stride]);
      __syncthreads();
    }
    // No cross-block barrier assumption: one atomic maximum per block.
    if (threadIdx.x==0) atomicMaxPositive(result,values[0]);
  } else {
    atomicMaxPositive(result,value);
  }
}

inline std::pair<int,double> solve(const qes::Domain &d,const WINDSGeneralData *w,
                                  const std::vector<double> &rhs,std::vector<double> &potential,
                                  double omega,double tolerance,int limit,bool sharedReduction)
{
  const size_t count=potential.size(), plane=size_t(d.nx()-1)*(d.ny()-1);
  Buffer<double> current(count),previous(count),r(count),error(1);
  Buffer<float> e(count),f(count),g(count),h(count),m(count),n(count);
  current.upload(potential); r.upload(rhs);
  e.upload(w->e); f.upload(w->f); g.upload(w->g); h.upload(w->h); m.upload(w->m); n.upload(w->n);
  const Coefficients coefficients{e.data,f.data,g.data,h.data,m.data,n.data};
  const unsigned blocks=(count+blockSize-1)/blockSize;
  int iterations=0;
  double update=std::numeric_limits<double>::infinity();
  while (iterations<limit && update>tolerance) {
    check(cudaMemcpy(previous.data,current.data,count*sizeof(double),cudaMemcpyDeviceToDevice));
    relax<<<blocks,blockSize>>>(current.data,r.data,coefficients,d.nx(),d.ny(),d.nz(),omega,0);
    check(cudaGetLastError());
    relax<<<blocks,blockSize>>>(current.data,r.data,coefficients,d.nx(),d.ny(),d.nz(),omega,1);
    check(cudaGetLastError());
    mirrorGround<<<(plane+blockSize-1)/blockSize,blockSize>>>(current.data,plane);
    check(cudaGetLastError());
    check(cudaMemset(error.data,0,sizeof(double)));
    if (sharedReduction) maxUpdate<true><<<blocks,blockSize>>>(current.data,previous.data,count,error.data);
    else maxUpdate<false><<<blocks,blockSize>>>(current.data,previous.data,count,error.data);
    check(cudaGetLastError());
    check(cudaMemcpy(&update,error.data,sizeof(double),cudaMemcpyDeviceToHost));
    if (!std::isfinite(update)) throw std::runtime_error("Non-finite CUDA solver update");
    ++iterations;
  }
  check(cudaMemcpy(potential.data(),current.data,count*sizeof(double),cudaMemcpyDeviceToHost));
  return {iterations,update};
}
} // namespace qesgpu

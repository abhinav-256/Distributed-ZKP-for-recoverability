#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include <string.h>
#include <pbc.h>
#include <time.h>

#include <cuda_runtime_api.h>
#include <cuda.h>

using namespace std;

pairing_t ctx;
element_t g1, g2, eg1g2;
element_pp_t ppg1, ppg2;

void bbkeygen(element_t _sk, element_t pk, pairing_t ctx) {
    element_random(_sk);
    element_pow_zn(pk, g2, _sk);
}

void bbsign(element_t sig, element_t m, element_t _sk, pairing_t ctx) {
    element_t tmp1, tmp2;
    element_init_Zr(tmp1, ctx);
    element_init_Zr(tmp2, ctx);

    element_add(tmp1, m, _sk);
    element_invert(tmp2, tmp1);
    element_pp_pow_zn(sig, tmp2, ppg1);
}

int bbverify(element_t sig, element_t m, element_t pk, pairing_t ctx) {
    element_t tmp1, tmp2, tmp3;
    element_init_G2(tmp1, ctx);
    element_init_G2(tmp2, ctx);
    element_init_GT(tmp3, ctx);

    element_pp_pow_zn(tmp1, m, ppg2);
    element_mul(tmp2, pk, tmp1);
    element_pairing(tmp3, sig, tmp2);

    return element_cmp(tmp3, eg1g2);
}

void initset(element_t* phi, int n) {
  for (int i=0; i<n; i++) {
    element_random(phi[i]);
  }
}

__global__ void bbsign_kernel(element_t* phi, element_t* sigs, element_t _sk, int n, pairing_t ctx) {
  int i = blockDim.x + blockIdx.x + threadIdx.x;
  if (i<n) {
    element_t tmp1, tmp2;
    element_init_Zr(tmp1, ctx);
    element_init_Zr(tmp2, ctx);

    element_add(tmp1, phi[i], _sk);
    element_invert(tmp2, tmp1);
    element_pp_pow_zn(sigs[i], tmp2, ppg1);
  }
}

void verfsigs(element_t pk, element_t* sigs, element_t* phi, int n, pairing_t ctx) {
  // declare private and secret keys
  element_t _sk;
  element_init_Zr(_sk, ctx);
  bbkeygen(_sk, pk, ctx);

  element_t* d_sigs;
  element_t* d_phi;

  cudaMalloc(&d_sigs, n * sizeof(element_t));
  cudaMalloc(&d_phi, n * sizeof(element_t));

  cudaMemcpy(d_phi, phi, n * sizeof(element_t), cudaMemcpyHostToDevice);
  int threadsPerBlock = 256;
  int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;

  bbsign_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_phi, d_sigs, _sk, n, ctx);

  cudaMemcpy(d_sigs, sigs, n * sizeof(element_t), cudaMemcpyDeviceToHost);

  cudaFree(d_sigs);
  cudaFree(d_phi);

  // Sequential code
  // for (int i=0; i<n; i++) {
  //   bbsign(sigs[i], phi[i], _sk, ctx);
  // }
}


int main(int argc, char** argv) {
    // Generic time markers
    clock_t begin, end;

    // initialise pairing context from the curve parameters file
    char param[1024];
    char* paramfilename = argv[1];
    FILE *fp = fopen(paramfilename, "r");
    size_t count = fread(param, 1, 1024, fp);
    if (!count) pbc_die("input error");
    pairing_init_set_buf(ctx, param, count);

    // initialise generators and their pre-computations
    element_init_G1(g1, ctx);
    element_init_G2(g2, ctx);
    element_init_GT(eg1g2, ctx);
    element_random(g1);
    element_random(g2);
    element_pairing(eg1g2, g1, g2);
    element_pp_init(ppg1, g1);
    element_pp_init(ppg2, g2);

    // declare set
    int n = 1000;
    element_t* phi = (element_t*) malloc(n * sizeof(element_t));
    for (int i=0; i<n; i++) {
        element_init_Zr(phi[i], ctx);
    }
    initset(phi, n);

    // obtain verifier signatures on all elements of the set
    begin = clock();
    element_t pk;
    element_init_G2(pk, ctx);
    element_t* sigs = (element_t*) malloc(n * sizeof(element_t));
    for (int i=0; i<n; i++) {
        element_init_G1(sigs[i], ctx);
    }
    verfsigs(pk, sigs, phi, n, ctx);
    end = clock();
    printf("Generating %d verifier signatures: %lf sec\n", n, (double)(end-begin)/CLOCKS_PER_SEC);

    free(phi);
    free(sigs);
    element_pp_clear(ppg1);

    return 0;
}

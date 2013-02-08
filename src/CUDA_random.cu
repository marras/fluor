/* CUDA-based free diffusion simulator
   To be used with "fluor" program
   (c) Marek Waligórski 2009-2011       */

//Include some standard headers (mixed C++ and CUDA code)
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cutil.h>
//The CUDA functions are called from a class method
#include "fluorescence.h"

#define PI 3.1415926535897932384626

/* Random number generator
   source: http://http.developer.nvidia.com/GPUGems3/gpugems3_ch37.html
   
   For molecule movement, we need Gaussian-distributed numbers,
   so we take a Combined Tausworthe generator (linear),
   plus a Box-Muller transform into Gaussian distribution.
   
   NOTE: NOT PERFECT! */

/* Tausworthe generator step: S1, S2, S3, and M are all constants, 
   z is part of the private per-thread generator state (~seed).  */
__device__ unsigned TausStep(unsigned &z, int S1, int S2, int S3, unsigned M) {  
    unsigned b=(((z << S1) ^ z) >> S2);
    return z = (((z & M) << S3) ^ b);  
}  

/* Linear Congruential Generator step */
__device__ unsigned LCGStep(unsigned &z, unsigned A, unsigned C)  {  
    return z=(A*z+C);
}

/* Linear generator consists of XORing 3 Tausworthe-generated numbers with one LCG.*/
__device__ float HybridTaus(unsigned *z) {  
   // Combined period is lcm(p1,p2,p3,p4)~ 2^121  
    return 2.3283064365387e-10 * (              // Periods  
     TausStep(z[0], 13, 19, 12, 4294967294UL) ^  // p1=2^31-1  
     TausStep(z[1], 2, 25, 4, 4294967288UL) ^    // p2=2^30-1  
     TausStep(z[2], 3, 11, 17, 4294967280UL) ^   // p3=2^28-1  
     LCGStep(z[3], 1664525, 1013904223UL)        // p4=2^32  
    );  
}  


/* Box-Muller transform: get 2 normally-distributed variables from 2 linear vars. */
__device__ float2 BoxMuller(unsigned *z) {  
   float u0=HybridTaus (z), u1=HybridTaus (z);  
   float r=sqrt(-2*log(u0));	//radius
   float theta=2*PI*u1;  	//angle
   return make_float2(r*sin(theta),r*cos(theta));  //here we get the actual (x,y) coordinates
}  

/* Simulate diffusion of molecules on GPU:
	da - positions of molecules
	ds - randomizer seed
	natoms - how many molecules to movement
	cudaDIFF_STEP - avg. movement distance
	size - simulation box size
*/
__global__ void DoTheMovement (float4 *da, unsigned *ds, int natoms, float cudaDIFF_STEP, float3 size) {
    unsigned z[4]={0,0,0,0}; //What the hell? Gdyby zadeklarowac ta zmienna "na zewnatrz" funkcji, to pomimo __device__ i tak bylaby tez widoczna z poziomu hosta!!!

    int moj_nr = blockIdx.x * blockDim.x + threadIdx.x; 	//this number is different for every thread
    z[0]=ds[moj_nr*4];z[1]=ds[1+moj_nr*4];z[2]=ds[2+moj_nr*4];z[3]=ds[3+moj_nr*4]; //faster than a 'for' loop
    
    float2 para = BoxMuller(z);  //get a pair of normally-distributed values

    //move molecules, enable periodic boundary conditions
    da[moj_nr].x += cudaDIFF_STEP * para.x; if (da[moj_nr].x > size.x) da[moj_nr].x -= size.x; else if (da[moj_nr].x < 0) da[moj_nr].x += size.x;
    da[moj_nr].y += cudaDIFF_STEP * para.y; if (da[moj_nr].y > size.y) da[moj_nr].y -= size.y; else if (da[moj_nr].y < 0) da[moj_nr].y += size.y;
    
    para = BoxMuller(z);  //generate next pair of numbers (one of them will get wated :-( )
    da[moj_nr].z += cudaDIFF_STEP * para.x; if (da[moj_nr].z > size.z) da[moj_nr].z -= size.z; else if (da[moj_nr].z < 0) da[moj_nr].z += size.z;
 
    ds[moj_nr*4]=z[0];ds[1+moj_nr*4]=z[1];ds[2+moj_nr*4]=z[2];ds[3+moj_nr*4]=z[3];	 //we have to save the random generator state
}

/* Check excitation of molecules:
	da - positions of molecules
	ds - randomizer seed
	natoms - how many molecules to movement
	dstates - array containing states of each molecule
	epsRAZYdT - probability of excitation at center of confocal volume during this timestep
	F1 - focal point
	SQR_KAPPA - structure parameter of confocal volume (squared)
	WXY - confocal volume radius (XY plane)
*/
__global__ void ExciteMeBaby (float4 *da, unsigned *ds, int natoms, enumStates *dstates, float epsRAZYdT, float3 F1, float SQR_KAPPA, float WXY) { //NOTE we'll try to use the da array already loaded with molecule positions and ds with seeds :)
    int moj_nr = blockIdx.x * blockDim.x + threadIdx.x;
    float4 pos = da[moj_nr];
    
    unsigned z[4]={0,0,0,0}; 
    z[0]=ds[moj_nr*4];z[1]=ds[1+moj_nr*4];z[2]=ds[2+moj_nr*4];z[3]=ds[3+moj_nr*4];
    
    float prob_ex = epsRAZYdT * exp(-2*((pos.x-F1.x)*(pos.x-F1.x)	//probability of excitation (Gaussian profile) [Winkler]
				+(pos.y-F1.y)*(pos.y-F1.y)
				+(pos.z-F1.z)*(pos.z-F1.z)/SQR_KAPPA ) / (WXY*WXY));
    float r1 = HybridTaus(z); 	//here we need a random (0,1)
    if (r1 < prob_ex && dstates[moj_nr] == MS_GROUND) dstates[moj_nr] = MS_EXC_1; //NOTE only excite if molecule is in ground state!

    ds[moj_nr*4]=z[0];ds[1+moj_nr*4]=z[1];ds[2+moj_nr*4]=z[2];ds[3+moj_nr*4]=z[3];	 //save seeds. Optimize somehow?
}


/************************************************************************/
/* CUDA Randomise                                                       */
/************************************************************************/

unsigned *ds; //zarodek randomizera (na karcie graficznej)
float4 *pos, *dpos; //wyniki - pozycje cząstek zmienne z rozkl.norm
enumStates *dstates; //stany czasteczek - na karcie

void Fluorescence :: MoveMolecules_GPU () {
    if (types != 1) {LOG ("!GPU diffusion is not implemented for more than one kind of dye."); return;}
    const float cudaDIFF_STEP = mol[0].DIFFUSION_STEP (); //NOTE THIS ASSUMES ALL MOLECULES HAVE THE SAME SPEED! TODO
    const float3 cudaBox_size = make_float3 (SIZE[0],SIZE[1],SIZE[2]);
    
    DoTheMovement<<<NBlocks,NThreads>>> (dpos, ds, natoms, cudaDIFF_STEP, cudaBox_size); //call the CUDA diffusion kernel
    CUDA_SAFE_CALL( cudaThreadSynchronize() ); 		//wait for GPU to finish calculations
    CUT_CHECK_ERROR("Kernel execution failed\n");
    CUDA_SAFE_CALL( cudaMemcpy(pos, dpos, sizeof(float4) * natoms, cudaMemcpyDeviceToHost));  //copy data from GPU to "pos" array on CPU

    for (int m=0; m<natoms; m++) {
	mol[m].x[0] = pos[m].x; //update info in Molecule class (TODO optimize)
	mol[m].x[1] = pos[m].y;
	mol[m].x[2] = pos[m].z;
    }
}

void Fluorescence:: TryExcitationALL_GPU () {
  assert (profil == PR_GAUSS && num_foci == 1); // The GPU version is only implemented for 1-focus setup with Gauss excitation.

  const float epsRAZYdT = EPSILON * dT;
  const float3 dF1 = make_float3 (F1[0],F1[1],F1[2]);
  
  CUDA_SAFE_CALL( cudaMemcpy(dstates, state, sizeof(enumStates) * natoms, cudaMemcpyHostToDevice));  
  
  ExciteMeBaby<<<NBlocks,NThreads>>> (dpos, ds, natoms, dstates, epsRAZYdT, dF1, SQR_KAPPA, WXY); //check which molecules get excited
  CUDA_SAFE_CALL( cudaThreadSynchronize() );
  CUT_CHECK_ERROR("Kernel execution failed\n");
  
  CUDA_SAFE_CALL( cudaMemcpy(state, dstates, sizeof(enumStates) * natoms, cudaMemcpyDeviceToHost)); //read back mol. states from GPU
}

void seed_rng (unsigned *s, int howmany) {	//s- pointer to seed array
    FILE* urandom = fopen( "/dev/urandom", "rb" );
    if (urandom) {
	    register int i = howmany;
	    register bool success = true;
	    while( success && i-- )
		    success = fread( s++, sizeof(unsigned), 1, urandom );
	    fclose(urandom);
    }
    else printf ("ERROR: Failed to open /dev/urandom!\n"); 
}

#ifdef ENABLE_GPU

/* Initialize CUDA */
void Fluorescence :: InitCUDA () {
    int count = 0;
    int i = 0;

    cudaGetDeviceCount(&count);		//how many devices are available?
    if(count == 0) {
	    fprintf(stderr, "There is no device.\n");
	    LOG ("!Error initializing CUDA.");
    }

    for(i = 0; i < count; i++) {
	    cudaDeviceProp prop;
	    if(cudaGetDeviceProperties(&prop, i) == cudaSuccess) {
		    if(prop.major >= 1) break;
	    }
    }
    if(i == count) {
	    fprintf(stderr, "There is no device supporting CUDA.\n");
	    LOG ("!Error initializing CUDA.");
    }
    cudaSetDevice(i);

    printf("CUDA initialized.\n");

    assert (NBlocks*NThreads == natoms); //1 molecule per thread
    
    seed = new unsigned [4*NBlocks*NThreads];
    seed_rng (seed, 4*NBlocks*NThreads); // Generate Threads*4 random numbers (externally, e.g. /dev/urandom) 

    pos = new float4 [natoms]; //host pointer to array of molecule positions
    CUDA_SAFE_CALL (cudaMalloc ((void**)&dpos, natoms*sizeof(float4))); //allocate mem for molecule positions (device)
    CUDA_SAFE_CALL (cudaMalloc ((void**)&ds, 4*NBlocks*NThreads*sizeof(unsigned))); //and for seeds array (device)

    CUDA_SAFE_CALL (cudaMalloc ((void**)&dstates, natoms*sizeof(enumStates))); //zarezerwuj miejsce na stany atomow (device)
    
    for (int m=0; m<natoms; m++) {
	pos[m] = make_float4 (mol[m].x[0],mol[m].x[1],mol[m].x[2],0); //wrap mol. positions into a convenient structure
    } //(float4 is supposed to work faster than float3)
    
    CUDA_SAFE_CALL( cudaMemcpy(ds, seed, sizeof(unsigned) * 4*NBlocks*NThreads, cudaMemcpyHostToDevice)); //copy seeds to GPU
    CUDA_SAFE_CALL( cudaMemcpy(dpos, pos, sizeof(float4) * natoms, cudaMemcpyHostToDevice)); //copy starting pos. of molecules
   
    LOG ("---- CUDA ----\nBlocks: %d, Threads: %d.",NBlocks,NThreads);
}

/* Clean up */
void Fluorescence :: DeInitCUDA () {
     CUDA_SAFE_CALL( cudaFree(dpos));   
     CUDA_SAFE_CALL( cudaFree(ds));   
     delete [] pos;
     delete [] seed;
}


#else //!ENABLE_GPU

void Fluorescence :: InitCUDA () {
	LOG ("*Using CPU only, no CUDA init required.");
}

void Fluorescence :: DeInitCUDA () {
	LOG ("*Using CPU only, no CUDA deinitialization required.");
}
#endif //ENABLE_GPU

/*
Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

This file is part of the Open Porous Media project (OPM).
OPM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
OPM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <cuda_profiler_api.h>
#include <algorithm>
#include <cmath>

using namespace std;

#define NBLOCKS 1024
#define NTHREADS 256

__global__ void dvc_FiniteVolumeMethod_Init(bool* boolfield, double affinity, double * Phi, double *DenA, double * DenB, int N,  double * VFmask, int Nx, int Ny, int Nz) {
double max = 0;

size_t finish = N;
size_t start = 0;

size_t n;
double denA,denB;

int nn,m;
int S = N/NBLOCKS/NTHREADS + 1;
for (int s=0; s<S; s++) {
//........Get 1-D index for this thread....................
n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
if (n < finish) {
if (boolfield[n] == true) {

if (VFmask[n] > 0 && VFmask[n] < 1) {
max = 0;
if (DenA[n+1] > DenA[n-1]) max = DenA[n+1];
else max = DenA[n-1];
// max = DenA[n+1], DenA[n-1]);
if (DenA[n+Nx] > max) max = DenA[n+Nx];
else max = max;
//max = std::max(DenA[n+Nx],max);
if (DenA[n-Nx] > max) max = DenA[n-Nx];
else max = max;
//max = std::max(max,DenA[n-Nx]);
if (DenA[n+Nx*Ny] > max) max = DenA[n+Nx*Ny];
else max = max;
//max = std::max(DenA[n+Nx*Ny],max);
if (DenA[n-Nx*Ny] > max) max = DenA[n-Nx*Ny];
else max = max;
//max = std::max(DenA[n-Nx*Ny],max);
DenA[n] = max;
if (max == 1) {
DenB[n] = 0;
}
else {
DenB[n] = 1;
}
Phi[n] = (DenA[n] - DenB[n]) / (DenA[n] + DenB[n]) * (1-VFmask[n]) + affinity * VFmask[n];
}
}
}
}

}

__global__ void dvc_MassAportionAndPhaseUpdate(bool * boolfield, char * id, double affinity, double * Phi, double *DenA, double * DenB, int N, double * VFmask, int Nx, int Ny, int Nz) {



double max = 0;
double nnvol = 0;
int count;

double vfs;

/* Determine the neighbor fluid phase fraction based upon neighbors */
int strideY = Nx;
int strideZ = Nx*Ny;

size_t finish = N;
size_t start = 0;

size_t n;
double denA,denB;

int nn,ijk,m;
int S = N/NBLOCKS/NTHREADS + 1;
for (int s=0; s<S; s++) {
//........Get 1-D index for this thread....................
ijk =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
if (ijk < finish) {
if (boolfield[ijk] == true) {
count = 0;
m = id[ijk];
vfs = VFmask[ijk];
if (m == 0 && vfs < 1 ) {
denA = denB = 0.0;
nn = ijk-1;
if (id[nn] != 0){
denA += DenA[nn];
denB += DenB[nn];
count++;
}


nn = ijk+1;
if (id[nn] != 0){
denA += DenA[nn];
denB += DenB[nn];
count++;
}


nn = ijk-strideY;
if (id[nn] != 0){
denA += DenA[nn];
denB += DenB[nn];
count++;
}


nn = ijk+strideY;
if (id[nn] != 0){
denA += DenA[nn];
denB += DenB[nn];
count++;
}


nn = ijk-strideZ;
if (id[nn] != 0){
denA += DenA[nn];
denB += DenB[nn];
count++;
}


nn = ijk+strideZ;
if (id[nn] != 0){
denA += DenA[nn];
denB += DenB[nn];
count++;
}

if (count == 0) count = 1;
denA /= double(count);
denB /= double(count);

DenA[ijk] = denA;
DenB[ijk] = denB;
if (denA == 0 && denB == 0) denA = 1.0;
Phi[ijk] = (denA-denB)/(denA+denB) * (1-vfs) + affinity * vfs;
}
}
else {

}
}
}
}

__global__  void dvc_ScaLBL_D3Q7_PhaseField_LIBB(int *neighborList,  int *Map, double *Aq, double *Bq, double *DenA, double *DenB,  double *Phi, int start, int finish, int Np, int N, int strideY, int strideZ)
{

int idx,n,nn,nread;
double fAq,fBq,nA,nB;
double fA1,fA2,fA3,fA4,fA5,fA6;
double fB1,fB2,fB3,fB4,fB5,fB6;
double denA,denB;
int nr1,nr2,nr3,nr4,nr5,nr6;

double rho_early;
double rho_diff;
double rho;
double delta;


int S = Np/NBLOCKS/NTHREADS + 1;
for (int s=0; s<S; s++) {
//........Get 1-D index for this thread....................
n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
if (n<finish) {

idx = Map[n];

/* Neighbors */
nr1 = neighborList[n];
nr2 = neighborList[n+Np];
nr3 = neighborList[n+2*Np];
nr4 = neighborList[n+3*Np];
nr5 = neighborList[n+4*Np];
nr6 = neighborList[n+5*Np];

//        fAq = Aq[n];
//        fA1 = Aq[n+Np];
//        fA2 = Aq[n+2*Np];
//        fA3 = Aq[n+3*Np];
//        fA4 = Aq[n+4*Np];
//        fA5 = Aq[n+5*Np];
//        fA6 = Aq[n+6*Np];
//        rho_early = fAq + fA1 + fA2 + fA3 + fA4 + fA5 + fA6;


/* MTE fluid A */
fAq = Aq[n];
fA1 = Aq[nr1];
fA2 = Aq[nr2];
fA3 = Aq[nr3];
fA4 = Aq[nr4];
fA5 = Aq[nr5];
fA6 = Aq[nr6];
// rho_early = fAq + fA1 + fA2 + fA3 + fA4 + fA5 + fA6;

//        fBq = Bq[n];
//        fB1 = Bq[n+Np];
//        fB2 = Bq[n+2*Np];
//        fB3 = Bq[n+3*Np];
//        fB4 = Bq[n+4*Np];
//        fB5 = Bq[n+5*Np];
//        fB6 = Bq[n+6*Np];

/* MTE fluid B */
fBq = Bq[n];
fB1 = Bq[nr1];
fB2 = Bq[nr2];
fB3 = Bq[nr3];
fB4 = Bq[nr4];
fB5 = Bq[nr5];
fB6 = Bq[nr6];
// rho_early = fBq + fB1 + fB2 + fB3 + fB4 + fB5 + fB6;

//        momA = fA1 - fA2 + fA3 - fA4 + fA5 - fA6;
nA = fAq + fA1 + fA2 + fA3 + fA4 + fA5 + fA6;  // rho = fq + f1 + f2 + f3 + f4 + f5 + f6;
//        rho_diff = rho - rho_early;

//        momB = fB1 - fB2 + fB3 - fB4 + fB5 - fB6;
nB = fBq + fB1 + fB2 + fB3 + fB4 + fB5 + fB6;  // rho = fq + f1 + f2 + f3 + f4 + f5 + f6;
//        rho_diff = rho - rho_early;
//        // Bq[n] -= rho_diff;

//        double nAB = nA + nB;
//        double mAB = momA + momB;
//        double err_nAB = 1.0-nAB;
// if (err_nAB > 1E-5) std::cout << "err_nAB=" << err_nAB << "    mAB=" << mAB << std::endl; // Sensible for diffusive problem but inadequate for CG active
// std::cout << "nA=" << nA << ",nB=" << nB << ",momA=" << momA << ",momB=" << momB << "|  " << nAB << "," << mAB << std::endl;

// save the number densities
DenA[idx] = nA;//+err_nAB*0.5;
DenB[idx] = nB;//+err_nAB*0.5;

// save the phase indicator field
Phi[idx] = (nA-nB)/(nA+nB);

}
}
}

__global__ void dvc_ScaLBL_D3Q19_Color_LIBB(int * interpolationList, int *neighborList, int *Map, double *dist, double *dist2, double *savedfq, double *Aq, double *Bq, double *DenA, double *DenB, double *Phi, double *Vel, double *Press, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta, double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np, int N, double * LIBBqA, double * LIBBqBC, double * LIBBqD) {



int n,nn,ijk,nread;
int nr1,nr2,nr3,nr4,nr5,nr6;
int nr7,nr8,nr9,nr10;
int nr11,nr12,nr13,nr14;
int nr15,nr16,nr17,nr18;

double fq;
// conserved momemnts
double rho,jx,jy,jz;
// non-conserved moments
double m0,m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
double m3,m5,m7;
double nA,nB; // number density
double a1,b1,a2,b2,nAB,delta;
double C,nx,ny,nz; //color gradient magnitude and direction
double ux,uy,uz;
double phi,tau,rho0,rlx_setA,rlx_setB;

double f0;
double f1,f2,f3,f4,f5,f6,f7,f8,f9;
double f10,f11,f12,f13,f14,f15,f16,f17,f18;

double g1,g2,g3,g4,g5,g6,g7,g8,g9;
double g10,g11,g12,g13,g14,g15,g16,g17,g18;

double q;
//q = 0.5;

double rho_early;


double A1,A2,A3,A4,A5,A6;
double A7,A8,A9,A10,A11,A12;
double A13,A14,A15,A16,A17,A18;

double BC1,BC2,BC3,BC4,BC5,BC6;
double BC7,BC8,BC9,BC10,BC11,BC12;
double BC13,BC14,BC15,BC16,BC17,BC18;

double D1,D2,D3,D4,D5,D6;
double D7,D8,D9,D10,D11,D12;
double D13,D14,D15,D16,D17,D18;

double r1,r2,r3,r4,r5,r6,r7,r8,r9;

int lr1,lr2,lr3,lr4,lr5,lr6;
int lr7,lr8,lr9,lr10;
int lr11,lr12,lr13,lr14;
int lr15,lr16,lr17,lr18;

double minAx, minBx;
double minAy, minBy;
double minAz, minBz;
double minABx, minABy, minABz;

double spx1, spy1, spz1;
double spx2, spy2, spz2;

double magC1,magC2;

// double nx,ny,nz;

double tmpx1, tmpy1, tmpz1;
double tmpx2, tmpy2, tmpz2;



int S = Np/NBLOCKS/NTHREADS + 1;
for (int s=0; s<S; s++) {
//........Get 1-D index for this thread....................
n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
if (n<finish) {


ijk = Map[n];
/* Accessed in same way as neighborList */
lr1 = interpolationList[n];
lr2 = interpolationList[n+Np];
lr3 = interpolationList[n+2*Np];
lr4 = interpolationList[n+3*Np];
lr5 = interpolationList[n+4*Np];
lr6 = interpolationList[n+5*Np];
lr7 = interpolationList[n+6*Np];
lr8 = interpolationList[n+7*Np];
lr9 = interpolationList[n+8*Np];
lr10 = interpolationList[n+9*Np];
lr11 = interpolationList[n+10*Np];
lr12 = interpolationList[n+11*Np];
lr13 = interpolationList[n+12*Np];
lr14 = interpolationList[n+13*Np];
lr15 = interpolationList[n+14*Np];
lr16 = interpolationList[n+15*Np];
lr17 = interpolationList[n+16*Np];
lr18 = interpolationList[n+17*Np];

nr1 = neighborList[n];
nr2 = neighborList[n+Np];
nr3 = neighborList[n+2*Np];
nr4 = neighborList[n+3*Np];
nr5 = neighborList[n+4*Np];
nr6 = neighborList[n+5*Np];
nr7 = neighborList[n+6*Np];
nr8 = neighborList[n+7*Np];
nr9 = neighborList[n+8*Np];
nr10 = neighborList[n+9*Np];
nr11 = neighborList[n+10*Np];
nr12 = neighborList[n+11*Np];
nr13 = neighborList[n+12*Np];
nr14 = neighborList[n+13*Np];
nr15 = neighborList[n+14*Np];
nr16 = neighborList[n+15*Np];
nr17 = neighborList[n+16*Np];
nr18 = neighborList[n+17*Np];



A2 = LIBBqA[ijk + 0*N];
BC2 = LIBBqBC[ijk + 0*N];
D2 = LIBBqD[ijk + 0*N];

A1 = LIBBqA[ijk + 1*N];
BC1 = LIBBqBC[ijk + 1*N];
D1 = LIBBqD[ijk + 1*N];

A4 = LIBBqA[ijk + 2*N];
BC4 = LIBBqBC[ijk + 2*N];
D4 = LIBBqD[ijk + 2*N];

A3 = LIBBqA[ijk + 3*N];
BC3 = LIBBqBC[ijk + 3*N];
D3 = LIBBqD[ijk + 3*N];

A6 = LIBBqA[ijk + 4*N];
BC6 = LIBBqBC[ijk + 4*N];
D6 = LIBBqD[ijk + 4*N];

A5 = LIBBqA[ijk + 5*N];
BC5 = LIBBqBC[ijk + 5*N];
D5 = LIBBqD[ijk + 5*N];

A8 = LIBBqA[ijk + 6*N];
BC8 = LIBBqBC[ijk + 6*N];
D8 = LIBBqD[ijk + 6*N];

A7 = LIBBqA[ijk + 7*N];
BC7 = LIBBqBC[ijk + 7*N];
D7 = LIBBqD[ijk + 7*N];

A10 = LIBBqA[ijk + 8*N];
BC10 = LIBBqBC[ijk + 8*N];
D10 = LIBBqD[ijk + 8*N];

A9 = LIBBqA[ijk + 9*N];
BC9 = LIBBqBC[ijk + 9*N];
D9= LIBBqD[ijk + 9*N];

A12 = LIBBqA[ijk + 10*N];
BC12 = LIBBqBC[ijk + 10*N];
D12 = LIBBqD[ijk + 10*N];

A11 = LIBBqA[ijk + 11*N];
BC11 = LIBBqBC[ijk + 11*N];
D11 = LIBBqD[ijk + 11*N];

A14 = LIBBqA[ijk + 12*N];
BC14 = LIBBqBC[ijk + 12*N];
D14 = LIBBqD[ijk + 12*N];

A13 = LIBBqA[ijk + 13*N];
BC13 = LIBBqBC[ijk + 13*N];
D13 = LIBBqD[ijk + 13*N];

A16 = LIBBqA[ijk + 14*N];
BC16 = LIBBqBC[ijk + 14*N];
D16 = LIBBqD[ijk + 14*N];

A15 = LIBBqA[ijk + 15*N];
BC15 = LIBBqBC[ijk + 15*N];
D15 = LIBBqD[ijk + 15*N];

A18 = LIBBqA[ijk + 16*N];
BC18 = LIBBqBC[ijk + 16*N];
D18 = LIBBqD[ijk + 16*N];

A17 = LIBBqA[ijk + 17*N];
BC17 = LIBBqBC[ijk + 17*N];
D17 = LIBBqD[ijk + 17*N];

// q = 1
fq = dist[n];
f1 = dist[nr1];
// Autogenerating LBM block...
f2 = dist[nr2];
f3 = dist[nr3];
f4 = dist[nr4];
f5 = dist[nr5];
f6 = dist[nr6];
f7 = dist[nr7];
f8 = dist[nr8];
f9 = dist[nr9];
f10 = dist[nr10];
f11 = dist[nr11];
f12 = dist[nr12];
f13 = dist[nr13];
f14 = dist[nr14];
f15 = dist[nr15];
f16 = dist[nr16];
f17 = dist[nr17];
f18 = dist[nr18];    //done...
rho_early = fq + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9
+ f10 + f11 + f12 + f13 + f14 + f15 + f16 + f17 + f18;


g1 = savedfq[lr1];
// Autogenerating LBM block...
g2 = savedfq[lr2];
g3 = savedfq[lr3];
g4 = savedfq[lr4];
g5 = savedfq[lr5];
g6 = savedfq[lr6];
g7 = savedfq[lr7];
g8 = savedfq[lr8];
g9 = savedfq[lr9];
g10 = savedfq[lr10];
g11 = savedfq[lr11];
g12 = savedfq[lr12];
g13 = savedfq[lr13];
g14 = savedfq[lr14];
g15 = savedfq[lr15];
g16 = savedfq[lr16];
g17 = savedfq[lr17];
g18 = savedfq[lr18];    //done...


nA = DenA[ijk];
nB = DenB[ijk];
phi=(nA-nB)/(nA+nB);

// local density
rho0=1.0;// + 0.5*(1.0-phi)*(rhoB-rhoA);
// local relaxation time
tau=tauA ;//+ 0.5*(1.0-phi)*(tauB-tauA);
rlx_setA = 1.0/tau;
rlx_setB = 8.*(2.-rlx_setA)/(8.-rlx_setA);

m0 = Phi[ijk];

nn = ijk-1;                            // neighbor index (get convention)
m1 = Phi[nn];                        // get neighbor for phi - 1
//........................................................................
nn = ijk+1;                            // neighbor index (get convention)
m2 = Phi[nn];                        // get neighbor for phi - 2
//........................................................................
nn = ijk-strideY;                            // neighbor index (get convention)
m3 = Phi[nn];                    // get neighbor for phi - 3
//........................................................................
nn = ijk+strideY;                            // neighbor index (get convention)
m4 = Phi[nn];                    // get neighbor for phi - 4
//........................................................................
nn = ijk-strideZ;                        // neighbor index (get convention)
m5 = Phi[nn];                    // get neighbor for phi - 5
//........................................................................
nn = ijk+strideZ;                        // neighbor index (get convention)
m6 = Phi[nn];                    // get neighbor for phi - 6
//........................................................................
nn = ijk-strideY-1;                        // neighbor index (get convention)
m7 = Phi[nn];                    // get neighbor for phi - 7
//........................................................................
nn = ijk+strideY+1;                        // neighbor index (get convention)
m8 = Phi[nn];                    // get neighbor for phi - 8
//........................................................................
nn = ijk+strideY-1;                        // neighbor index (get convention)
m9 = Phi[nn];                    // get neighbor for phi - 9
//........................................................................
nn = ijk-strideY+1;                        // neighbor index (get convention)
m10 = Phi[nn];                    // get neighbor for phi - 10
//........................................................................
nn = ijk-strideZ-1;                        // neighbor index (get convention)
m11 = Phi[nn];                    // get neighbor for phi - 11
//........................................................................
nn = ijk+strideZ+1;                        // neighbor index (get convention)
m12 = Phi[nn];                    // get neighbor for phi - 12
//........................................................................
nn = ijk+strideZ-1;                        // neighbor index (get convention)
m13 = Phi[nn];                    // get neighbor for phi - 13
//........................................................................
nn = ijk-strideZ+1;                        // neighbor index (get convention)
m14 = Phi[nn];                    // get neighbor for phi - 14
//........................................................................
nn = ijk-strideZ-strideY;                    // neighbor index (get convention)
m15 = Phi[nn];                    // get neighbor for phi - 15
//........................................................................
nn = ijk+strideZ+strideY;                    // neighbor index (get convention)
m16 = Phi[nn];                    // get neighbor for phi - 16
//........................................................................
nn = ijk+strideZ-strideY;                    // neighbor index (get convention)
m17 = Phi[nn];                    // get neighbor for phi - 17
//........................................................................
nn = ijk-strideZ+strideY;                    // neighbor index (get convention)
m18 = Phi[nn];                    // get neighbor for phi - 18
// ............Compute the Color Gradient...................................
nx = -(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
ny = -(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
nz = -(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));

// ...........Normalize the Color Gradient.................................
C = sqrt(nx*nx+ny*ny+nz*nz);
if (C==0.0) C=1.0;
nx = nx/C;
ny = ny/C;
nz = nz/C;


f1*= BC1; f1 += (D1*dist[nr2] + A1*g1);
// Autogenerating LBM block...
f2*= BC2; f2 += (D2*dist[nr1] + A2*g2);
f3*= BC3; f3 += (D3*dist[nr4] + A3*g3);
f4*= BC4; f4 += (D4*dist[nr3] + A4*g4);
f5*= BC5; f5 += (D5*dist[nr6] + A5*g5);
f6*= BC6; f6 += (D6*dist[nr5] + A6*g6);
f7*= BC7; f7 += (D7*dist[nr8] + A7*g7);
f8*= BC8; f8 += (D8*dist[nr7] + A8*g8);
f9*= BC9; f9 += (D9*dist[nr10] + A9*g9);
f10*= BC10; f10 += (D10*dist[nr9] + A10*g10);
f11*= BC11; f11 += (D11*dist[nr12] + A11*g11);
f12*= BC12; f12 += (D12*dist[nr11] + A12*g12);
f13*= BC13; f13 += (D13*dist[nr14] + A13*g13);
f14*= BC14; f14 += (D14*dist[nr13] + A14*g14);
f15*= BC15; f15 += (D15*dist[nr16] + A15*g15);
f16*= BC16; f16 += (D16*dist[nr15] + A16*g16);
f17*= BC17; f17 += (D17*dist[nr18] + A17*g17);
f18*= BC18; f18 += (D18*dist[nr17] + A18*g18);    //done...


// q=0
rho = fq;
m1  = -30.0*fq;
m2  = 12.0*fq;

// q = 1
rho += f1;
m1 -= 11.0*f1;
m2 -= 4.0*f1;
jx = f1;
m4 = -4.0*f1;
m9 = 2.0*f1;
m10 = -4.0*f1;

// q = 2
rho += f2;
m1 -= 11.0*(f2);
m2 -= 4.0*(f2);
jx -= f2;
m4 += 4.0*(f2);
m9 += 2.0*(f2);
m10 -= 4.0*(f2);

// q = 3
rho += f3;
m1 -= 11.0*f3;
m2 -= 4.0*f3;
jy = f3;
m6 = -4.0*f3;
m9 -= f3;
m10 += 2.0*f3;
m11 = f3;
m12 = -2.0*f3;

// q = 4
rho+= f4;
m1 -= 11.0*f4;
m2 -= 4.0*f4;
jy -= f4;
m6 += 4.0*f4;
m9 -= f4;
m10 += 2.0*f4;
m11 += f4;
m12 -= 2.0*f4;

// q=5
rho += f5;
m1 -= 11.0*f5;
m2 -= 4.0*f5;
jz = f5;
m8 = -4.0*f5;
m9 -= f5;
m10 += 2.0*f5;
m11 -= f5;
m12 += 2.0*f5;

// q = 6
rho+= f6;
m1 -= 11.0*f6;
m2 -= 4.0*f6;
jz -= f6;
m8 += 4.0*f6;
m9 -= f6;
m10 += 2.0*f6;
m11 -= f6;
m12 += 2.0*f6;

// q=7
rho += f7;
m1 += 8.0*f7;
m2 += f7;
jx += f7;
m4 += f7;
jy += f7;
m6 += f7;
m9  += f7;
m10 += f7;
m11 += f7;
m12 += f7;
m13 = f7;
m16 = f7;
m17 = -f7;

// q = 8
rho += f8;
m1 += 8.0*f8;
m2 += f8;
jx -= f8;
m4 -= f8;
jy -= f8;
m6 -= f8;
m9 += f8;
m10 += f8;
m11 += f8;
m12 += f8;
m13 += f8;
m16 -= f8;
m17 += f8;

// q=9
rho += f9;
m1 += 8.0*f9;
m2 += f9;
jx += f9;
m4 += f9;
jy -= f9;
m6 -= f9;
m9 += f9;
m10 += f9;
m11 += f9;
m12 += f9;
m13 -= f9;
m16 += f9;
m17 += f9;

// q = 10
rho += f10;
m1 += 8.0*f10;
m2 += f10;
jx -= f10;
m4 -= f10;
jy += f10;
m6 += f10;
m9 += f10;
m10 += f10;
m11 += f10;
m12 += f10;
m13 -= f10;
m16 -= f10;
m17 -= f10;

// q=11
rho += f11;
m1 += 8.0*f11;
m2 += f11;
jx += f11;
m4 += f11;
jz += f11;
m8 += f11;
m9 += f11;
m10 += f11;
m11 -= f11;
m12 -= f11;
m15 = f11;
m16 -= f11;
m18 = f11;

// q=12
rho += f12;
m1 += 8.0*f12;
m2 += f12;
jx -= f12;
m4 -= f12;
jz -= f12;
m8 -= f12;
m9 += f12;
m10 += f12;
m11 -= f12;
m12 -= f12;
m15 += f12;
m16 += f12;
m18 -= f12;

// q=13
rho += f13;
m1 += 8.0*f13;
m2 += f13;
jx += f13;
m4 += f13;
jz -= f13;
m8 -= f13;
m9 += f13;
m10 += f13;
m11 -= f13;
m12 -= f13;
m15 -= f13;
m16 -= f13;
m18 -= f13;

// q=14
rho += f14;
m1 += 8.0*f14;
m2 += f14;
jx -= f14;
m4 -= f14;
jz += f14;
m8 += f14;
m9 += f14;
m10 += f14;
m11 -= f14;
m12 -= f14;
m15 -= f14;
m16 += f14;
m18 += f14;

// q=15
rho += f15;
m1 += 8.0*f15;
m2 += f15;
jy += f15;
m6 += f15;
jz += f15;
m8 += f15;
m9 -= 2.0*f15;
m10 -= 2.0*f15;
m14 = f15;
m17 += f15;
m18 -= f15;

// q=16
rho += f16;
m1 += 8.0*f16;
m2 += f16;
jy -= f16;
m6 -= f16;
jz -= f16;
m8 -= f16;
m9 -= 2.0*f16;
m10 -= 2.0*f16;
m14 += f16;
m17 -= f16;
m18 += f16;

// q=17
rho += f17;
m1 += 8.0*f17;
m2 += f17;
jy += f17;
m6 += f17;
jz -= f17;
m8 -= f17;
m9 -= 2.0*f17;
m10 -= 2.0*f17;
m14 -= f17;
m17 += f17;
m18 += f17;

// q=18
rho += f18;
m1 += 8.0*f18;
m2 += f18;
jy -= f18;
m6 -= f18;
jz += f18;
m8 += f18;
m9 -= 2.0*f18;
m10 -= 2.0*f18;
m14 -= f18;
m17 -= f18;
m18 -= f18;
double rho_diff = (rho - rho_early);

Vel[n] = jx + 0.5*Fx;
Vel[Np+n] = jy + 0.5*Fy;
Vel[2*Np+n] = jz + 0.5*Fz;
Press[n] = rho/3.0;

// ENTER MOMENTUM INTO THE MOMENTS, NOT VELOCITY
//..........Toelke, Fruediger et. al. 2006................................
if (C == 0.0)    nx = ny = nz = 0.0;

//        C = nx = ny = nz = 0;
m1 = m1 + rlx_setA*((19*(jx*jx+jy*jy+jz*jz)/rho - 11*rho) -alpha*C - m1);
m2 = m2 + rlx_setA*((3*rho - 5.5*(jx*jx+jy*jy+jz*jz)/rho)- m2);
m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
// pxx
m9 = m9 + rlx_setA*(   ((2*jx*jx-jy*jy-jz*jz)/(rho)) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
// PIxx
m10 = m10 + rlx_setA*( -0.5*((2*jx*jx-jy*jy-jz*jz)/(rho))  - m10);
// pww
m11 = m11 + rlx_setA*(   ((jy*jy-jz*jz)/rho) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
// PIww
m12 = m12 + rlx_setA*( -0.5*((jy*jy-jz*jz)/rho)  - m12);
m13 = m13 + rlx_setA*( (jx*jy/rho) + 0.5*alpha*C*nx*ny - m13);
m14 = m14 + rlx_setA*( (jy*jz/rho) + 0.5*alpha*C*ny*nz - m14);
m15 = m15 + rlx_setA*( (jx*jz/rho) + 0.5*alpha*C*nx*nz - m15);
m16 = m16 + rlx_setB*( - m16);
m17 = m17 + rlx_setB*( - m17);
m18 = m18 + rlx_setB*( - m18);

/* Stokes flow - single phase */
//        C = nx = ny = nz = 0;
//        m1 = m1 + rlx_setA*(- 11*rho - m1);
//        m2 = m2 + rlx_setA*(3*rho - m2);
//        m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
//        m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
//        m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
//        // pxx
//        m9 = m9 + rlx_setA*( - m9);
//        // PIxx
//        m10 = m10 + rlx_setA*( - m10);
//        // pww
//        m11 = m11 + rlx_setA*( - m11);
//        // PIww
//        m12 = m12 + rlx_setA*( - m12);
//        m13 = m13 + rlx_setA*( - m13);
//        m14 = m14 + rlx_setA*( - m14);
//        m15 = m15 + rlx_setA*( - m15);
//        m16 = m16 + rlx_setB*( - m16);
//        m17 = m17 + rlx_setB*( - m17);
//        m18 = m18 + rlx_setB*( - m18);



//.................inverse transformation......................................................
// q = 0
fq = 0. - 0.012531328320802004*m1 + 0.047619047619047616*m2 + 0.05263157894736842*rho;
dist2[n] = fq  - rho_diff;

// q = 1
fq = 0. + 0.1*jx - 0.004594820384294068*m1 - 0.05555555555555555*m10 - 0.015873015873015872*m2 -
0.1*m4 + 0.05555555555555555*m9 + 0.05263157894736842*rho + 0.16666666666666667*Fx;
dist2[n+Np] = fq;

// q = 2
fq =  0. - 0.1*jx - 0.004594820384294068*m1 - 0.05555555555555555*m10 - 0.015873015873015872*m2 +
0.1*m4 + 0.05555555555555555*m9 + 0.05263157894736842*rho - 0.16666666666666667*Fx;
dist2[n+2*Np] = fq;

// q= 3
fq = 0. + 0.1*jy - 0.004594820384294068*m1 + 0.027777777777777776*m10 + 0.0833333333333333333*m11 -
0.0833333333333333333*m12 - 0.015873015873015872*m2 - 0.1*m6 - 0.027777777777777776*m9 +
0.05263157894736842*rho + 0.16666666666666667*Fy;
dist2[n+3*Np] = fq;

fq =  0. - 0.1*jy - 0.004594820384294068*m1 + 0.027777777777777776*m10 +
0.0833333333333333333*m11 - 0.0833333333333333333*m12 - 0.015873015873015872*m2 + 0.1*m6 -
0.027777777777777776*m9 + 0.05263157894736842*rho - 0.16666666666666667*Fy;
dist2[n+4*Np] = fq;


fq =  0. + 0.1*jz - 0.004594820384294068*m1 + 0.027777777777777776*m10 - 0.0833333333333333333*m11 +
0.0833333333333333333*m12 - 0.015873015873015872*m2 - 0.1*m8 - 0.027777777777777776*m9 +
0.05263157894736842*rho + 0.16666666666666667*Fz;
dist2[n+5*Np] = fq;


fq = 0. - 0.1*jz - 0.004594820384294068*m1 + 0.027777777777777776*m10 -
0.0833333333333333333*m11 + 0.0833333333333333333*m12 - 0.015873015873015872*m2 + 0.1*m8 -
0.027777777777777776*m9 + 0.05263157894736842*rho - 0.16666666666666667*Fz;
dist2[n+6*Np] = fq;

fq =  0. + 0.1*jx + 0.1*jy + 0.003341687552213868*m1 + 0.013888888888888888*m10 +
0.0833333333333333333*m11 + 0.041666666666666664*m12 + 0.25*m13 + 0.125*m16 - 0.125*m17 +
0.003968253968253968*m2 + 0.025*m4 + 0.025*m6 + 0.027777777777777776*m9 +
0.05263157894736842*rho + 0.08333333333333333*(Fx+Fy);
dist2[n+7*Np] = fq;

fq = 0. - 0.1*jx - 0.1*jy + 0.003341687552213868*m1 +
0.013888888888888888*m10 + 0.0833333333333333333*m11 + 0.041666666666666664*m12 + 0.25*m13 -
0.125*m16 + 0.125*m17 + 0.003968253968253968*m2 - 0.025*m4 - 0.025*m6 +
0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx+Fy);
dist2[n+8*Np] = fq;

fq = 0. + 0.1*jx - 0.1*jy + 0.003341687552213868*m1 + 0.013888888888888888*m10 +
0.0833333333333333333*m11 + 0.041666666666666664*m12 - 0.25*m13 + 0.125*m16 + 0.125*m17 +
0.003968253968253968*m2 + 0.025*m4 - 0.025*m6 + 0.027777777777777776*m9 +
0.05263157894736842*rho + 0.08333333333333333*(Fx-Fy);
dist2[n+9*Np] = fq;

fq = 0. - 0.1*jx + 0.1*jy + 0.003341687552213868*m1 +
0.013888888888888888*m10 + 0.0833333333333333333*m11 + 0.041666666666666664*m12 - 0.25*m13 -
0.125*m16 - 0.125*m17 + 0.003968253968253968*m2 - 0.025*m4 + 0.025*m6 +
0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx-Fy);
dist2[n+10*Np] = fq;

fq =   0. + 0.1*jx + 0.1*jz + 0.003341687552213868*m1 + 0.013888888888888888*m10 -
0.0833333333333333333*m11 - 0.041666666666666664*m12 + 0.25*m15 - 0.125*m16 + 0.125*m18 +
0.003968253968253968*m2 + 0.025*m4 + 0.025*m8 + 0.027777777777777776*m9 +
0.05263157894736842*rho + 0.08333333333333333*(Fx+Fz);
dist2[n+11*Np] = fq;

fq = 0. - 0.1*jx - 0.1*jz + 0.003341687552213868*m1 +
0.013888888888888888*m10 - 0.0833333333333333333*m11 - 0.041666666666666664*m12 + 0.25*m15 +
0.125*m16 - 0.125*m18 + 0.003968253968253968*m2 - 0.025*m4 - 0.025*m8 +
0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx+Fz);
dist2[n+12*Np] = fq;

fq =   0. + 0.1*jx - 0.1*jz + 0.003341687552213868*m1 + 0.013888888888888888*m10 -
0.0833333333333333333*m11 - 0.041666666666666664*m12 - 0.25*m15 - 0.125*m16 - 0.125*m18 +
0.003968253968253968*m2 + 0.025*m4 - 0.025*m8 + 0.027777777777777776*m9 +
0.05263157894736842*rho + 0.08333333333333333*(Fx-Fz);
dist2[n+13*Np] = fq;

fq = 0. - 0.1*jx + 0.1*jz + 0.003341687552213868*m1 +
0.013888888888888888*m10 - 0.0833333333333333333*m11 - 0.041666666666666664*m12 - 0.25*m15 +
0.125*m16 + 0.125*m18 + 0.003968253968253968*m2 - 0.025*m4 + 0.025*m8 +
0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx-Fz);
dist2[n+14*Np] = fq;

fq =  0. + 0.1*jy + 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 + 0.25*m14 +
0.125*m17 - 0.125*m18 + 0.003968253968253968*m2 + 0.025*m6 + 0.025*m8 -
0.05555555555555555*m9 + 0.05263157894736842*rho + 0.08333333333333333*(Fy+Fz);
dist2[n+15*Np] = fq;

fq = 0. - 0.1*jy - 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 + 0.25*m14 -
0.125*m17 + 0.125*m18 + 0.003968253968253968*m2 - 0.025*m6 - 0.025*m8 -
0.05555555555555555*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fy+Fz);
dist2[n+16*Np] = fq;

fq = 0. + 0.1*jy - 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 - 0.25*m14 +
0.125*m17 + 0.125*m18 + 0.003968253968253968*m2 + 0.025*m6 - 0.025*m8 -
0.05555555555555555*m9 + 0.05263157894736842*rho + 0.08333333333333333*(Fy-Fz);
dist2[n+17*Np] = fq;

fq = 0. - 0.1*jy + 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 - 0.25*m14 -
0.125*m17 - 0.125*m18 + 0.003968253968253968*m2 - 0.025*m6 + 0.025*m8 -
0.05555555555555555*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fy-Fz);
dist2[n+18*Np] = fq;

ux = jx / rho0;
uy = jy / rho0;
uz = jz / rho0;

nAB = 1.0/(nA+nB);
Aq[n] = 0.3333333333333333*nA;
Bq[n] = 0.3333333333333333*nB;

delta = beta*nA*nB*nAB*0.1111111111111111*nx;
if (!(nA*nB*nAB>0)) delta=0;
a1 = nA*(0.1111111111111111*(1+4.5*ux)) + delta;
b1 = nB*(0.1111111111111111*(1+4.5*ux)) - delta;
a2 = nA*(0.1111111111111111*(1-4.5*ux)) - delta;
b2 = nB*(0.1111111111111111*(1-4.5*ux)) + delta;
Aq[n+1*Np] = a1;
Bq[n+1*Np] = b1;
Aq[n+2*Np] = a2;
Bq[n+2*Np] = b2;

delta = beta*nA*nB*nAB*0.1111111111111111*ny;
if (!(nA*nB*nAB>0)) delta=0;
a1 = nA*(0.1111111111111111*(1+4.5*uy)) + delta;
b1 = nB*(0.1111111111111111*(1+4.5*uy)) - delta;
a2 = nA*(0.1111111111111111*(1-4.5*uy)) - delta;
b2 = nB*(0.1111111111111111*(1-4.5*uy)) + delta;
Aq[n+3*Np] = a1;
Bq[n+3*Np] = b1;
Aq[n+4*Np] = a2;
Bq[n+4*Np] = b2;

delta = beta*nA*nB*nAB*0.1111111111111111*nz;
if (!(nA*nB*nAB>0)) delta=0;
a1 = nA*(0.1111111111111111*(1+4.5*uz)) + delta;
b1 = nB*(0.1111111111111111*(1+4.5*uz)) - delta;
a2 = nA*(0.1111111111111111*(1-4.5*uz)) - delta;
b2 = nB*(0.1111111111111111*(1-4.5*uz)) + delta;
Aq[n+5*Np] = a1;
Bq[n+5*Np] = b1;
Aq[n+6*Np] = a2;
Bq[n+6*Np] = b2;







}
}
}









__global__ void dvc_ScaLBL_PhaseField_Init_LIBB(int *Map, double *Phi, double *DenA, double * DenB,  double *Aq, double *Bq, int start, int finish, int Np){
int idx,n;
double phi,nA,nB;

int S = Np/NBLOCKS/NTHREADS + 1;
for (int s=0; s<S; s++) {
//........Get 1-D index for this thread....................
idx =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
if (idx<finish) {

n = Map[idx];
phi = Phi[n];
if (phi > 0.0){
nA = 1.0; nB = 0.0;
}
else if (phi < 0.0){
nB = 1.0; nA = 0.0;
}
else{
nA=0.5*(phi+1.0);
nB=0.5*(1.0-phi);
}
DenA[n] = nA;
DenB[n] = nB;

Aq[idx] = 0.3333333333333333*nA;
Aq[1*Np+idx] = (1./9.)*nA;
Aq[2*Np+idx] = (1./9.)*nA;
Aq[3*Np+idx] = (1./9.)*nA;
Aq[4*Np+idx] = (1./9.)*nA;
Aq[5*Np+idx] = (1./9.)*nA;
Aq[6*Np+idx] = (1./9.)*nA;

Bq[idx] = 0.3333333333333333*nB;
Bq[1*Np+idx] = (1./9.)*nB;
Bq[2*Np+idx] = (1./9.)*nB;
Bq[3*Np+idx] = (1./9.)*nB;
Bq[4*Np+idx] = (1./9.)*nB;
Bq[5*Np+idx] = (1./9.)*nB;
Bq[6*Np+idx] = (1./9.)*nB;






}
}
}




__global__ void dvc_LIBB_INTERPOLATION_TEST(int *interpolationList,int *neighborList, int *Map, double *dist, int start, int finish, int Np, int N, double * LIBBqA, double * LIBBqBC, double * LIBBqD, int Nx, int Ny, int Nz,int rank,double * phi) {

int n,ijk,i,j,k;

int nr1,nr2,nr3,nr4,nr5,nr6;
int nr7,nr8,nr9,nr10;
int nr11,nr12,nr13,nr14;
int nr15,nr16,nr17,nr18;

int lr1,lr2,lr3,lr4,lr5,lr6;
int lr7,lr8,lr9,lr10;
int lr11,lr12,lr13,lr14;
int lr15,lr16,lr17,lr18;

double f1,f2,f3,f4,f5,f6,f7,f8,f9;
double f10,f11,f12,f13,f14,f15,f16,f17,f18;


double A = 0.0;
double BC = 0.0;
double D = 0.0;


int count = 0;
//  for (int n=start; n<finish; n++){

int S = Np/NBLOCKS/NTHREADS + 1;
for (int s=0; s<S; s++) {
//........Get 1-D index for this thread....................
n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
if (n<finish) {
ijk = Map[n];
// k = ijk/(Nx*Ny); j = (ijk-Nx*Ny*k)/Nx; i = ijk-Nx*Ny*k-Nx*j;

nr1 = neighborList[n];
nr2 = neighborList[n+Np];
nr3 = neighborList[n+2*Np];
nr4 = neighborList[n+3*Np];
nr5 = neighborList[n+4*Np];
nr6 = neighborList[n+5*Np];
nr7 = neighborList[n+6*Np];
nr8 = neighborList[n+7*Np];
nr9 = neighborList[n+8*Np];
nr10 = neighborList[n+9*Np];
nr11 = neighborList[n+10*Np];
nr12 = neighborList[n+11*Np];
nr13 = neighborList[n+12*Np];
nr14 = neighborList[n+13*Np];
nr15 = neighborList[n+14*Np];
nr16 = neighborList[n+15*Np];
nr17 = neighborList[n+16*Np];
nr18 = neighborList[n+17*Np];

lr1 = interpolationList[n];
lr2 = interpolationList[n+Np];
lr3 = interpolationList[n+2*Np];
lr4 = interpolationList[n+3*Np];
lr5 = interpolationList[n+4*Np];
lr6 = interpolationList[n+5*Np];
lr7 = interpolationList[n+6*Np];
lr8 = interpolationList[n+7*Np];
lr9 = interpolationList[n+8*Np];
lr10 = interpolationList[n+9*Np];
lr11 = interpolationList[n+10*Np];
lr12 = interpolationList[n+11*Np];
lr13 = interpolationList[n+12*Np];
lr14 = interpolationList[n+13*Np];
lr15 = interpolationList[n+14*Np];
lr16 = interpolationList[n+15*Np];
lr17 = interpolationList[n+16*Np];
lr18 = interpolationList[n+17*Np];

/*
// q = 1
A1 = LIBBqA[ijk];
BC1 = LIBBqBC[ijk];
D1 = LIBBqD[ijk];
f1 = BC1*dist[nr1] + D1*dist[nr2] + A1*dist[lr2];
*/

A = LIBBqA[ijk];
BC = LIBBqBC[ijk];
D = LIBBqD[ijk];
dist[nr1] *= BC;
dist[nr1] += (D*dist[lr2] + A*dist[nr2]);  // NOTE THAT LR2 AND NR2 ARE NOW ON D AND A, RESPECTIVELY

A = LIBBqA[ijk+N];
BC = LIBBqBC[ijk+N];
D = LIBBqD[ijk+N];
dist[nr2] *= BC;
dist[nr2] += (D*dist[lr1] + A*dist[nr1]);

A = LIBBqA[ijk+2*N];
BC = LIBBqBC[ijk+2*N];
D = LIBBqD[ijk+2*N];
dist[nr3] *= BC;
dist[nr3] += (D*dist[lr4] + A*dist[nr4]);

A = LIBBqA[ijk+3*N];
BC = LIBBqBC[ijk+3*N];
D = LIBBqD[ijk+3*N];
dist[nr4] *= BC;
dist[nr4] += (D*dist[lr3] + A*dist[nr3]);

A = LIBBqA[ijk+4*N];
BC = LIBBqBC[ijk+4*N];
D = LIBBqD[ijk+4*N];
dist[nr5] *= BC;
dist[nr5] += (D*dist[lr6] + A*dist[nr6]);

A = LIBBqA[ijk+5*N];
BC = LIBBqBC[ijk+5*N];
D = LIBBqD[ijk+5*N];
dist[nr6] *= BC;
dist[nr6] += (D*dist[lr5] + A*dist[nr5]);

A = LIBBqA[ijk+6*N];
BC = LIBBqBC[ijk+6*N];
D = LIBBqD[ijk+6*N];
dist[nr7] *= BC;
dist[nr7] += (D*dist[lr8] + A*dist[nr8]);

A = LIBBqA[ijk+7*N];
BC = LIBBqBC[ijk+7*N];
D = LIBBqD[ijk+7*N];
dist[nr8] *= BC;
dist[nr8] += (D*dist[lr7] + A*dist[nr7]);

A = LIBBqA[ijk+8*N];
BC = LIBBqBC[ijk+8*N];
D = LIBBqD[ijk+8*N];
dist[nr9] *= BC;
dist[nr9] += (D*dist[lr10] + A*dist[nr10]);

A = LIBBqA[ijk+9*N];
BC = LIBBqBC[ijk+9*N];
D = LIBBqD[ijk+9*N];
dist[nr10] *= BC;
dist[nr10] += (D*dist[lr9] + A*dist[nr9]);

A = LIBBqA[ijk+10*N];
BC = LIBBqBC[ijk+10*N];
D = LIBBqD[ijk+10*N];
dist[nr11] *= BC;
dist[nr11] += (D*dist[lr12] + A*dist[nr12]);

A = LIBBqA[ijk+11*N];
BC = LIBBqBC[ijk+11*N];
D = LIBBqD[ijk+11*N];
dist[nr12] *= BC;
dist[nr12] += (D*dist[lr11] + A*dist[nr11]);

A = LIBBqA[ijk+12*N];
BC = LIBBqBC[ijk+12*N];
D = LIBBqD[ijk+12*N];
dist[nr13] *= BC;
dist[nr13] += (D*dist[lr14] + A*dist[nr14]);

A = LIBBqA[ijk+13*N];
BC = LIBBqBC[ijk+13*N];
D = LIBBqD[ijk+13*N];
dist[nr14] *= BC;
dist[nr14] += (D*dist[lr13] + A*dist[nr13]);

A = LIBBqA[ijk+14*N];
BC = LIBBqBC[ijk+14*N];
D = LIBBqD[ijk+14*N];
dist[nr15] *= BC;
dist[nr15] += (D*dist[lr16] + A*dist[nr16]);

A = LIBBqA[ijk+15*N];
BC = LIBBqBC[ijk+15*N];
D = LIBBqD[ijk+15*N];
dist[nr16] *= BC;
dist[nr16] += (D*dist[lr15] + A*dist[nr15]);

A = LIBBqA[ijk+16*N];
BC = LIBBqBC[ijk+16*N];
D = LIBBqD[ijk+16*N];
dist[nr17] *= BC;
dist[nr17] += (D*dist[lr18] + A*dist[nr18]);  //std::cout << " dist[nr17]=" << dist[nr17] << std::endl;

A = LIBBqA[ijk+17*N];
BC = LIBBqBC[ijk+17*N];
D = LIBBqD[ijk+17*N];
dist[nr18] *= BC;
dist[nr18] += (D*dist[lr17] + A*dist[nr17]);

}

}
}




extern "C" void ScaLBL_PhaseField_Init_LIBB(int *Map, double *Phi, double *DenA, double* DenB, double *Aq, double *Bq, int start, int finish, int Np) {

dvc_ScaLBL_PhaseField_Init_LIBB<<<NBLOCKS,NTHREADS >>>(Map,Phi,DenA,DenB,Aq,Bq,start,finish,Np);
}


/*
extern "C" void ScaLBL_D3Q7_PhaseField_LIBB(int *neighborList, int *Map, double *Aq, double *Bq, double *DenA, double * DenB, double *Phi, int start, int finish, int Np, int N, int strideY, int strideZ)
*/

extern "C" void ScaLBL_D3Q7_PhaseField_LIBB(int *neighborList,  int *Map, double *Aq, double *Bq, double *DenA, double *DenB, double *Phi, int start, int finish, int Np, int N, int strideY, int strideZ) {

dvc_ScaLBL_D3Q7_PhaseField_LIBB<<<NBLOCKS,NTHREADS >>>(neighborList, Map, Aq, Bq, DenA, DenB, Phi, start, finish, Np, N, strideY, strideZ);
}



/*
extern "C" void ScaLBL_D3Q19_Color_LIBB(int * interpolationList, int *neighborList, int *Map, double *dist, double *dist2, double *savedfq, double *Aq, double *Bq, double *DenA, double *DenB, double *Phi, double *Vel, double *Press, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta, double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np, int N, double * LIBBqA, double * LIBBqBC, double * LIBBqD)
*/


extern "C" void ScaLBL_D3Q19_Color_LIBB(int * interpolationList, int *neighborList, int *Map, double *dist, double *dist2, double *savedfq, double *Aq, double *Bq, double *DenA, double *DenB, double *Phi, double *Vel, double *Press, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta, double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np, int N, double * LIBBqA, double * LIBBqBC, double * LIBBqD) {

dvc_ScaLBL_D3Q19_Color_LIBB<<<NBLOCKS,NTHREADS >>>(interpolationList,neighborList,
Map, dist, dist2, savedfq, Aq, Bq, DenA,DenB, Phi, Vel, Press, rhoA, rhoB, tauA, tauB, alpha, beta, Fx, Fy, Fz, strideY, strideZ, start,finish, Np, N, LIBBqA, LIBBqBC, LIBBqD);
}


extern "C" void MassAportionAndPhaseUpdate(bool * boolfield, char * id, double affinity, double * Phi, double *DenA, double * DenB, int N, double * VFmask, int Nx, int Ny, int Nz) {

dvc_MassAportionAndPhaseUpdate<<<NBLOCKS,NTHREADS >>>(boolfield, id, affinity, Phi, DenA, DenB, N, VFmask, Nx, Ny, Nz);
}

extern "C" void FiniteVolumeMethod_Init(bool * boolfield, double affinity, double * Phi, double *DenA, double * DenB, int N,  double * VFmask, int Nx, int Ny, int Nz) {

dvc_FiniteVolumeMethod_Init<<<NBLOCKS,NTHREADS >>>(boolfield, affinity,Phi, DenA, DenB, N, VFmask, Nx, Ny, Nz);
}



extern "C" void LIBB_INTERPOLATION_TEST(int *interpolationList,int *neighborList, int *Map, double *dist, int start, int finish, int Np, int N, double * LIBBqA, double * LIBBqBC, double * LIBBqD, int Nx, int Ny, int Nz,int rank,double * phi) {

dvc_LIBB_INTERPOLATION_TEST<<<NBLOCKS,NTHREADS >>>(interpolationList,neighborList, Map, dist, start, finish, Np, N,LIBBqA, LIBBqBC, LIBBqD, Nx, Ny, Nz, rank,phi);

}




__global__ void dvc_ScaLBL_D3Q19_AAodd_Color_LIBB(int * interpolationList, int *neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den,
double *Phi, double *Vel, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np, int N, double * LIBBqA, double * LIBBqBC, double * LIBBqD) {

int n,nn,ijk,nread;
int nr1,nr2,nr3,nr4,nr5,nr6;
int nr7,nr8,nr9,nr10;
int nr11,nr12,nr13,nr14;
int nr15,nr16,nr17,nr18;


double fqNR1,fqNR2,fqNR3,fqNR4,fqNR5,fqNR6,fqNR7,fqNR8,fqNR9;
double fqNR11,fqNR12,fqNR13,fqNR14,fqNR15,fqNR16,fqNR17,fqNR18,fqNR10;

double fq;
// conserved momemnts
double rho,jx,jy,jz;
// non-conserved moments
double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
double m3,m5,m7;
double nA,nB; // number density
double a1,b1,a2,b2,nAB,delta;
double C,nx,ny,nz; //color gradient magnitude and direction
double ux,uy,uz;
double phi,tau,rho0,rlx_setA,rlx_setB;

const double mrt_V1=0.05263157894736842;
const double mrt_V2=0.012531328320802;
const double mrt_V3=0.04761904761904762;
const double mrt_V4=0.004594820384294068;
const double mrt_V5=0.01587301587301587;
const double mrt_V6=0.0555555555555555555555555;
const double mrt_V7=0.02777777777777778;
const double mrt_V8=0.08333333333333333333333;
const double mrt_V9=0.003341687552213868;
const double mrt_V10=0.003968253968253968;
const double mrt_V11=0.01388888888888889;
const double mrt_V12=0.04166666666666666;
double f0;
double f1,f2,f3,f4,f5,f6,f7,f8,f9;
double f10,f11,f12,f13,f14,f15,f16,f17,f18;

double q;
//q = 0.5;

double A1,A2,A3,A4,A5,A6;
double A7,A8,A9,A10,A11,A12;
double A13,A14,A15,A16,A17,A18;

double BC1,BC2,BC3,BC4,BC5,BC6;
double BC7,BC8,BC9,BC10,BC11,BC12;
double BC13,BC14,BC15,BC16,BC17,BC18;

double D1,D2,D3,D4,D5,D6;
double D7,D8,D9,D10,D11,D12;
double D13,D14,D15,D16,D17,D18;


int lr1,lr2,lr3,lr4,lr5,lr6;
int lr7,lr8,lr9,lr10;
int lr11,lr12,lr13,lr14;
int lr15,lr16,lr17,lr18;


double A = 0.0;

int S = Np/NBLOCKS/NTHREADS + 1;
for (int s=0; s<S; s++) {
//........Get 1-D index for this thread....................
n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
if (n<finish) {
ijk = Map[n];

lr1 = interpolationList[n];
lr2 = interpolationList[n+Np];
lr3 = interpolationList[n+2*Np];
lr4 = interpolationList[n+3*Np];
lr5 = interpolationList[n+4*Np];
lr6 = interpolationList[n+5*Np];
lr7 = interpolationList[n+6*Np];
lr8 = interpolationList[n+7*Np];
lr9 = interpolationList[n+8*Np];
lr10 = interpolationList[n+9*Np];
lr11 = interpolationList[n+10*Np];
lr12 = interpolationList[n+11*Np];
lr13 = interpolationList[n+12*Np];
lr14 = interpolationList[n+13*Np];
lr15 = interpolationList[n+14*Np];
lr16 = interpolationList[n+15*Np];
lr17 = interpolationList[n+16*Np];
lr18 = interpolationList[n+17*Np];

nr1 = neighborList[n];
nr2 = neighborList[n+Np];
nr3 = neighborList[n+2*Np];
nr4 = neighborList[n+3*Np];
nr5 = neighborList[n+4*Np];
nr6 = neighborList[n+5*Np];
nr7 = neighborList[n+6*Np];
nr8 = neighborList[n+7*Np];
nr9 = neighborList[n+8*Np];
nr10 = neighborList[n+9*Np];
nr11 = neighborList[n+10*Np];
nr12 = neighborList[n+11*Np];
nr13 = neighborList[n+12*Np];
nr14 = neighborList[n+13*Np];
nr15 = neighborList[n+14*Np];
nr16 = neighborList[n+15*Np];
nr17 = neighborList[n+16*Np];
nr18 = neighborList[n+17*Np];

// read the component number densities
nA = Den[n];
nB = Den[Np + n];

// compute phase indicator field
phi=(nA-nB)/(nA+nB);

// local density
rho0=1 ;//+ 0.5*(1.0-phi)*(rhoB-rhoA);
// local relaxation time
tau=tauA ;//+ 0.5*(1.0-phi)*(tauB-tauA);
rlx_setA = 1.0/tau;
rlx_setB = 8.0*(2.0-rlx_setA)/(8.0-rlx_setA);

// Get the 1D index based on regular data layout
//        ijk = Map[n];
//                    COMPUTE THE COLOR GRADIENT
//........................................................................
//.................Read Phase Indicator Values............................
//........................................................................
nn = ijk-1;                            // neighbor index (get convention)
m1 = Phi[nn];                        // get neighbor for phi - 1
//........................................................................
nn = ijk+1;                            // neighbor index (get convention)
m2 = Phi[nn];                        // get neighbor for phi - 2
//........................................................................
nn = ijk-strideY;                            // neighbor index (get convention)
m3 = Phi[nn];                    // get neighbor for phi - 3
//........................................................................
nn = ijk+strideY;                            // neighbor index (get convention)
m4 = Phi[nn];                    // get neighbor for phi - 4
//........................................................................
nn = ijk-strideZ;                        // neighbor index (get convention)
m5 = Phi[nn];                    // get neighbor for phi - 5
//........................................................................
nn = ijk+strideZ;                        // neighbor index (get convention)
m6 = Phi[nn];                    // get neighbor for phi - 6
//........................................................................
nn = ijk-strideY-1;                        // neighbor index (get convention)
m7 = Phi[nn];                    // get neighbor for phi - 7
//........................................................................
nn = ijk+strideY+1;                        // neighbor index (get convention)
m8 = Phi[nn];                    // get neighbor for phi - 8
//........................................................................
nn = ijk+strideY-1;                        // neighbor index (get convention)
m9 = Phi[nn];                    // get neighbor for phi - 9
//........................................................................
nn = ijk-strideY+1;                        // neighbor index (get convention)
m10 = Phi[nn];                    // get neighbor for phi - 10
//........................................................................
nn = ijk-strideZ-1;                        // neighbor index (get convention)
m11 = Phi[nn];                    // get neighbor for phi - 11
//........................................................................
nn = ijk+strideZ+1;                        // neighbor index (get convention)
m12 = Phi[nn];                    // get neighbor for phi - 12
//........................................................................
nn = ijk+strideZ-1;                        // neighbor index (get convention)
m13 = Phi[nn];                    // get neighbor for phi - 13
//........................................................................
nn = ijk-strideZ+1;                        // neighbor index (get convention)
m14 = Phi[nn];                    // get neighbor for phi - 14
//........................................................................
nn = ijk-strideZ-strideY;                    // neighbor index (get convention)
m15 = Phi[nn];                    // get neighbor for phi - 15
//........................................................................
nn = ijk+strideZ+strideY;                    // neighbor index (get convention)
m16 = Phi[nn];                    // get neighbor for phi - 16
//........................................................................
nn = ijk+strideZ-strideY;                    // neighbor index (get convention)
m17 = Phi[nn];                    // get neighbor for phi - 17
//........................................................................
nn = ijk-strideZ+strideY;                    // neighbor index (get convention)
m18 = Phi[nn];                    // get neighbor for phi - 18
//............Compute the Color Gradient...................................
nx = -(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
ny = -(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
nz = -(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));

// ...........Normalize the Color Gradient.................................
C = sqrt(nx*nx+ny*ny+nz*nz);
if (C==0.0) C=1.0;
nx = nx/C;
ny = ny/C;
nz = nz/C;

// q=0
fq = dist[n];
rho = fq;
m1  = -30.0*fq;
m2  = 12.0*fq;

// q = 1
A1 = LIBBqA[ijk + 0*N];
BC1 = LIBBqBC[ijk + 0*N];
D1 = LIBBqD[ijk + 0*N];
fq = BC1*dist[nr1] + D1*dist[nr2] + A1*dist[lr1];
rho += fq;
m1 -= 11.0*fq;
m2 -= 4.0*fq;
jx = fq;
m4 = -4.0*fq;
m9 = 2.0*fq;
m10 = -4.0*fq;

// q = 2
A2 = LIBBqA[ijk + 1*N];
BC2 = LIBBqBC[ijk + 1*N];
D2 = LIBBqD[ijk + 1*N];
fq = BC2*dist[nr2] + D2*dist[nr1] + A2*dist[lr2];
rho += fq;
m1 -= 11.0*(fq);
m2 -= 4.0*(fq);
jx -= fq;
m4 += 4.0*(fq);
m9 += 2.0*(fq);
m10 -= 4.0*(fq);

// q = 3
A3 = LIBBqA[ijk + 2*N];
BC3 = LIBBqBC[ijk + 2*N];
D3 = LIBBqD[ijk + 2*N];
fq = BC3*dist[nr3] + D3*dist[nr4] + A3*dist[lr3];
rho += fq;
m1 -= 11.0*fq;
m2 -= 4.0*fq;
jy = fq;
m6 = -4.0*fq;
m9 -= fq;
m10 += 2.0*fq;
m11 = fq;
m12 = -2.0*fq;

// q = 4
A4 = LIBBqA[ijk + 3*N];
BC4 = LIBBqBC[ijk + 3*N];
D4 = LIBBqD[ijk + 3*N];
fq = BC4*dist[nr4] + D4*dist[nr3] + A4*dist[lr4];
rho+= fq;
m1 -= 11.0*fq;
m2 -= 4.0*fq;
jy -= fq;
m6 += 4.0*fq;
m9 -= fq;
m10 += 2.0*fq;
m11 += fq;
m12 -= 2.0*fq;

// q=5
A5 = LIBBqA[ijk + 4*N];
BC5 = LIBBqBC[ijk + 4*N];
D5 = LIBBqD[ijk + 4*N];
fq = BC5*dist[nr5] + D5*dist[nr6] + A5*dist[lr5];
rho += fq;
m1 -= 11.0*fq;
m2 -= 4.0*fq;
jz = fq;
m8 = -4.0*fq;
m9 -= fq;
m10 += 2.0*fq;
m11 -= fq;
m12 += 2.0*fq;

// q = 6
A6 = LIBBqA[ijk + 5*N];
BC6 = LIBBqBC[ijk + 5*N];
D6 = LIBBqD[ijk + 5*N];
fq = BC6*dist[nr6] + D6*dist[nr5] + A6*dist[lr6];
rho+= fq;
m1 -= 11.0*fq;
m2 -= 4.0*fq;
jz -= fq;
m8 += 4.0*fq;
m9 -= fq;
m10 += 2.0*fq;
m11 -= fq;
m12 += 2.0*fq;

// q=7
A7 = LIBBqA[ijk + 6*N];
BC7 = LIBBqBC[ijk + 6*N];
D7 = LIBBqD[ijk + 6*N];
fq = BC7*dist[nr7] + D7*dist[nr8] + A7*dist[lr7];
rho += fq;
m1 += 8.0*fq;
m2 += fq;
jx += fq;
m4 += fq;
jy += fq;
m6 += fq;
m9  += fq;
m10 += fq;
m11 += fq;
m12 += fq;
m13 = fq;
m16 = fq;
m17 = -fq;

// q = 8
A8 = LIBBqA[ijk + 7*N];
BC8 = LIBBqBC[ijk + 7*N];
D8 = LIBBqD[ijk + 7*N];
fq = BC8*dist[nr8] + D8*dist[nr7] + A8*dist[lr8];
rho += fq;
m1 += 8.0*fq;
m2 += fq;
jx -= fq;
m4 -= fq;
jy -= fq;
m6 -= fq;
m9 += fq;
m10 += fq;
m11 += fq;
m12 += fq;
m13 += fq;
m16 -= fq;
m17 += fq;

// q=9
A9 = LIBBqA[ijk + 8*N];
BC9 = LIBBqBC[ijk + 8*N];
D9 = LIBBqD[ijk + 8*N];
fq = BC9*dist[nr9] + D9*dist[nr10] + A9*dist[lr9];
rho += fq;
m1 += 8.0*fq;
m2 += fq;
jx += fq;
m4 += fq;
jy -= fq;
m6 -= fq;
m9 += fq;
m10 += fq;
m11 += fq;
m12 += fq;
m13 -= fq;
m16 += fq;
m17 += fq;

// q = 10
A10 = LIBBqA[ijk + 9*N];
BC10 = LIBBqBC[ijk + 9*N];
D10= LIBBqD[ijk + 9*N];
fq = BC10*dist[nr10] + D10*dist[nr9] + A10*dist[lr10];
rho += fq;
m1 += 8.0*fq;
m2 += fq;
jx -= fq;
m4 -= fq;
jy += fq;
m6 += fq;
m9 += fq;
m10 += fq;
m11 += fq;
m12 += fq;
m13 -= fq;
m16 -= fq;
m17 -= fq;

// q=11
A11 = LIBBqA[ijk + 10*N];
BC11 = LIBBqBC[ijk + 10*N];
D11 = LIBBqD[ijk + 10*N];
fq = BC11*dist[nr11] + D11*dist[nr12] + A11*dist[lr11];
rho += fq;
m1 += 8.0*fq;
m2 += fq;
jx += fq;
m4 += fq;
jz += fq;
m8 += fq;
m9 += fq;
m10 += fq;
m11 -= fq;
m12 -= fq;
m15 = fq;
m16 -= fq;
m18 = fq;

// q=12
A12 = LIBBqA[ijk + 11*N];
BC12 = LIBBqBC[ijk + 11*N];
D12 = LIBBqD[ijk + 11*N];
fq = BC12*dist[nr12] + D12*dist[nr11] + A12*dist[lr12];
rho += fq;
m1 += 8.0*fq;
m2 += fq;
jx -= fq;
m4 -= fq;
jz -= fq;
m8 -= fq;
m9 += fq;
m10 += fq;
m11 -= fq;
m12 -= fq;
m15 += fq;
m16 += fq;
m18 -= fq;

// q=13
A13 = LIBBqA[ijk + 12*N];
BC13 = LIBBqBC[ijk + 12*N];
D13 = LIBBqD[ijk + 12*N];
fq = BC13*dist[nr13] + D13*dist[nr14] + A13*dist[lr13];
rho += fq;
m1 += 8.0*fq;
m2 += fq;
jx += fq;
m4 += fq;
jz -= fq;
m8 -= fq;
m9 += fq;
m10 += fq;
m11 -= fq;
m12 -= fq;
m15 -= fq;
m16 -= fq;
m18 -= fq;

// q=14
A14 = LIBBqA[ijk + 13*N];
BC14 = LIBBqBC[ijk + 13*N];
D14 = LIBBqD[ijk + 13*N];
fq = BC14*dist[nr14] + D14*dist[nr13] + A14*dist[lr14];
rho += fq;
m1 += 8.0*fq;
m2 += fq;
jx -= fq;
m4 -= fq;
jz += fq;
m8 += fq;
m9 += fq;
m10 += fq;
m11 -= fq;
m12 -= fq;
m15 -= fq;
m16 += fq;
m18 += fq;

// q=15
A15 = LIBBqA[ijk + 14*N];
BC15 = LIBBqBC[ijk + 14*N];
D15 = LIBBqD[ijk + 14*N];
fq = BC15*dist[nr15] + D15*dist[nr16] + A15*dist[lr15];
rho += fq;
m1 += 8.0*fq;
m2 += fq;
jy += fq;
m6 += fq;
jz += fq;
m8 += fq;
m9 -= 2.0*fq;
m10 -= 2.0*fq;
m14 = fq;
m17 += fq;
m18 -= fq;

// q=16
A16 = LIBBqA[ijk + 15*N];
BC16 = LIBBqBC[ijk + 15*N];
D16 = LIBBqD[ijk + 15*N];
fq = BC16*dist[nr16] + D16*dist[nr15] + A16*dist[lr16];
rho += fq;
m1 += 8.0*fq;
m2 += fq;
jy -= fq;
m6 -= fq;
jz -= fq;
m8 -= fq;
m9 -= 2.0*fq;
m10 -= 2.0*fq;
m14 += fq;
m17 -= fq;
m18 += fq;

// q=17
A17 = LIBBqA[ijk + 16*N];
BC17 = LIBBqBC[ijk + 16*N];
D17 = LIBBqD[ijk + 16*N];
fq = BC17*dist[nr17] + D17*dist[nr18] + A17*dist[lr17];
rho += fq;
m1 += 8.0*fq;
m2 += fq;
jy += fq;
m6 += fq;
jz -= fq;
m8 -= fq;
m9 -= 2.0*fq;
m10 -= 2.0*fq;
m14 -= fq;
m17 += fq;
m18 += fq;

// q=18
A18 = LIBBqA[ijk + 17*N];
BC18 = LIBBqBC[ijk + 17*N];
D18 = LIBBqD[ijk + 17*N];
fq = BC18*dist[nr18] + D18*dist[nr17] + A18*dist[lr18];
rho += fq;
m1 += 8.0*fq;
m2 += fq;
jy -= fq;
m6 -= fq;
jz += fq;
m8 += fq;
m9 -= 2.0*fq;
m10 -= 2.0*fq;
m14 -= fq;
m17 -= fq;
m18 -= fq;




//........................................................................
//..............carry out relaxation process..............................
//..........Toelke, Fruediger et. al. 2006................................
if (C == 0.0)    nx = ny = nz = 0.0;
//        C = nx = ny = nz = 0;
//        m1 = m1 + rlx_setA*((19*(jx*jx+jy*jy+jz*jz)/rho - 11*rho) -alpha*C - m1);
//        m2 = m2 + rlx_setA*((3*rho - 5.5*(jx*jx+jy*jy+jz*jz)/rho)- m2);
//        m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
//        m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
//        m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
//        // pxx
//        m9 = m9 + rlx_setA*(   ((2*jx*jx-jy*jy-jz*jz)/(rho)) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
//        // PIxx
//        m10 = m10 + rlx_setA*( -0.5*((2*jx*jx-jy*jy-jz*jz)/(rho))  - m10);
//        // pww
//        m11 = m11 + rlx_setA*(   ((jy*jy-jz*jz)/rho) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
//        // PIww
//        m12 = m12 + rlx_setA*( -0.5*((jy*jy-jz*jz)/rho)  - m12);
//        m13 = m13 + rlx_setA*( (jx*jy/rho) + 0.5*alpha*C*nx*ny - m13);
//        m14 = m14 + rlx_setA*( (jy*jz/rho) + 0.5*alpha*C*ny*nz - m14);
//        m15 = m15 + rlx_setA*( (jx*jz/rho) + 0.5*alpha*C*nx*nz - m15);
//        m16 = m16 + rlx_setB*( - m16);
//        m17 = m17 + rlx_setB*( - m17);
//        m18 = m18 + rlx_setB*( - m18);
C = nx = ny = nz = 0;
m1 = m1 + rlx_setA*(- 11*rho - m1);
m2 = m2 + rlx_setA*(3*rho - m2);
m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
// pxx
m9 = m9 + rlx_setA*( - m9);
// PIxx
m10 = m10 + rlx_setA*( - m10);
// pww
m11 = m11 + rlx_setA*( - m11);
// PIww
m12 = m12 + rlx_setA*( - m12);
m13 = m13 + rlx_setA*( - m13);
m14 = m14 + rlx_setA*( - m14);
m15 = m15 + rlx_setA*( - m15);
m16 = m16 + rlx_setB*( - m16);
m17 = m17 + rlx_setB*( - m17);
m18 = m18 + rlx_setB*( - m18);
//.................inverse transformation......................................................

fq = 0. - 0.012531328320802004*m1 + 0.047619047619047616*m2 + 0.05263157894736842*rho;
dist[n] = fq;

fq = 0. + 0.1*jx - 0.004594820384294068*m1 - 0.05555555555555555*m10 - 0.015873015873015872*m2 -
0.1*m4 + 0.05555555555555555*m9 + 0.05263157894736842*rho + 0.16666666666666667*Fx;
dist[nr2] = BC2*fq;

fq =  0. - 0.1*jx - 0.004594820384294068*m1 - 0.05555555555555555*m10 - 0.015873015873015872*m2 +
0.1*m4 + 0.05555555555555555*m9 + 0.05263157894736842*rho - 0.16666666666666667*Fx;
dist[nr1] = BC1*fq;

fq = 0. + 0.1*jy - 0.004594820384294068*m1 + 0.027777777777777776*m10 + 0.0833333333333333333*m11 -
0.0833333333333333333*m12 - 0.015873015873015872*m2 - 0.1*m6 - 0.027777777777777776*m9 +
0.05263157894736842*rho + 0.16666666666666667*Fy;
dist[nr4] = BC4*fq;

fq =  0. - 0.1*jy - 0.004594820384294068*m1 + 0.027777777777777776*m10 +
0.0833333333333333333*m11 - 0.0833333333333333333*m12 - 0.015873015873015872*m2 + 0.1*m6 -
0.027777777777777776*m9 + 0.05263157894736842*rho - 0.16666666666666667*Fy;
dist[nr3] = BC3*fq;

fq =  0. + 0.1*jz - 0.004594820384294068*m1 + 0.027777777777777776*m10 - 0.0833333333333333333*m11 +
0.0833333333333333333*m12 - 0.015873015873015872*m2 - 0.1*m8 - 0.027777777777777776*m9 +
0.05263157894736842*rho + 0.16666666666666667*Fz;
dist[nr6] = BC6*fq;

fq = 0. - 0.1*jz - 0.004594820384294068*m1 + 0.027777777777777776*m10 -
0.0833333333333333333*m11 + 0.0833333333333333333*m12 - 0.015873015873015872*m2 + 0.1*m8 -
0.027777777777777776*m9 + 0.05263157894736842*rho - 0.16666666666666667*Fz;
dist[nr5] = BC5*fq;

fq =  0. + 0.1*jx + 0.1*jy + 0.003341687552213868*m1 + 0.013888888888888888*m10 +
0.0833333333333333333*m11 + 0.041666666666666664*m12 + 0.25*m13 + 0.125*m16 - 0.125*m17 +
0.003968253968253968*m2 + 0.025*m4 + 0.025*m6 + 0.027777777777777776*m9 +
0.05263157894736842*rho + 0.08333333333333333*(Fx+Fy);
dist[nr8] = BC8*fq;

fq = 0. - 0.1*jx - 0.1*jy + 0.003341687552213868*m1 +
0.013888888888888888*m10 + 0.0833333333333333333*m11 + 0.041666666666666664*m12 + 0.25*m13 -
0.125*m16 + 0.125*m17 + 0.003968253968253968*m2 - 0.025*m4 - 0.025*m6 +
0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx+Fy);
dist[nr7] = BC7*fq;

fq = 0. + 0.1*jx - 0.1*jy + 0.003341687552213868*m1 + 0.013888888888888888*m10 +
0.0833333333333333333*m11 + 0.041666666666666664*m12 - 0.25*m13 + 0.125*m16 + 0.125*m17 +
0.003968253968253968*m2 + 0.025*m4 - 0.025*m6 + 0.027777777777777776*m9 +
0.05263157894736842*rho + 0.08333333333333333*(Fx-Fy);
dist[nr10] = BC10*fq;

fq = 0. - 0.1*jx + 0.1*jy + 0.003341687552213868*m1 +
0.013888888888888888*m10 + 0.0833333333333333333*m11 + 0.041666666666666664*m12 - 0.25*m13 -
0.125*m16 - 0.125*m17 + 0.003968253968253968*m2 - 0.025*m4 + 0.025*m6 +
0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx-Fy);
dist[nr9] = BC9*fq;

fq =   0. + 0.1*jx + 0.1*jz + 0.003341687552213868*m1 + 0.013888888888888888*m10 -
0.0833333333333333333*m11 - 0.041666666666666664*m12 + 0.25*m15 - 0.125*m16 + 0.125*m18 +
0.003968253968253968*m2 + 0.025*m4 + 0.025*m8 + 0.027777777777777776*m9 +
0.05263157894736842*rho + 0.08333333333333333*(Fx+Fz);
dist[nr12] = BC12*fq;

fq = 0. - 0.1*jx - 0.1*jz + 0.003341687552213868*m1 +
0.013888888888888888*m10 - 0.0833333333333333333*m11 - 0.041666666666666664*m12 + 0.25*m15 +
0.125*m16 - 0.125*m18 + 0.003968253968253968*m2 - 0.025*m4 - 0.025*m8 +
0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx+Fz);
dist[nr11] = BC11*fq;

fq =   0. + 0.1*jx - 0.1*jz + 0.003341687552213868*m1 + 0.013888888888888888*m10 -
0.0833333333333333333*m11 - 0.041666666666666664*m12 - 0.25*m15 - 0.125*m16 - 0.125*m18 +
0.003968253968253968*m2 + 0.025*m4 - 0.025*m8 + 0.027777777777777776*m9 +
0.05263157894736842*rho + 0.08333333333333333*(Fx-Fz);
dist[nr14] = BC14*fq;

fq = 0. - 0.1*jx + 0.1*jz + 0.003341687552213868*m1 +
0.013888888888888888*m10 - 0.0833333333333333333*m11 - 0.041666666666666664*m12 - 0.25*m15 +
0.125*m16 + 0.125*m18 + 0.003968253968253968*m2 - 0.025*m4 + 0.025*m8 +
0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx-Fz);
dist[nr13] = BC13*fq;

fq =  0. + 0.1*jy + 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 + 0.25*m14 +
0.125*m17 - 0.125*m18 + 0.003968253968253968*m2 + 0.025*m6 + 0.025*m8 -
0.05555555555555555*m9 + 0.05263157894736842*rho + 0.08333333333333333*(Fy+Fz);
dist[nr16] = BC16*fq;

fq =  0. - 0.1*jy - 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 + 0.25*m14 -
0.125*m17 + 0.125*m18 + 0.003968253968253968*m2 - 0.025*m6 - 0.025*m8 -
0.05555555555555555*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fy+Fz);
dist[nr15] = BC15*fq;

fq =   0. + 0.1*jy - 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 - 0.25*m14 +
0.125*m17 + 0.125*m18 + 0.003968253968253968*m2 + 0.025*m6 - 0.025*m8 -
0.05555555555555555*m9 + 0.05263157894736842*rho + 0.08333333333333333*(Fy-Fz);
dist[nr18] = BC18*fq;

fq =   0. - 0.1*jy + 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 - 0.25*m14 -
0.125*m17 - 0.125*m18 + 0.003968253968253968*m2 - 0.025*m6 + 0.025*m8 -
0.05555555555555555*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fy-Fz);
dist[nr17] = BC17*fq;


ux = jx/rho0;
ux = jy/rho0;
ux = jz/rho0;

// Instantiate mass transport distributions
// Stationary value - distribution 0
nAB = 1.0/(nA+nB);
Aq[n] = 0.3333333333333333*nA;
Bq[n] = 0.3333333333333333*nB;

//...............................................
// q = 0,2,4
// Cq = {1,0,0}, {0,1,0}, {0,0,1}
delta = beta*nA*nB*nAB*0.1111111111111111*nx;
if (!(nA*nB*nAB>0)) delta=0;
a1 = nA*(0.1111111111111111*(1+4.5*ux))+delta;
b1 = nB*(0.1111111111111111*(1+4.5*ux))-delta;
a2 = nA*(0.1111111111111111*(1-4.5*ux))-delta;
b2 = nB*(0.1111111111111111*(1-4.5*ux))+delta;

// q = 1
Aq[nr2] = a1;
Bq[nr2] = b1;
// q=2
Aq[nr1] = a2;
Bq[nr1] = b2;

//...............................................
// Cq = {0,1,0}
delta = beta*nA*nB*nAB*0.1111111111111111*ny;
if (!(nA*nB*nAB>0)) delta=0;
a1 = nA*(0.1111111111111111*(1+4.5*uy))+delta;
b1 = nB*(0.1111111111111111*(1+4.5*uy))-delta;
a2 = nA*(0.1111111111111111*(1-4.5*uy))-delta;
b2 = nB*(0.1111111111111111*(1-4.5*uy))+delta;

// q = 3
Aq[nr4] = a1;
Bq[nr4] = b1;
// q = 4
Aq[nr3] = a2;
Bq[nr3] = b2;

//...............................................
// q = 4
// Cq = {0,0,1}
delta = beta*nA*nB*nAB*0.1111111111111111*nz;
if (!(nA*nB*nAB>0)) delta=0;
a1 = nA*(0.1111111111111111*(1+4.5*uz))+delta;
b1 = nB*(0.1111111111111111*(1+4.5*uz))-delta;
a2 = nA*(0.1111111111111111*(1-4.5*uz))-delta;
b2 = nB*(0.1111111111111111*(1-4.5*uz))+delta;

// q = 5
Aq[nr6] = a1;
Bq[nr6] = b1;
// q = 6
Aq[nr5] = a2;
Bq[nr5] = b2;
//...............................................

}
}
}


extern "C" void ScaLBL_D3Q19_AAodd_Color_LIBB(int * interpolationList, int *neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den,
double *Phi, double *Vel, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np, int N, double * LIBBqA, double * LIBBqBC, double * LIBBqD)  {


dvc_ScaLBL_D3Q19_AAodd_Color_LIBB<<<NBLOCKS,NTHREADS >>>(interpolationList, neighborList, Map, dist, Aq, Bq, Den,
Phi, Vel,  rhoA,  rhoB,  tauA,  tauB,  alpha,  beta, Fx,  Fy,  Fz, strideY,  strideZ,  start,  finish,  Np,  N,  LIBBqA, LIBBqBC,  LIBBqD);


}

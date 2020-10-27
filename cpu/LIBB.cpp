#include <math.h>
#include <iostream>



extern "C" void FiniteVolumeMethod(double *DenA, double * DenB, double * MassA, double * MassB, int N, double * FVA /* Finite Volume Area */, int strideY, int strideZ, double * VFmask, double * FVDenA, double * FVDenB, int Nx, int Ny, int Nz) {}

extern "C" void ComputeFictitiousPhi(int *Map, char * id, double* Phi_fictitious, double * Phi, int start, int finish, int Np, int N, double * D_coefs, double * DD_coefs, int strideY, int strideZ) {}


extern "C" void FiniteVolumeMethod_Init(bool * boolfield, double affinity, double * Phi, double *DenA, double * DenB, int N,  double * VFmask, int Nx, int Ny, int Nz) {}


extern "C" void ScaLBL_D3Q19_Color_LIBB(int * interpolationList, int *neighborList, int *Map, double *dist, double *dist2, double *savedfq, double *Aq, double *Bq, double *DenA, double *DenB, double *Phi, double *Vel, double *Press, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta, double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np, int N, double * LIBBqA, double * LIBBqBC, double * LIBBqD, double*  GradPhiX, double*GradPhiY, double* GradPhiZ) {
    
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
    
    for (int n=start; n<finish; n++){
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
        
        nn = ijk-1;
        m1 = Phi[nn];                        // get neighbor for phi - 1
        //........................................................................
        nn = ijk+1;
        m2 = Phi[nn];                        // get neighbor for phi - 2
        //........................................................................
        nn = ijk-strideY;
        m3 = Phi[nn];                    // get neighbor for phi - 3
        //........................................................................
        nn = ijk+strideY;
        m4 = Phi[nn];                    // get neighbor for phi - 4
        //........................................................................
        nn = ijk-strideZ;
        m5 = Phi[nn];                    // get neighbor for phi - 5
        //........................................................................
        nn = ijk+strideZ;
        m6 = Phi[nn];                    // get neighbor for phi - 6
        //........................................................................
        nn = ijk-strideY-1;
        m7 = Phi[nn];                    // get neighbor for phi - 7
        //........................................................................
        nn = ijk+strideY+1;
        m8 = Phi[nn];                    // get neighbor for phi - 8
        //........................................................................
        nn = ijk+strideY-1;
        m9 = Phi[nn];                    // get neighbor for phi - 9
        //........................................................................
        nn = ijk-strideY+1;
        m10 = Phi[nn];                    // get neighbor for phi - 10
        //........................................................................
        nn = ijk-strideZ-1;
        m11 = Phi[nn];                    // get neighbor for phi - 11
        //........................................................................
        nn = ijk+strideZ+1;
        m12 = Phi[nn];                    // get neighbor for phi - 12
        //........................................................................
        nn = ijk+strideZ-1;
        m13 = Phi[nn];                    // get neighbor for phi - 13
        //........................................................................
        nn = ijk-strideZ+1;
        m14 = Phi[nn];                    // get neighbor for phi - 14
        //........................................................................
        nn = ijk-strideZ-strideY;
        m15 = Phi[nn];                    // get neighbor for phi - 15
        //........................................................................
        nn = ijk+strideZ+strideY;
        m16 = Phi[nn];                    // get neighbor for phi - 16
        //........................................................................
        nn = ijk+strideZ-strideY;
        m17 = Phi[nn];                    // get neighbor for phi - 17
        //........................................................................
        nn = ijk-strideZ+strideY;
        m18 = Phi[nn];                    // get neighbor for phi - 18
        // ............Compute the Color Gradient...................................
        nx = -(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
        ny = -(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
        nz = -(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
        
        GradPhiX[ijk] = nx;
        GradPhiY[ijk] = ny;
        GradPhiZ[ijk] = nz;
        
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



//        nAB = 1.0/(nA+nB);
//        Aq[n] = 0.3333333333333333*nA;
//        Bq[n] = 0.3333333333333333*nB;
//
//        a1 = nA*(0.1111111111111111*(1+4.5*ux));
//        b1 = nB*(0.1111111111111111*(1+4.5*ux));
//        a2 = nA*(0.1111111111111111*(1-4.5*ux));
//        b2 = nB*(0.1111111111111111*(1-4.5*ux));
//        minAx = std::min(a1,a2);  minBx = std::min(b1,b2);  minABx = std::min(minAx,minBx);
//        tmpx1 = nx*minABx;
//        Aq[n+1*Np] = a1 + tmpx1;
//        Bq[n+1*Np] = b1 - tmpx1;
//        Aq[n+2*Np] = a2 - tmpx1;
//        Bq[n+2*Np] = b2 + tmpx1;
//
//        //...............................................
//        a1 = nA*(0.1111111111111111*(1+4.5*uy));
//        b1 = nB*(0.1111111111111111*(1+4.5*uy));
//        a2 = nA*(0.1111111111111111*(1-4.5*uy));
//        b2 = nB*(0.1111111111111111*(1-4.5*uy));
//        minAy = std::min(a1,a2);  minBy = std::min(b1,b2);  minABy = std::min(minAy,minBy);
//        tmpy1 = ny*minABy;
//        Aq[n+3*Np] = a1 + tmpy1;
//        Bq[n+3*Np] = b1 - tmpy1;
//        Aq[n+4*Np] = a2 - tmpy1;
//        Bq[n+4*Np] = b2 + tmpy1;
//
//        //...............................................
//        a1 = nA*(0.1111111111111111*(1+4.5*uz));
//        b1 = nB*(0.1111111111111111*(1+4.5*uz));
//        a2 = nA*(0.1111111111111111*(1-4.5*uz));
//        b2 = nB*(0.1111111111111111*(1-4.5*uz));
//        minAz = std::min(a1,a2);  minBz = std::min(b1,b2);  minABz = std::min(minAz,minBz);
//        tmpz1 = nz*minABz;
//        Aq[n+5*Np] = a1 + tmpz1;
//        Bq[n+5*Np] = b1 - tmpz1;
//        Aq[n+6*Np] = a2 - tmpz1;
//        Bq[n+6*Np] = b2 + tmpz1;
//        //...............................................



extern "C" void ScaLBL_D3Q7_PhaseField_LIBB(int *neighborList, int *Map, double *Aq, double *Bq, double *DenA, double * DenB, double *Phi, int start, int finish, int Np, int N, int strideY, int strideZ) {
    
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
    
    for (int n=start; n<finish; n++) {
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


extern "C" void ScaLBL_D3Q7_AAodd_PhaseField_LIBB(int *neighborList, int *Map, double *Aq, double *Bq,
                                             double *DenA, double * DenB, double *Phi, int start, int finish, int Np){
    
    int idx,n,nread;
        double fAq,fBq,nA,nB;
    
    double gA1,gA2,gA3,gA4,gA5,gA6;
    double gB1,gB2,gB3,gB4,gB5,gB6;
    double fA1,fA2,fA3,fA4,fA5,fA6;
    double fB1,fB2,fB3,fB4,fB5,fB6;
    
    int nr1,nr2,nr3,nr4,nr5,nr6;
    
    for (int n=start; n<finish; n++){
        idx = Map[n];
        //..........Compute the number density for component A............
        /* Neighbors */
           nr1 = neighborList[n];
           nr2 = neighborList[n+Np];
           nr3 = neighborList[n+2*Np];
           nr4 = neighborList[n+3*Np];
           nr5 = neighborList[n+4*Np];
           nr6 = neighborList[n+5*Np];
           
           /* MTE fluid A */
           fAq = Aq[n];
           fA1 = Aq[nr1];
           fA2 = Aq[nr2];
           fA3 = Aq[nr3];
           fA4 = Aq[nr4];
           fA5 = Aq[nr5];
           fA6 = Aq[nr6];
           // rho_early = fAq + fA1 + fA2 + fA3 + fA4 + fA5 + fA6;
           
           /* MTE fluid B */
           fBq = Bq[n];
           fB1 = Bq[nr1];
           fB2 = Bq[nr2];
           fB3 = Bq[nr3];
           fB4 = Bq[nr4];
           fB5 = Bq[nr5];
           fB6 = Bq[nr6];
        
        nA = fAq + fA1 + fA2 + fA3 + fA4 + fA5 + fA6;  // rho = fq + f1 + f2 + f3 + f4 + f5 + f6;
        //        rho_diff = rho - rho_early;
        
        //        momB = fB1 - fB2 + fB3 - fB4 + fB5 - fB6;
        nB = fBq + fB1 + fB2 + fB3 + fB4 + fB5 + fB6;
        
        // save the number densities
        DenA[idx] = nA;
        DenB[idx] = nB;
        
        // save the phase indicator field
        
        Phi[idx] = (nA-nB)/(nA+nB);
    }
}

extern "C" void ScaLBL_D3Q7_AAeven_PhaseField_LIBB(int *Map, double *Aq, double *Bq, double *DenA, double* DenB, double *Phi,  int start, int finish, int Np) {
    int idx,n,nread;
        double fAq,fBq,nA,nB;
    
    double gA1,gA2,gA3,gA4,gA5,gA6;
    double gB1,gB2,gB3,gB4,gB5,gB6;
    double fA1,fA2,fA3,fA4,fA5,fA6;
    double fB1,fB2,fB3,fB4,fB5,fB6;
    
    
    for (int n=start; n<finish; n++){
        idx = Map[n];
        // compute number density for component A
        // q=0
        fAq = Aq[n];
        fA1 = Aq[2*Np+n];
        fA2 = Aq[1*Np+n];
        fA3 = Aq[4*Np+n];
        fA4 = Aq[3*Np+n];
        fA5 = Aq[6*Np+n];
        fA6 = Aq[5*Np+n];
        
        fBq = Bq[n];
        fB1 = Bq[2*Np+n];
        fB2 = Bq[1*Np+n];
        fB3 = Bq[4*Np+n];
        fB4 = Bq[3*Np+n];
        fB5 = Bq[6*Np+n];
        fB6 = Bq[5*Np+n];
        
        nA = fAq + fA1 + fA2 + fA3 + fA4 + fA5 + fA6;
        
        nB = fBq + fB1 + fB2 + fB3 + fB4 + fB5 + fB6;
        
        // save the number densities
        DenA[idx] = nA;
        DenB[idx] = nB;
        
        // save the phase indicator field
        
        Phi[idx] = (nA-nB)/(nA+nB);
    }
}



// Tolke recoloring method for phase separation
//        fA1 += tmpx1;
//        fA2 += tmpx2;
//        fA3 += tmpy1;
//        fA4 += tmpy2;
//        fA5 += tmpz1;
//        fA6 += tmpz2;
//
//        fB1 -= tmpx1;
//        fB2 -= tmpx2;
//        fB3 -= tmpy1;
//        fB4 -= tmpy2;
//        fB5 -= tmpz1;
//        fB6 -= tmpz2;




// Leclaire method of phase separation
//        fA1 += rc1;
//        fA2 += rc2;
//        fA3 += rc3;
//        fA4 += rc4;
//        fA5 += rc5;
//        fA6 += rc6;
//
//        fB1 -= rc1;
//        fB2 -= rc2;
//        fB3 -= rc3;
//        fB4 -= rc4;
//        fB5 -= rc5;
//        fB6 -= rc6;

extern "C" void ScaLBL_D3Q7_PhaseField_LIBB_VERSION2(int* interpolationList, int *neighborList, int *Map, char * id, double *Aq, double *Bq, double *DenA, double * DenB, double *Phi, int start, int finish, int Np, int N, double * LIBBqA, double * LIBBqBC, double* LIBBqD, double * CG1, double * CG2, double * CG3, double * Velocity, int strideY, int strideZ) {}


extern "C" void LIBB_STREAM_TEST(int * interpolationList, int *neighborList, int *Map, double *dist, double * dist2, int start, int finish, int Np, int N, double * LIBBqA, double * LIBBqBC, double * LIBBqD) {}



extern "C" void LIBB_COLLISION_TEST(int *Map, double *dist, int start, int finish, int Np) {}


extern "C" void ScaLBL_D3Q19_AAodd_Color_LIBB(int * interpolationList, int *neighborList, int *Map, double *dist, double *savedfq, double *Aq, double *Bq, double *DenA, double *DenB,
                                              double *Phi, double *Vel, double* Press, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
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
    
    double m0;
    
    double f0;
    double f1,f2,f3,f4,f5,f6,f7,f8,f9;
    double f10,f11,f12,f13,f14,f15,f16,f17,f18;
    
    double g1,g2,g3,g4,g5,g6,g7,g8,g9;
    double g10,g11,g12,g13,g14,g15,g16,g17,g18;
    
    double rho_early;
    
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
    
    for (int n=start; n<finish; n++){
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
        
        nn = ijk-1;
        m1 = Phi[nn];                        // get neighbor for phi - 1
        //........................................................................
        nn = ijk+1;
        m2 = Phi[nn];                        // get neighbor for phi - 2
        //........................................................................
        nn = ijk-strideY;
        m3 = Phi[nn];                    // get neighbor for phi - 3
        //........................................................................
        nn = ijk+strideY;
        m4 = Phi[nn];                    // get neighbor for phi - 4
        //........................................................................
        nn = ijk-strideZ;
        m5 = Phi[nn];                    // get neighbor for phi - 5
        //........................................................................
        nn = ijk+strideZ;
        m6 = Phi[nn];                    // get neighbor for phi - 6
        //........................................................................
        nn = ijk-strideY-1;
        m7 = Phi[nn];                    // get neighbor for phi - 7
        //........................................................................
        nn = ijk+strideY+1;
        m8 = Phi[nn];                    // get neighbor for phi - 8
        //........................................................................
        nn = ijk+strideY-1;
        m9 = Phi[nn];                    // get neighbor for phi - 9
        //........................................................................
        nn = ijk-strideY+1;
        m10 = Phi[nn];                    // get neighbor for phi - 10
        //........................................................................
        nn = ijk-strideZ-1;
        m11 = Phi[nn];                    // get neighbor for phi - 11
        //........................................................................
        nn = ijk+strideZ+1;
        m12 = Phi[nn];                    // get neighbor for phi - 12
        //........................................................................
        nn = ijk+strideZ-1;
        m13 = Phi[nn];                    // get neighbor for phi - 13
        //........................................................................
        nn = ijk-strideZ+1;
        m14 = Phi[nn];                    // get neighbor for phi - 14
        //........................................................................
        nn = ijk-strideZ-strideY;
        m15 = Phi[nn];                    // get neighbor for phi - 15
        //........................................................................
        nn = ijk+strideZ+strideY;
        m16 = Phi[nn];                    // get neighbor for phi - 16
        //........................................................................
        nn = ijk+strideZ-strideY;
        m17 = Phi[nn];                    // get neighbor for phi - 17
        //........................................................................
        nn = ijk-strideZ+strideY;
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
//                C = nx = ny = nz = 0;
//                m1 = m1 + rlx_setA*(- 11*rho - m1);
//                m2 = m2 + rlx_setA*(3*rho - m2);
//                m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
//                m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
//                m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
//                // pxx
//                m9 = m9 + rlx_setA*( - m9);
//                // PIxx
//                m10 = m10 + rlx_setA*( - m10);
//                // pww
//                m11 = m11 + rlx_setA*( - m11);
//                // PIww
//                m12 = m12 + rlx_setA*( - m12);
//                m13 = m13 + rlx_setA*( - m13);
//                m14 = m14 + rlx_setA*( - m14);
//                m15 = m15 + rlx_setA*( - m15);
//                m16 = m16 + rlx_setB*( - m16);
//                m17 = m17 + rlx_setB*( - m17);
//                m18 = m18 + rlx_setB*( - m18);
        //.................inverse transformation......................................................
        
        
        // Compute rho_early prior to stream and send it to AAeven
        
        fq = 0. - 0.012531328320802004*m1 + 0.047619047619047616*m2 + 0.05263157894736842*rho;
        dist[n] = fq;
        
        fq = 0. + 0.1*jx - 0.004594820384294068*m1 - 0.05555555555555555*m10 - 0.015873015873015872*m2 -
        0.1*m4 + 0.05555555555555555*m9 + 0.05263157894736842*rho + 0.16666666666666667*Fx;
        dist[nr2] = fq;
        
        fq =  0. - 0.1*jx - 0.004594820384294068*m1 - 0.05555555555555555*m10 - 0.015873015873015872*m2 +
        0.1*m4 + 0.05555555555555555*m9 + 0.05263157894736842*rho - 0.16666666666666667*Fx;
        dist[nr1] = fq;
        
        fq = 0. + 0.1*jy - 0.004594820384294068*m1 + 0.027777777777777776*m10 + 0.0833333333333333333*m11 -
        0.0833333333333333333*m12 - 0.015873015873015872*m2 - 0.1*m6 - 0.027777777777777776*m9 +
        0.05263157894736842*rho + 0.16666666666666667*Fy;
        dist[nr4] = fq;
        
        fq =  0. - 0.1*jy - 0.004594820384294068*m1 + 0.027777777777777776*m10 +
        0.0833333333333333333*m11 - 0.0833333333333333333*m12 - 0.015873015873015872*m2 + 0.1*m6 -
        0.027777777777777776*m9 + 0.05263157894736842*rho - 0.16666666666666667*Fy;
        dist[nr3] = fq;
        
        fq =  0. + 0.1*jz - 0.004594820384294068*m1 + 0.027777777777777776*m10 - 0.0833333333333333333*m11 +
        0.0833333333333333333*m12 - 0.015873015873015872*m2 - 0.1*m8 - 0.027777777777777776*m9 +
        0.05263157894736842*rho + 0.16666666666666667*Fz;
        dist[nr6] = fq;
        
        fq = 0. - 0.1*jz - 0.004594820384294068*m1 + 0.027777777777777776*m10 -
        0.0833333333333333333*m11 + 0.0833333333333333333*m12 - 0.015873015873015872*m2 + 0.1*m8 -
        0.027777777777777776*m9 + 0.05263157894736842*rho - 0.16666666666666667*Fz;
        dist[nr5] = fq;
        
        fq =  0. + 0.1*jx + 0.1*jy + 0.003341687552213868*m1 + 0.013888888888888888*m10 +
        0.0833333333333333333*m11 + 0.041666666666666664*m12 + 0.25*m13 + 0.125*m16 - 0.125*m17 +
        0.003968253968253968*m2 + 0.025*m4 + 0.025*m6 + 0.027777777777777776*m9 +
        0.05263157894736842*rho + 0.08333333333333333*(Fx+Fy);
        dist[nr8] = fq;
        
        fq = 0. - 0.1*jx - 0.1*jy + 0.003341687552213868*m1 +
        0.013888888888888888*m10 + 0.0833333333333333333*m11 + 0.041666666666666664*m12 + 0.25*m13 -
        0.125*m16 + 0.125*m17 + 0.003968253968253968*m2 - 0.025*m4 - 0.025*m6 +
        0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx+Fy);
        dist[nr7] = fq;
        
        fq = 0. + 0.1*jx - 0.1*jy + 0.003341687552213868*m1 + 0.013888888888888888*m10 +
        0.0833333333333333333*m11 + 0.041666666666666664*m12 - 0.25*m13 + 0.125*m16 + 0.125*m17 +
        0.003968253968253968*m2 + 0.025*m4 - 0.025*m6 + 0.027777777777777776*m9 +
        0.05263157894736842*rho + 0.08333333333333333*(Fx-Fy);
        dist[nr10] = fq;
        
        fq = 0. - 0.1*jx + 0.1*jy + 0.003341687552213868*m1 +
        0.013888888888888888*m10 + 0.0833333333333333333*m11 + 0.041666666666666664*m12 - 0.25*m13 -
        0.125*m16 - 0.125*m17 + 0.003968253968253968*m2 - 0.025*m4 + 0.025*m6 +
        0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx-Fy);
        dist[nr9] = fq;
        
        fq =   0. + 0.1*jx + 0.1*jz + 0.003341687552213868*m1 + 0.013888888888888888*m10 -
        0.0833333333333333333*m11 - 0.041666666666666664*m12 + 0.25*m15 - 0.125*m16 + 0.125*m18 +
        0.003968253968253968*m2 + 0.025*m4 + 0.025*m8 + 0.027777777777777776*m9 +
        0.05263157894736842*rho + 0.08333333333333333*(Fx+Fz);
        dist[nr12] = fq;
        
        fq = 0. - 0.1*jx - 0.1*jz + 0.003341687552213868*m1 +
        0.013888888888888888*m10 - 0.0833333333333333333*m11 - 0.041666666666666664*m12 + 0.25*m15 +
        0.125*m16 - 0.125*m18 + 0.003968253968253968*m2 - 0.025*m4 - 0.025*m8 +
        0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx+Fz);
        dist[nr11] = fq;
        
        fq =   0. + 0.1*jx - 0.1*jz + 0.003341687552213868*m1 + 0.013888888888888888*m10 -
        0.0833333333333333333*m11 - 0.041666666666666664*m12 - 0.25*m15 - 0.125*m16 - 0.125*m18 +
        0.003968253968253968*m2 + 0.025*m4 - 0.025*m8 + 0.027777777777777776*m9 +
        0.05263157894736842*rho + 0.08333333333333333*(Fx-Fz);
        dist[nr14] = fq;
        
        fq = 0. - 0.1*jx + 0.1*jz + 0.003341687552213868*m1 +
        0.013888888888888888*m10 - 0.0833333333333333333*m11 - 0.041666666666666664*m12 - 0.25*m15 +
        0.125*m16 + 0.125*m18 + 0.003968253968253968*m2 - 0.025*m4 + 0.025*m8 +
        0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx-Fz);
        dist[nr13] = fq;
        
        fq =  0. + 0.1*jy + 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 + 0.25*m14 +
        0.125*m17 - 0.125*m18 + 0.003968253968253968*m2 + 0.025*m6 + 0.025*m8 -
        0.05555555555555555*m9 + 0.05263157894736842*rho + 0.08333333333333333*(Fy+Fz);
        dist[nr16] = fq;
        
        fq =  0. - 0.1*jy - 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 + 0.25*m14 -
        0.125*m17 + 0.125*m18 + 0.003968253968253968*m2 - 0.025*m6 - 0.025*m8 -
        0.05555555555555555*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fy+Fz);
        dist[nr15] = fq;
        
        fq =   0. + 0.1*jy - 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 - 0.25*m14 +
        0.125*m17 + 0.125*m18 + 0.003968253968253968*m2 + 0.025*m6 - 0.025*m8 -
        0.05555555555555555*m9 + 0.05263157894736842*rho + 0.08333333333333333*(Fy-Fz);
        dist[nr18] = fq;
        
        fq =   0. - 0.1*jy + 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 - 0.25*m14 -
        0.125*m17 - 0.125*m18 + 0.003968253968253968*m2 - 0.025*m6 + 0.025*m8 -
        0.05555555555555555*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fy-Fz);
        dist[nr17] = fq;
        
        
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



extern "C" void ScaLBL_D3Q19_AAeven_Color_LIBB(int * interpolationList, int *neighborList, int *Map, double *dist, double *savedfq, double *Aq, double *Bq, double *DenA, double* DenB, double *Phi,
                                               double *Vel, double* Press, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
                                               double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np, int N, double * LIBBqA, double * LIBBqBC, double * LIBBqD) {
    
    int ijk,nn,n;
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
    
    double f1,f2,f3,f4,f5,f6,f7,f8,f9;
    double f10,f11,f12,f13,f14,f15,f16,f17,f18;
       
    double g1,g2,g3,g4,g5,g6,g7,g8,g9;
    double g10,g11,g12,g13,g14,g15,g16,g17,g18;
    
    double A1,A2,A3,A4,A5,A6;
    double A7,A8,A9,A10,A11,A12;
    double A13,A14,A15,A16,A17,A18;
    
    double BC1,BC2,BC3,BC4,BC5,BC6;
    double BC7,BC8,BC9,BC10,BC11,BC12;
    double BC13,BC14,BC15,BC16,BC17,BC18;
    
    double D1,D2,D3,D4,D5,D6;
    double D7,D8,D9,D10,D11,D12;
    double D13,D14,D15,D16,D17,D18;
    
    int nr1,nr2,nr3,nr4,nr5,nr6;
    int nr7,nr8,nr9,nr10;
    int nr11,nr12,nr13,nr14;
    int nr15,nr16,nr17,nr18;
    
    int lr1,lr2,lr3,lr4,lr5,lr6;
    int lr7,lr8,lr9,lr10;
    int lr11,lr12,lr13,lr14;
    int lr15,lr16,lr17,lr18;
    
    
    for (int n=start; n<finish; n++){
        ijk = Map[n];
        
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
        
        f1 = dist[n+1*Np];
            // Autogenerating LBM block...
        f2 = dist[n+2*Np];
        f3 = dist[n+3*Np];
        f4 = dist[n+4*Np];
        f5 = dist[n+5*Np];
        f6 = dist[n+6*Np];
        f7 = dist[n+7*Np];
        f8 = dist[n+8*Np];
        f9 = dist[n+9*Np];
        f10 = dist[n+10*Np];
        f11 = dist[n+11*Np];
        f12 = dist[n+12*Np];
        f13 = dist[n+13*Np];
        f14 = dist[n+14*Np];
        f15 = dist[n+15*Np];
        f16 = dist[n+16*Np];
        f17 = dist[n+17*Np];
        f18 = dist[n+18*Np];    //done...
        
        g1 = savedfq[nr1];
        // Autogenerating LBM block...
        g2 = savedfq[nr2];
        g3 = savedfq[nr3];
        g4 = savedfq[nr4];
        g5 = savedfq[nr5];
        g6 = savedfq[nr6];
        g7 = savedfq[nr7];
        g8 = savedfq[nr8];
        g9 = savedfq[nr9];
        g10 = savedfq[nr10];
        g11 = savedfq[nr11];
        g12 = savedfq[nr12];
        g13 = savedfq[nr13];
        g14 = savedfq[nr14];
        g15 = savedfq[nr15];
        g16 = savedfq[nr16];
        g17 = savedfq[nr17];
        g18 = savedfq[nr18];    //done...
        
        nA = DenA[ijk];
        nB = DenB[ijk];
        
        // compute phase indicator field
        phi=(nA-nB)/(nA+nB);
        
        // local density
        rho0=1.0;// + 0.5*(1.0-phi)*(rhoB-rhoA);  // 1.0 for single phase
        // local relaxation time
        tau=tauA;// + 0.5*(1.0-phi)*(tauB-tauA); // tau = tauA for single phase
        rlx_setA = 1.0/tau;
        rlx_setB = 8.*(2.-rlx_setA)/(8.-rlx_setA); // 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
        

        //                    COMPUTE THE COLOR GRADIENT
        //........................................................................
        //.................Read Phase Indicator Values............................
        //........................................................................
        nn = ijk-1;
        m1 = Phi[nn];                        // get neighbor for phi - 1
        //........................................................................
        nn = ijk+1;
        m2 = Phi[nn];                        // get neighbor for phi - 2
        //........................................................................
        nn = ijk-strideY;
        m3 = Phi[nn];                    // get neighbor for phi - 3
        //........................................................................
        nn = ijk+strideY;
        m4 = Phi[nn];                    // get neighbor for phi - 4
        //........................................................................
        nn = ijk-strideZ;
        m5 = Phi[nn];                    // get neighbor for phi - 5
        //........................................................................
        nn = ijk+strideZ;
        m6 = Phi[nn];                    // get neighbor for phi - 6
        //........................................................................
        nn = ijk-strideY-1;
        m7 = Phi[nn];                    // get neighbor for phi - 7
        //........................................................................
        nn = ijk+strideY+1;
        m8 = Phi[nn];                    // get neighbor for phi - 8
        //........................................................................
        nn = ijk+strideY-1;
        m9 = Phi[nn];                    // get neighbor for phi - 9
        //........................................................................
        nn = ijk-strideY+1;
        m10 = Phi[nn];                    // get neighbor for phi - 10
        //........................................................................
        nn = ijk-strideZ-1;
        m11 = Phi[nn];                    // get neighbor for phi - 11
        //........................................................................
        nn = ijk+strideZ+1;
        m12 = Phi[nn];                    // get neighbor for phi - 12
        //........................................................................
        nn = ijk+strideZ-1;
        m13 = Phi[nn];                    // get neighbor for phi - 13
        //........................................................................
        nn = ijk-strideZ+1;
        m14 = Phi[nn];                    // get neighbor for phi - 14
        //........................................................................
        nn = ijk-strideZ-strideY;
        m15 = Phi[nn];                    // get neighbor for phi - 15
        //........................................................................
        nn = ijk+strideZ+strideY;
        m16 = Phi[nn];                    // get neighbor for phi - 16
        //........................................................................
        nn = ijk+strideZ-strideY;
        m17 = Phi[nn];                    // get neighbor for phi - 17
        //........................................................................
        nn = ijk-strideZ+strideY;
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
        
        f2*= BC1; f2 += (D1*dist[n+1*Np] + A1*g2);
        // Autogenerating LBM block...
        f1*= BC2; f1 += (D2*dist[n+2*Np] + A2*g1);
        f4*= BC3; f4 += (D3*dist[n+3*Np] + A3*g4);
        f3*= BC4; f3 += (D4*dist[n+4*Np] + A4*g3);
        f6*= BC5; f6 += (D5*dist[n+5*Np] + A5*g6);
        f5*= BC6; f5 += (D6*dist[n+6*Np] + A6*g5);
        f8*= BC7; f8 += (D7*dist[n+7*Np] + A7*g8);
        f7*= BC8; f7 += (D8*dist[n+8*Np] + A8*g7);
        f10*= BC9; f10 += (D9*dist[n+9*Np] + A9*g10);
        f9*= BC10; f9 += (D10*dist[n+10*Np] + A10*g9);
        f12*= BC11; f12 += (D11*dist[n+11*Np] + A11*g12);
        f11*= BC12; f11 += (D12*dist[n+12*Np] + A12*g11);
        f14*= BC13; f14 += (D13*dist[n+13*Np] + A13*g14);
        f13*= BC14; f13 += (D14*dist[n+14*Np] + A14*g13);
        f16*= BC15; f16 += (D15*dist[n+15*Np] + A15*g16);
        f15*= BC16; f15 += (D16*dist[n+16*Np] + A16*g15);
        f18*= BC17; f18 += (D17*dist[n+17*Np] + A17*g18);
        f17*= BC18; f17 += (D18*dist[n+18*Np] + A18*g17);    //done...
        
        
        // q = 0
        fq = dist[n];
        rho = fq;
        m1  = -30.0*fq;
        m2  = 12.0*fq;
        
        // q = 1
        rho += f2;
        m1 -= 11.0*f2;
        m2 -= 4.0*f2;
        jx = f2;
        m4 = -4.0*f2;
        m9 = 2.0*f2;
        m10 = -4.0*f2;
        
        // q = 2
        rho += f1;
        m1 -= 11.0*(f1);
        m2 -= 4.0*(f1);
        jx -= f1;
        m4 += 4.0*(f1);
        m9 += 2.0*(f1);
        m10 -= 4.0*(f1);
        
        // q = 3
        rho += f4;
        m1 -= 11.0*f4;
        m2 -= 4.0*f4;
        jy = f4;
        m6 = -4.0*f4;
        m9 -= f4;
        m10 += 2.0*f4;
        m11 = f4;
        m12 = -2.0*f4;
        
        // q = 4
        rho+= f3;
        m1 -= 11.0*f3;
        m2 -= 4.0*f3;
        jy -= f3;
        m6 += 4.0*f3;
        m9 -= f3;
        m10 += 2.0*f3;
        m11 += f3;
        m12 -= 2.0*f3;
        
        // q=5
        rho += f6;
        m1 -= 11.0*f6;
        m2 -= 4.0*f6;
        jz = f6;
        m8 = -4.0*f6;
        m9 -= f6;
        m10 += 2.0*f6;
        m11 -= f6;
        m12 += 2.0*f6;
        
        // q = 6
        rho+= f5;
        m1 -= 11.0*f5;
        m2 -= 4.0*f5;
        jz -= f5;
        m8 += 4.0*f5;
        m9 -= f5;
        m10 += 2.0*f5;
        m11 -= f5;
        m12 += 2.0*f5;
        
        // q=7
        rho += f8;
        m1 += 8.0*f8;
        m2 += f8;
        jx += f8;
        m4 += f8;
        jy += f8;
        m6 += f8;
        m9  += f8;
        m10 += f8;
        m11 += f8;
        m12 += f8;
        m13 = f8;
        m16 = f8;
        m17 = -f8;
        
        // q = 8
        rho += f7;
        m1 += 8.0*f7;
        m2 += f7;
        jx -= f7;
        m4 -= f7;
        jy -= f7;
        m6 -= f7;
        m9 += f7;
        m10 += f7;
        m11 += f7;
        m12 += f7;
        m13 += f7;
        m16 -= f7;
        m17 += f7;
        
        // q=9
        rho += f10;
        m1 += 8.0*f10;
        m2 += f10;
        jx += f10;
        m4 += f10;
        jy -= f10;
        m6 -= f10;
        m9 += f10;
        m10 += f10;
        m11 += f10;
        m12 += f10;
        m13 -= f10;
        m16 += f10;
        m17 += f10;
        
        // q = 10
        rho += f9;
        m1 += 8.0*f9;
        m2 += f9;
        jx -= f9;
        m4 -= f9;
        jy += f9;
        m6 += f9;
        m9 += f9;
        m10 += f9;
        m11 += f9;
        m12 += f9;
        m13 -= f9;
        m16 -= f9;
        m17 -= f9;
        
        // q=11
        rho += f12;
        m1 += 8.0*f12;
        m2 += f12;
        jx += f12;
        m4 += f12;
        jz += f12;
        m8 += f12;
        m9 += f12;
        m10 += f12;
        m11 -= f12;
        m12 -= f12;
        m15 = f12;
        m16 -= f12;
        m18 = f12;
        
        // q=12
        rho += f11;
        m1 += 8.0*f11;
        m2 += f11;
        jx -= f11;
        m4 -= f11;
        jz -= f11;
        m8 -= f11;
        m9 += f11;
        m10 += f11;
        m11 -= f11;
        m12 -= f11;
        m15 += f11;
        m16 += f11;
        m18 -= f11;
        
        // q=13
        rho += f14;
        m1 += 8.0*f14;
        m2 += f14;
        jx += f14;
        m4 += f14;
        jz -= f14;
        m8 -= f14;
        m9 += f14;
        m10 += f14;
        m11 -= f14;
        m12 -= f14;
        m15 -= f14;
        m16 -= f14;
        m18 -= f14;
        
        // q=14
        rho += f13;
        m1 += 8.0*f13;
        m2 += f13;
        jx -= f13;
        m4 -= f13;
        jz += f13;
        m8 += f13;
        m9 += f13;
        m10 += f13;
        m11 -= f13;
        m12 -= f13;
        m15 -= f13;
        m16 += f13;
        m18 += f13;
        
        // q=15
        rho += f16;
        m1 += 8.0*f16;
        m2 += f16;
        jy += f16;
        m6 += f16;
        jz += f16;
        m8 += f16;
        m9 -= 2.0*f16;
        m10 -= 2.0*f16;
        m14 = f16;
        m17 += f16;
        m18 -= f16;
        
        // q=16
        rho += f15;
        m1 += 8.0*f15;
        m2 += f15;
        jy -= f15;
        m6 -= f15;
        jz -= f15;
        m8 -= f15;
        m9 -= 2.0*f15;
        m10 -= 2.0*f15;
        m14 += f15;
        m17 -= f15;
        m18 += f15;
        
        // q=17
        rho += f18;
        m1 += 8.0*f18;
        m2 += f18;
        jy += f18;
        m6 += f18;
        jz -= f18;
        m8 -= f18;
        m9 -= 2.0*f18;
        m10 -= 2.0*f18;
        m14 -= f18;
        m17 += f18;
        m18 += f18;
        
        // q=18
        rho += f17;
        m1 += 8.0*f17;
        m2 += f17;
        jy -= f17;
        m6 -= f17;
        jz += f17;
        m8 += f17;
        m9 -= 2.0*f17;
        m10 -= 2.0*f17;
        m14 -= f17;
        m17 -= f17;
        m18 -= f17;

        
        
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
        
        
        fq = 0. - 0.012531328320802004*m1 + 0.047619047619047616*m2 + 0.05263157894736842*rho;
        dist[n] = fq;//  - rho_diff;
        
        // q = 1
        fq = 0. + 0.1*jx - 0.004594820384294068*m1 - 0.05555555555555555*m10 - 0.015873015873015872*m2 -
        0.1*m4 + 0.05555555555555555*m9 + 0.05263157894736842*rho + 0.16666666666666667*Fx;
        dist[n+Np] = fq;
        
        // q = 2
        fq =  0. - 0.1*jx - 0.004594820384294068*m1 - 0.05555555555555555*m10 - 0.015873015873015872*m2 +
        0.1*m4 + 0.05555555555555555*m9 + 0.05263157894736842*rho - 0.16666666666666667*Fx;
        dist[n+2*Np] = fq;
        
        // q= 3
        fq = 0. + 0.1*jy - 0.004594820384294068*m1 + 0.027777777777777776*m10 + 0.0833333333333333333*m11 -
        0.0833333333333333333*m12 - 0.015873015873015872*m2 - 0.1*m6 - 0.027777777777777776*m9 +
        0.05263157894736842*rho + 0.16666666666666667*Fy;
        dist[n+3*Np] = fq;
        
        fq =  0. - 0.1*jy - 0.004594820384294068*m1 + 0.027777777777777776*m10 +
        0.0833333333333333333*m11 - 0.0833333333333333333*m12 - 0.015873015873015872*m2 + 0.1*m6 -
        0.027777777777777776*m9 + 0.05263157894736842*rho - 0.16666666666666667*Fy;
        dist[n+4*Np] = fq;
        
        
        fq =  0. + 0.1*jz - 0.004594820384294068*m1 + 0.027777777777777776*m10 - 0.0833333333333333333*m11 +
        0.0833333333333333333*m12 - 0.015873015873015872*m2 - 0.1*m8 - 0.027777777777777776*m9 +
        0.05263157894736842*rho + 0.16666666666666667*Fz;
        dist[n+5*Np] = fq;
        
        
        fq = 0. - 0.1*jz - 0.004594820384294068*m1 + 0.027777777777777776*m10 -
        0.0833333333333333333*m11 + 0.0833333333333333333*m12 - 0.015873015873015872*m2 + 0.1*m8 -
        0.027777777777777776*m9 + 0.05263157894736842*rho - 0.16666666666666667*Fz;
        dist[n+6*Np] = fq;
        
        fq =  0. + 0.1*jx + 0.1*jy + 0.003341687552213868*m1 + 0.013888888888888888*m10 +
        0.0833333333333333333*m11 + 0.041666666666666664*m12 + 0.25*m13 + 0.125*m16 - 0.125*m17 +
        0.003968253968253968*m2 + 0.025*m4 + 0.025*m6 + 0.027777777777777776*m9 +
        0.05263157894736842*rho + 0.08333333333333333*(Fx+Fy);
        dist[n+7*Np] = fq;
        
        fq = 0. - 0.1*jx - 0.1*jy + 0.003341687552213868*m1 +
        0.013888888888888888*m10 + 0.0833333333333333333*m11 + 0.041666666666666664*m12 + 0.25*m13 -
        0.125*m16 + 0.125*m17 + 0.003968253968253968*m2 - 0.025*m4 - 0.025*m6 +
        0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx+Fy);
        dist[n+8*Np] = fq;
        
        fq = 0. + 0.1*jx - 0.1*jy + 0.003341687552213868*m1 + 0.013888888888888888*m10 +
        0.0833333333333333333*m11 + 0.041666666666666664*m12 - 0.25*m13 + 0.125*m16 + 0.125*m17 +
        0.003968253968253968*m2 + 0.025*m4 - 0.025*m6 + 0.027777777777777776*m9 +
        0.05263157894736842*rho + 0.08333333333333333*(Fx-Fy);
        dist[n+9*Np] = fq;
        
        fq = 0. - 0.1*jx + 0.1*jy + 0.003341687552213868*m1 +
        0.013888888888888888*m10 + 0.0833333333333333333*m11 + 0.041666666666666664*m12 - 0.25*m13 -
        0.125*m16 - 0.125*m17 + 0.003968253968253968*m2 - 0.025*m4 + 0.025*m6 +
        0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx-Fy);
        dist[n+10*Np] = fq;
        
        fq =   0. + 0.1*jx + 0.1*jz + 0.003341687552213868*m1 + 0.013888888888888888*m10 -
        0.0833333333333333333*m11 - 0.041666666666666664*m12 + 0.25*m15 - 0.125*m16 + 0.125*m18 +
        0.003968253968253968*m2 + 0.025*m4 + 0.025*m8 + 0.027777777777777776*m9 +
        0.05263157894736842*rho + 0.08333333333333333*(Fx+Fz);
        dist[n+11*Np] = fq;
        
        fq = 0. - 0.1*jx - 0.1*jz + 0.003341687552213868*m1 +
        0.013888888888888888*m10 - 0.0833333333333333333*m11 - 0.041666666666666664*m12 + 0.25*m15 +
        0.125*m16 - 0.125*m18 + 0.003968253968253968*m2 - 0.025*m4 - 0.025*m8 +
        0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx+Fz);
        dist[n+12*Np] = fq;
        
        fq =   0. + 0.1*jx - 0.1*jz + 0.003341687552213868*m1 + 0.013888888888888888*m10 -
        0.0833333333333333333*m11 - 0.041666666666666664*m12 - 0.25*m15 - 0.125*m16 - 0.125*m18 +
        0.003968253968253968*m2 + 0.025*m4 - 0.025*m8 + 0.027777777777777776*m9 +
        0.05263157894736842*rho + 0.08333333333333333*(Fx-Fz);
        dist[n+13*Np] = fq;
        
        fq = 0. - 0.1*jx + 0.1*jz + 0.003341687552213868*m1 +
        0.013888888888888888*m10 - 0.0833333333333333333*m11 - 0.041666666666666664*m12 - 0.25*m15 +
        0.125*m16 + 0.125*m18 + 0.003968253968253968*m2 - 0.025*m4 + 0.025*m8 +
        0.027777777777777776*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fx-Fz);
        dist[n+14*Np] = fq;
        
        fq =  0. + 0.1*jy + 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 + 0.25*m14 +
        0.125*m17 - 0.125*m18 + 0.003968253968253968*m2 + 0.025*m6 + 0.025*m8 -
        0.05555555555555555*m9 + 0.05263157894736842*rho + 0.08333333333333333*(Fy+Fz);
        dist[n+15*Np] = fq;
        
        fq = 0. - 0.1*jy - 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 + 0.25*m14 -
        0.125*m17 + 0.125*m18 + 0.003968253968253968*m2 - 0.025*m6 - 0.025*m8 -
        0.05555555555555555*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fy+Fz);
        dist[n+16*Np] = fq;
        
        fq = 0. + 0.1*jy - 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 - 0.25*m14 +
        0.125*m17 + 0.125*m18 + 0.003968253968253968*m2 + 0.025*m6 - 0.025*m8 -
        0.05555555555555555*m9 + 0.05263157894736842*rho + 0.08333333333333333*(Fy-Fz);
        dist[n+17*Np] = fq;
        
        fq = 0. - 0.1*jy + 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 - 0.25*m14 -
        0.125*m17 - 0.125*m18 + 0.003968253968253968*m2 - 0.025*m6 + 0.025*m8 -
        0.05555555555555555*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fy-Fz);
        dist[n+18*Np] = fq;
        
        
        ux = jx/rho0;
        ux = jy/rho0;
        ux = jz/rho0;
        
        nAB = 1.0/(nA+nB);
        Aq[n] = 0.3333333333333333*nA;
        Bq[n] = 0.3333333333333333*nB;
        
        delta = beta*nA*nB*nAB*0.1111111111111111*nx;
        if (!(nA*nB*nAB>0)) delta=0;
        a1 = nA*(0.1111111111111111*(1+4.5*ux))+delta;
        b1 = nB*(0.1111111111111111*(1+4.5*ux))-delta;
        a2 = nA*(0.1111111111111111*(1-4.5*ux))-delta;
        b2 = nB*(0.1111111111111111*(1-4.5*ux))+delta;
        
        Aq[1*Np+n] = a1;
        Bq[1*Np+n] = b1;
        Aq[2*Np+n] = a2;
        Bq[2*Np+n] = b2;
        
        delta = beta*nA*nB*nAB*0.1111111111111111*ny;
        if (!(nA*nB*nAB>0)) delta=0;
        a1 = nA*(0.1111111111111111*(1+4.5*uy))+delta;
        b1 = nB*(0.1111111111111111*(1+4.5*uy))-delta;
        a2 = nA*(0.1111111111111111*(1-4.5*uy))-delta;
        b2 = nB*(0.1111111111111111*(1-4.5*uy))+delta;
        
        Aq[3*Np+n] = a1;
        Bq[3*Np+n] = b1;
        Aq[4*Np+n] = a2;
        Bq[4*Np+n] = b2;

        delta = beta*nA*nB*nAB*0.1111111111111111*nz;
        if (!(nA*nB*nAB>0)) delta=0;
        a1 = nA*(0.1111111111111111*(1+4.5*uz))+delta;
        b1 = nB*(0.1111111111111111*(1+4.5*uz))-delta;
        a2 = nA*(0.1111111111111111*(1-4.5*uz))-delta;
        b2 = nB*(0.1111111111111111*(1-4.5*uz))+delta;
        
        Aq[5*Np+n] = a1;
        Bq[5*Np+n] = b1;
        Aq[6*Np+n] = a2;
        Bq[6*Np+n] = b2;
    }
}


extern "C" void ScaLBL_DensityInterpolation(int *Map, char* ID, double * VFmask, double *PHI, double *DenA, double *DenB, double *laggedDenA, double *laggedDenB, double * analysisDenA, double * analysisDenB, int start, int finish, int Ni, int strideY, int strideZ, double lower_threshold,double upper_threshold) {
    int ijk,n,N,i,j,k,nn;
    // distributions
    double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
    double f10,f11,f12,f13,f14,f15,f16,f17,f18;
    
    double m0,m1,m2,m3,m4,m5,m6,m7,m8,m9;
    double m10,m11,m12,m13,m14,m15,m16,m17,m18;
    double nx,ny,nz;
    
    int count = 0;
    double vfs = 0;
    double new_densityA, new_densityB;
    
    double weight;
    int id;
    double w0 =  1.0;
    double w1to6 = 1.0;
    double w7to19 =  1.0;
    
    double phi;
    
    for (int n=start; n<finish; n++) {
        ijk = Map[n];
        count = 0;
        m0 = 0;
        f0 = w0*laggedDenB[ijk];
        weight = w0;
        count++;
//        weight = 0;
        
        vfs = VFmask[ijk];
//        if (vfs > lower_threshold && vfs <= upper_threshold) {

            nn = ijk-1;
            phi = PHI[nn]; vfs = VFmask[nn];
            if (phi > 0) {
                m0 += phi;
//                m1
                weight += w1to6;
                //std::cout << "!: m1=" << m1 << " f1=" << f1 << std::endl;
                count++;
            }
//            else if (vfs <=upper_threshold) {
//                m1 = w1to6*laggedDenA[nn];
//                f1 = w1to6*laggedDenB[nn];
//                weight += w1to6;
//                //std::cout << "x: m1=" << m1 << " f1=" << f1 << std::endl;
//                count++;
//            }
            else {
                m1 = 0.0;
                f1 = 0.0;
//                m0 += phi;
                count--;
            }
            //........................................................................
            nn = ijk+1;
            phi = PHI[nn]; vfs = VFmask[nn];
            if (phi > 0) {
//                m0 += phi;
//                m2 = w1to6*DenA[nn];
//                f2 = w1to6*DenB[nn];
//                weight += w1to6;
                count++;
            }
//            else if (vfs <=upper_threshold) {
//                m2 = w1to6*laggedDenA[nn];
//                f2 = w1to6*laggedDenB[nn];
//                weight += w1to6;
//                count++;
//            }
            else {
//                m0 += phi;
                m2 = 0.0;
                f2 = 0.0;
                count--;
            }
            //........................................................................
            nn = ijk-strideY;
            phi = PHI[nn]; vfs = VFmask[nn];
            if (phi > 0) {
                m3 = w1to6*DenA[nn];
                f3 = w1to6*DenB[nn];
                weight += w1to6;
                count++;
                m0 += phi;
            }
//            else if (vfs <=upper_threshold) {
//                m3 = w1to6*laggedDenA[nn];
//                f3 = w1to6*laggedDenB[nn];
//                weight += w1to6;
//                count++;
//
//            }
            else {
//                m0 += phi;
                m3 = 0.0;
                f3 = 0.0;
                count--;
            }
            //........................................................................
            nn = ijk+strideY;
            phi = PHI[nn]; vfs = VFmask[nn];
            if (phi > 0) {
                m0 += phi;
                m4 = w1to6*DenA[nn];
                f4 = w1to6*DenB[nn];
                weight += w1to6;
                count++;
            }
//            else if (vfs <=upper_threshold) {
//                m4 = w1to6*laggedDenA[nn];
//                f4 = w1to6*laggedDenB[nn];
//                weight += w1to6;
//                count++;
//            }
            else {
//                m0 += phi;
                m4 = 0.0;
                f4 = 0.0;
                count--;
            }
            //........................................................................
            nn = ijk-strideZ;
            phi = PHI[nn]; vfs = VFmask[nn];
            if (phi > 0) {
                m5 = w1to6*DenA[nn];
                f5 = w1to6*DenB[nn];
                weight += w1to6;
                count++;
                m0 += phi;
            }
//            else if (vfs <=upper_threshold) {
//                m5 = w1to6*laggedDenA[nn];
//                f5 = w1to6*laggedDenB[nn];
//                weight += w1to6;
//                count++;
//            }
            else {
//                m0 += phi;
                count--;
                m5 = 0.0;
                f5 = 0.0;
            }
            //........................................................................
            nn = ijk+strideZ;
            phi = PHI[nn]; vfs = VFmask[nn];
            if (phi > 0) {
                m6 = w1to6*DenA[nn];
                f6 = w1to6*DenB[nn];
                weight += w1to6;
                count++;
                m0 += phi;
            }

            else {
//                m0 += phi;
                count--;
            }


//
//            //........................................................................
//            nn = ijk-strideY-1;
//            id = ID[nn]; vfs = VFmask[nn];
//            if (id > 0) {
//                m7 = w7to19*DenA[nn];
//                f7 = w7to19*DenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else if (vfs <=upper_threshold) {
//                m7 = w7to19*laggedDenA[nn];
//                f7 = w7to19*laggedDenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else {
//                m7 = 0.0;
//                f7 = 0.0;
//            }
//            //........................................................................
//            nn = ijk+strideY+1;
//            id = ID[nn]; vfs = VFmask[nn];
//            if (id > 0) {
//                m8 = w7to19*DenA[nn];
//                f8 = w7to19*DenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else if (vfs <=upper_threshold) {
//                m8 = w7to19*laggedDenA[nn];
//                f8 = w7to19*laggedDenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else {
//                m8 = 0.0;
//                f8 = 0.0;
//            }
//            //........................................................................
//            nn = ijk+strideY-1;
//            id = ID[nn]; vfs = VFmask[nn];
//            if (id > 0) {
//                m9 = w7to19*DenA[nn];
//                f9 = w7to19*DenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else if (vfs <=upper_threshold) {
//                m9 = w7to19*laggedDenA[nn];
//                f9 = w7to19*laggedDenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else {
//                m9 = 0.0;
//                f9 = 0.0;
//            }
//            //........................................................................
//            nn = ijk-strideY+1;
//            id = ID[nn]; vfs = VFmask[nn];
//            if (id > 0) {
//                m10 = w7to19*DenA[nn];
//                f10 = w7to19*DenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else if (vfs <=upper_threshold) {
//                m10 = w7to19*laggedDenA[nn];
//                f10 = w7to19*laggedDenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else {
//                m10 = 0.0;
//                f10 = 0.0;
//            }
//            //........................................................................
//            nn = ijk-strideZ-1;
//            id = ID[nn]; vfs = VFmask[nn];
//            if (id > 0) {
//                m11 = w7to19*DenA[nn];
//                f11 = w7to19*DenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else if (vfs <=upper_threshold) {
//                m11 = w7to19*laggedDenA[nn];
//                f11 = w7to19*laggedDenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else {
//                m11 = 0.0;
//                f11 = 0.0;
//            }
//            //........................................................................
//            nn = ijk+strideZ+1;
//            id = ID[nn]; vfs = VFmask[nn];
//            if (id > 0) {
//                m12 = w7to19*DenA[nn];
//                f12 = w7to19*DenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else if (vfs <=upper_threshold) {
//                m12 = w7to19*laggedDenA[nn];
//                f12 = w7to19*laggedDenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else {
//                m12 = 0.0;
//                f12 = 0.0;
//            }
//            //........................................................................
//            nn = ijk+strideZ-1;
//            id = ID[nn]; vfs = VFmask[nn];
//            if (id > 0) {
//                m13 = w7to19*DenA[nn];
//                f13 = w7to19*DenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else if (vfs <=upper_threshold) {
//                m13 = w7to19*laggedDenA[nn];
//                f13 = w7to19*laggedDenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else {
//                m13 = 0.0;
//                f13 = 0.0;
//            }
//            //........................................................................
//            nn = ijk-strideZ+1;
//            id = ID[nn]; vfs = VFmask[nn];
//            if (id > 0) {
//                m14 = w7to19*DenA[nn];
//                f14 = w7to19*DenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else if (vfs <=upper_threshold) {
//                m14 = w7to19*laggedDenA[nn];
//                f14 = w7to19*laggedDenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else {
//                m14 = 0.0;
//                f14 = 0.0;
//            }
//            //........................................................................
//            nn = ijk-strideZ-strideY;
//            id = ID[nn]; vfs = VFmask[nn];
//            if (id > 0) {
//                m15 = w7to19*DenA[nn];
//                f15 = w7to19*DenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else if (vfs <=upper_threshold) {
//                m15 = w7to19*laggedDenA[nn];
//                f15 = w7to19*laggedDenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else {
//                m15 = 0.0;
//                f15 = 0.0;
//            }
//            //........................................................................
//            nn = ijk+strideZ+strideY;
//            id = ID[nn]; vfs = VFmask[nn];
//            if (id > 0) {
//                m16 = w7to19*DenA[nn];
//                f16 = w7to19*DenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else if (vfs <=upper_threshold) {
//                m16 = w7to19*laggedDenA[nn];
//                f16 = w7to19*laggedDenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else {
//                m16 = 0.0;
//                f16 = 0.0;
//            }
//            //........................................................................
//            nn = ijk+strideZ-strideY;
//            id = ID[nn]; vfs = VFmask[nn];
//            if (id > 0) {
//                m17 = w7to19*DenA[nn];
//                f17 = w7to19*DenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else if (vfs <=upper_threshold) {
//                m17 = w7to19*laggedDenA[nn];
//                f17 = w7to19*laggedDenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else {
//                m17 = 0.0;
//                f17 = 0.0;
//            }
//            //........................................................................
//            nn = ijk-strideZ+strideY;
//            id = ID[nn]; vfs = VFmask[nn];
//            if (id > 0) {
//                m18 = w7to19*DenA[nn];
//                f18 = w7to19*DenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else if (vfs <=upper_threshold) {
//                m18 = w7to19*laggedDenA[nn];
//                f18 = w7to19*laggedDenB[nn];
//                weight += w7to19;
//                count++;
//            }
//            else {
//                m18 = 0.0;
//                f18 = 0.0;
//            }


//            if (count == 0) count = 1;
//
//            if (count <= 0) {new_densityA = m0/2.0; new_densityB = (1.0-new_densityA)/2.0; }
//            else {new_densityA = (1.0-m0)/2.0; new_densityB = (1.0-new_densityA)/2.0;}
            
//            if (count <= 0) {new_densityA = 1.0-m0; new_densityB = m0; }
//            else{ new_densityA = m0; new_densityB = 1.0-m0; }
//            new_densityA = (m1+m2+m3+m4+m5+m6)/double(count);
////                            m7+m8+m9+m10+m11+m12+m13 +
////                            m14+m15+m16+m17+m18)/double(count);
//            new_densityB = (f1+f2+f3+f4+f5+f6)/double(count);
//            f7+f8+f9+f10+f11+f12+f13 +
//                            f14+f15+f16+f17+f18)/double(count);

            
        analysisDenA[ijk] = -0.9;
        analysisDenB[ijk] = -0.9;
            
//        }
        
    }
}



extern "C" void ScaLBL_InitializeDensityInterpolation(int *Map, char * ID, double *phi, double *DenA, double *DenB, double *laggedDenA, double *laggedDenB, double * analysisDenA, double * analysisDenB, int start, int finish, int Ni, int strideY, int strideZ) {
    int ijk,n,N,i,j,k,nn;
    // distributions
    double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
    double f10,f11,f12,f13,f14,f15,f16,f17,f18;
    
    double m0,m1,m2,m3,m4,m5,m6,m7,m8,m9;
    double m10,m11,m12,m13,m14,m15,m16,m17,m18;
    double nx,ny,nz;
    
    int count = 0;
    
    double new_densityA, new_densityB;
    
    for (int n=start; n<finish; n++) {
        ijk = Map[n];
        count = 0;
        m1 = f1 = 0.0;
        

        nn = ijk-1;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk+1;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk-strideY;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk+strideY;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk-strideZ;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk+strideZ;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk-strideY-1;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk+strideY+1;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk+strideY-1;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk-strideY+1;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk-strideZ-1;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk+strideZ+1;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk+strideZ-1;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk-strideZ+1;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk-strideZ-strideY;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk+strideZ+strideY;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk+strideZ-strideY;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }
        //........................................................................
        nn = ijk-strideZ+strideY;
        if (ID[nn] > 0) {
            m1 += DenA[nn];
            f1 += DenB[nn];
            count++;
        }

        
        if (count == 0) count = 1;
        new_densityA = m1/double(count);
        new_densityB = f1/double(count);
        
        if (new_densityA >= 0.5) { new_densityA = 1; new_densityB = 0; }
        else { new_densityB = 1; new_densityA = 0; }
        
        
        laggedDenA[ijk] = new_densityA;
        laggedDenB[ijk] = new_densityB;
        
        analysisDenA[ijk] = new_densityA;
        analysisDenB[ijk] = new_densityB;
        
        
        
    }
}


/* =========================================================================== */



extern "C" void Interpolation(int *interpolationList,int *neighborList, int *Map, double *dist, int start, int finish, int Np, int N, double * LIBBqA, double * LIBBqD) {}

extern "C" void LIBB_TWO_STREAM_TEST(int * interpolationList, int *neighborList, int *Map, double *dist, int start, int finish, int Np, int N, double * LIBBqA, double * LIBBqBC, double * LIBBqD,int rank){}

extern "C" void LIBB_INTERPOLATION_TEST(int *interpolationList,int *neighborList, int *Map, double *dist, int start, int finish, int Np, int N, double * LIBBqA, double * LIBBqBC, double * LIBBqD, int Nx, int Ny, int Nz,int rank,double * phi) {}



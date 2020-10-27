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
#include <math.h>


#include <iostream>

#define STOKES


extern "C" void save_densities(int *Map, double *dist, double* saved_dist, int start, int finish, int Np){
    int ijk;
    for (size_t n = start; n < finish; n++){
        ijk = Map[n];
        saved_dist[ijk] = dist[ijk];
    }
}

extern "C" void save_state(double *dist, double* saved_dist, int start, int finish, int q, int Np){
    for (int a = 0; a < q; a++) {
        for (size_t n = start; n < finish; n++){
            saved_dist[n + a*Np] = dist[n + a*Np];
        }
    }
}


extern "C" void ScaLBL_D3Q7_PhaseField(int *neighborList, int *Map, double *Aq, double *Bq,
                                       double *Den, double *Phi, int start, int finish, int Np){
    
    int idx,n,nread;
    double fq,nA,nB;
    
    for (int n=start; n<finish; n++) {
        //..........Compute the number density for component A............
        
        fq = Aq[n];
        nA = fq;
        
        nread = neighborList[n];
        fq = Aq[nread];
        nA += fq;
        
        nread = neighborList[n+1*Np];
        fq = Aq[nread];
        nA += fq;
        // Autogenerating LBM block...
        nread = neighborList[n+2*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+3*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+4*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+5*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+6*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+7*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+8*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+9*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+10*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+11*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+12*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+13*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+14*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+15*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+16*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+17*Np];
        fq = Aq[nread];
        nA += fq;

        
        
        //..........Compute the number density for component B............
        fq = Bq[n];
        nB = fq;
        
        nread = neighborList[n];
        fq = Bq[nread];
        nB += fq;
        
        nread = neighborList[n+1*Np];
        fq = Bq[nread];
        nB += fq;
        // Autogenerating LBM block...
        nread = neighborList[n+2*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+3*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+4*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+5*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+6*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+7*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+8*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+9*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+10*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+11*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+12*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+13*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+14*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+15*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+16*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+17*Np];
        fq = Bq[nread];
        nB += fq;

        
        
      
        // save the number densities
        Den[n] = nA;
        Den[Np+n] = nB;
        
        // save the phase indicator field
        idx = Map[n];
        Phi[idx] = (nA-nB)/(nA+nB);
    }
}





extern "C" void ScaLBL_D3Q19_Color(int *neighborList, int *Map, double *dist, double* dist2, double *Aq, double *Bq, double *Den,double *Phi, double *Vel, double *Press, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta, double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np){
    
    int n,nn,ijk,nread;
    int nr1,nr2,nr3,nr4,nr5,nr6;
    int nr7,nr8,nr9,nr10;
    int nr11,nr12,nr13,nr14;
    int nr15,nr16,nr17,nr18;
    double fq;
    // conserved momemnts
    double rho,jx,jy,jz;
    double jxp, jyp, jzp;
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
    const double mrt_V6=0.055555555555555556;
    const double mrt_V7=0.02777777777777778;
    const double mrt_V8=0.08333333333333333333333;
    const double mrt_V9=0.003341687552213868;
    const double mrt_V10=0.003968253968253968;
    const double mrt_V11=0.01388888888888889;
    const double mrt_V12=0.04166666666666666;
    
    
  
    
    for (int n=start; n<finish; n++) {
        
        nA = Den[n];
        nB = Den[Np + n];
        
        // compute phase indicator field
        phi=(nA-nB)/(nA+nB);
        
        // local density
        rho0=1.0 + 0.5*(1.0-phi)*(rhoB-rhoA);  // 1.0 for single phase
        // local relaxation time
        tau=tauA + 0.5*(1.0-phi)*(tauB-tauA); // tau = tauA for single phase
        rlx_setA = 1.0/tau;
        rlx_setB = 8.*(2.-rlx_setA)/(8.-rlx_setA); // 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
        
        //           // Get the 1D index based on regular data layout
        ijk = Map[n];
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
        
        //...........Normalize the Color Gradient.................................
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
        
        // q = 1 // Opposite read
        nr1 = neighborList[n];
        fq = dist[nr1];
        rho += fq;
        m1 -= 11.0*fq;
        m2 -= 4.0*fq;
        jx = fq;
        m4 = -4.0*fq;
        m9 = 2.0*fq;
        m10 = -4.0*fq;
        
        // q = 2
        nr2 = neighborList[n+Np];
        fq = dist[nr2];
        rho += fq;
        m1 -= 11.0*(fq);
        m2 -= 4.0*(fq);
        jx -= fq;
        m4 += 4.0*(fq);
        m9 += 2.0*(fq);
        m10 -= 4.0*(fq);
        
        // q = 3
        nr3 = neighborList[n+2*Np];
        fq = dist[nr3];
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
        nr4 = neighborList[n+3*Np];
        fq = dist[nr4];
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
        nr5 = neighborList[n+4*Np];
        fq = dist[nr5];
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
        nr6 = neighborList[n+5*Np];
        fq = dist[nr6];
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
        nr7 = neighborList[n+6*Np];
        fq = dist[nr7];
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
        nr8 = neighborList[n+7*Np];
        fq = dist[nr8];
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
        nr9 = neighborList[n+8*Np];
        fq = dist[nr9];
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
        nr10 = neighborList[n+9*Np];
        fq = dist[nr10];
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
        nr11 = neighborList[n+10*Np];
        fq = dist[nr11];
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
        nr12 = neighborList[n+11*Np];
        fq = dist[nr12];
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
        nr13 = neighborList[n+12*Np];
        fq = dist[nr13];
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
        nr14 = neighborList[n+13*Np];
        fq = dist[nr14];
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
        nr15 = neighborList[n+14*Np];
        fq = dist[nr15];
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
        nr16 = neighborList[n+15*Np];
        fq = dist[nr16];
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
        nr17 = neighborList[n+16*Np];
        fq = dist[nr17];
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
        nr18 = neighborList[n+17*Np];
        fq = dist[nr18];
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
        
        
        Press[n] = (1./3.)*rho;
        Vel[n] = jx + 0.5*Fx;
        Vel[Np+n] = jy + 0.5*Fy;
        Vel[2*Np+n] = jz + 0.5*Fz;
        
        jxp = jx + 0.5*Fx;
        jyp = jy + 0.5*Fy;
        jzp = jz + 0.5*Fz;
        
        
        // ENTER MOMENTUM INTO THE MOMENTS, NOT VELOCITY
        //..........Toelke, Fruediger et. al. 2006................................
        if (C == 0.0)    nx = ny = nz = 0.0;
        //C = nx = ny = nz = 0;
        m1 = m1 + rlx_setA*((19*(jxp*jxp+jyp*jyp+jzp*jzp)/rho0 - 11*rho) -alpha*C - m1);
        m2 = m2 + rlx_setA*((3*rho - 5.5*(jxp*jxp+jyp*jyp+jzp*jzp)/rho0)- m2);
        m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
        m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
        m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
        // pxx
        m9 = m9 + rlx_setA*(   ((2*jxp*jxp-jyp*jyp-jzp*jzp)/(3.0*rho0)) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
        // PIxx
        m10 = m10 + rlx_setA*( -0.5*((2*jxp*jxp-jyp*jyp-jzp*jzp)/(3.0*rho0))  - m10);
        // pww
        m11 = m11 + rlx_setA*(   ((jyp*jyp-jzp*jzp)/rho0) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
        // PIww
        m12 = m12 + rlx_setA*( -0.5*((jyp*jyp-jzp*jzp)/rho0)  - m12);
        m13 = m13 + rlx_setA*( (jxp*jyp/rho0) + 0.5*alpha*C*nx*ny - m13);
        m14 = m14 + rlx_setA*( (jyp*jzp/rho0) + 0.5*alpha*C*ny*nz - m14);
        m15 = m15 + rlx_setA*( (jxp*jzp/rho0) + 0.5*alpha*C*nx*nz - m15);
        m16 = m16 + rlx_setB*( - m16);
        m17 = m17 + rlx_setB*( - m17);
        m18 = m18 + rlx_setB*( - m18);
        
        /* Stokes flow */
        //                     C = nx = ny = nz = 0;
        //                     m1 = m1 + rlx_setA*(- 11*rho - m1);
        //                     m2 = m2 + rlx_setA*(3*rho - m2);
        //                     m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
        //                     m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
        //                     m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
        //                     // pxx
        //                     m9 = m9 + rlx_setA*( - m9);
        //                     // PIxx
        //                     m10 = m10 + rlx_setA*( - m10);
        //                     // pww
        //                     m11 = m11 + rlx_setA*( - m11);
        //                     // PIww
        //                     m12 = m12 + rlx_setA*( - m12);
        //                     m13 = m13 + rlx_setA*( - m13);
        //                     m14 = m14 + rlx_setA*( - m14);
        //                     m15 = m15 + rlx_setA*( - m15);
        //                     m16 = m16 + rlx_setB*( - m16);
        //                     m17 = m17 + rlx_setB*( - m17);
        //                     m18 = m18 + rlx_setB*( - m18);
        //.................inverse transformation......................................................
        
        fq = 0. - 0.012531328320802004*m1 + 0.047619047619047616*m2 + 0.05263157894736842*rho;
        dist2[n] = fq;  // 0
        
        fq = 0. + 0.1*jx - 0.004594820384294068*m1 - 0.05555555555555555*m10 - 0.015873015873015872*m2 -
        0.1*m4 + 0.05555555555555555*m9 + 0.05263157894736842*rho + 0.16666666666666667*Fx;
        dist2[n+Np] = fq; // 1
        
        fq =  0. - 0.1*jx - 0.004594820384294068*m1 - 0.05555555555555555*m10 - 0.015873015873015872*m2 +
        0.1*m4 + 0.05555555555555555*m9 + 0.05263157894736842*rho - 0.16666666666666667*Fx;
        dist2[n+2*Np] = fq; // 2
        
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
        
        fq =  0. - 0.1*jy - 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 + 0.25*m14 -
        0.125*m17 + 0.125*m18 + 0.003968253968253968*m2 - 0.025*m6 - 0.025*m8 -
        0.05555555555555555*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fy+Fz);
        dist2[n+16*Np] = fq;
        
        fq =   0. + 0.1*jy - 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 - 0.25*m14 +
        0.125*m17 + 0.125*m18 + 0.003968253968253968*m2 + 0.025*m6 - 0.025*m8 -
        0.05555555555555555*m9 + 0.05263157894736842*rho + 0.08333333333333333*(Fy-Fz);
        dist2[n+17*Np] = fq;
        
        fq =   0. - 0.1*jy + 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 - 0.25*m14 -
        0.125*m17 - 0.125*m18 + 0.003968253968253968*m2 - 0.025*m6 + 0.025*m8 -
        0.05555555555555555*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fy-Fz);
        dist2[n+18*Np] = fq;
        
        
        ux = jx / rho0;
        uy = jy / rho0;
        uz = jz / rho0;

        delta = 0;
        Aq[n] = nA/3.;
        Bq[n] = nB/3.;
        nAB = 1./(nA+nB);
        
        delta = (beta*nA*nB*nx)*nAB;
        if (!(nA*nB*nAB>0)) delta=0;
        Aq[n+1*Np] = (delta  + nA*(1  + 1*ux))/18.;
        Bq[n+1*Np] = (-delta + nB*(1  + 1*ux))/18.;
        Aq[n+2*Np] = (-delta + nA*(1  - 1*ux))/18.;
        Bq[n+2*Np] = (delta  + nB*(1  - 1*ux))/18.;
        
        delta = (beta*nA*nB*ny)*nAB;
        if (!(nA*nB*nAB>0)) delta=0;
        Aq[n+3*Np] = (delta  + nA*(1  + 1*uy))/18.;
        Bq[n+3*Np] = (-delta + nB*(1  + 1*uy))/18.;
        Aq[n+4*Np] = (-delta + nA*(1  - 1*uy))/18.;
        Bq[n+4*Np] = (delta  + nB*(1  - 1*uy))/18.;
        
        delta = (beta*nA*nB*nz)*nAB;
        if (!(nA*nB*nAB>0)) delta=0;
        Aq[n+5*Np] = (delta  + nA*(1  + 1*uz))/18.;
        Bq[n+5*Np] = (-delta + nB*(1  + 1*uz))/18.;
        Aq[n+6*Np] = (-delta + nA*(1  - 1*uz))/18.;
        Bq[n+6*Np] = (delta  + nB*(1  - 1*uz))/18.;
        
        delta = (beta*nA*nB*(nx + ny))*nAB;
        if (!(nA*nB*nAB>0)) delta=0;
        Aq[n+7*Np] = (delta  + nA*(1  + 1*(ux + uy)))/36.;
        Bq[n+7*Np] = (-delta + nB*(1  + 1*(ux + uy)))/36.;
        Aq[n+8*Np] = (-delta + nA*(1  + 1*(-ux - uy)))/36.;
        Bq[n+8*Np] = (delta  + nB*(1  + 1*(-ux - uy)))/36.;
        
        delta = (beta*nA*nB*(nx - ny))*nAB;
        if (!(nA*nB*nAB>0)) delta=0;
        Aq[n+9*Np]  = (delta  + nA*(1  + 1*(ux - uy)))/36.;
        Bq[n+9*Np]  = (-delta + nB*(1  + 1*(ux - uy)))/36.;
        Aq[n+10*Np] = (-delta + nA*(1  + 1*(-ux + uy)))/36.;
        Bq[n+10*Np] = (delta  + nB*(1  + 1*(-ux + uy)))/36.;
        
        delta = (beta*nA*nB*(nx + nz))*nAB;
        if (!(nA*nB*nAB>0)) delta=0;
        Aq[n+11*Np] = (delta  + nA*(1  + 1*(ux + uz)))/36.;
        Bq[n+11*Np] = (-delta + nB*(1  + 1*(ux + uz)))/36.;
        Aq[n+12*Np] = (-delta + nA*(1  + 1*(-ux - uz)))/36.;
        Bq[n+12*Np] = (delta  + nB*(1  + 1*(-ux - uz)))/36.;
        
        
        delta = (beta*nA*nB*(nx - nz))*nAB;
        if (!(nA*nB*nAB>0)) delta=0;
        Aq[n+13*Np] = (delta  + nA*(1  + 1*(ux - uz)))/36.;
        Bq[n+13*Np] = (-delta + nB*(1  + 1*(ux - uz)))/36.;
        Aq[n+14*Np] = (-delta + nA*(1  + 1*(-ux + uz)))/36.;
        Bq[n+14*Np] = (delta  + nB*(1  + 1*(-ux + uz)))/36.;
        
        
        delta = (beta*nA*nB*(ny + nz))*nAB;
        if (!(nA*nB*nAB>0)) delta=0;
        Aq[n+15*Np] = (delta  + nA*(1  + 1*(uy + uz)))/36.;
        Bq[n+15*Np] = (-delta + nB*(1  + 1*(uy + uz)))/36.;
        Aq[n+16*Np] = (-delta + nA*(1  + 1*(-uy - uz)))/36.;
        Bq[n+16*Np] = (delta  + nB*(1  + 1*(-uy - uz)))/36.;
        
        
        delta = (beta*nA*nB*(ny - nz))*nAB;
        if (!(nA*nB*nAB>0)) delta=0;
        Aq[n+17*Np] = (delta  + nA*(1  + 1*(uy - uz)))/36.;
        Bq[n+17*Np] = (-delta + nB*(1  + 1*(uy - uz)))/36.;
        Aq[n+18*Np] = (-delta + nA*(1  + 1*(-uy + uz)))/36.;
        Bq[n+18*Np] = (delta  + nB*(1  + 1*(-uy + uz)))/36.;
    }
}





extern "C" void ScaLBL_Color_Init(char *ID, double *Den, double *Phi, double das, double dbs, int Nx, int Ny, int Nz)
{
    int n,N;
    
    N = Nx*Ny*Nz;
    
    for (n=0; n<N; n++){
        
        if ( ID[n] == 1){
            Den[n] = 1.0;
            Den[N+n] = 0.0;
            Phi[n] = 1.0;
        }
        else if ( ID[n] == 2){
            Den[n] = 0.0;
            Den[N+n] = 1.0;
            Phi[n] = -1.0;
        }
        else{
            Den[n] = das;
            Den[N+n] = dbs;
            Phi[n] = (das-dbs)/(das+dbs);
        }
    }
}


extern "C" void ScaLBL_Color_InitDistance(char *ID, double *Den, double *Phi, double *Distance,
                                          double das, double dbs, double beta, double xp, int Nx, int Ny, int Nz)
{
    int n,N;
    double d;
    
    N = Nx*Ny*Nz;
    
    for (n=0; n<N; n++){
        
        if ( ID[n] == 1){
            Den[n] = 1.0;
            Den[N+n] = 0.0;
            Phi[n] = 1.0;
        }
        else if ( ID[n] == 2){
            Den[n] = 0.0;
            Den[N+n] = 1.0;
            Phi[n] = -1.0;
        }
        else{
            Den[n] = das;
            Den[N+n] = dbs;
            Phi[n] = (das-dbs)/(das+dbs);
            d = fabs(Distance[n]);
            Phi[n] = (2.f*(exp(-2.f*beta*(d+xp)))/(1.f+exp(-2.f*beta*(d+xp))) - 1.f);
        }
    }
}

//
//extern "C" void ScaLBL_Color_InitDistance(char *ID, double *Den, double *Phi, double *Distance,
//                                          double das, double dbs, double beta, double xp, int Nx, int Ny, int Nz)
//{
//    int n,N;
//    double d;
//
//    N = Nx*Ny*Nz;
//
//    for (n=0; n<N; n++){
//
//        if ( ID[n] == 1){
//            Den[n] = 1.0;
//            Den[N+n] = 0.0;
//            Phi[n] = 1.0;
//        }
//        else if ( ID[n] == 2){
//            Den[n] = 0.0;
//            Den[N+n] = 1.0;
//            Phi[n] = -1.0;
//        }
//        else{
//            Den[n] = das;
//            Den[N+n] = dbs;
//            Phi[n] = (das-dbs)/(das+dbs);
//            d = fabs(Distance[n]);
//            Phi[n] = (2.f*(exp(-2.f*beta*(d+xp)))/(1.f+exp(-2.f*beta*(d+xp))) - 1.f);
//        }
//    }
//}
//


//*************************************************************************

//*************************************************************************
extern "C" void ScaLBL_Color_BC(int *list, int *Map, double *Phi, double *Den, double vA, double vB, int count, int Np)
{
    int idx,n,nm;
    // Fill the outlet with component b
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        Den[n] = vA;
        Den[Np+n] = vB;
        
        nm = Map[n];
        Phi[nm] = (vA-vB)/(vA+vB);
    }
}

extern "C" void ScaLBL_Color_BC_z(int *list, int *Map, double *Phi, double *Den, double vA, double vB, int count, int Np)
{
    int idx,n,nm;
    // Fill the outlet with component b
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        Den[n] = vA;
        Den[Np+n] = vB;
        //double valB = Den[Np+n]; // mass that reaches inlet is conserved
        
        nm = Map[n];
        Phi[nm] = (vA-vB)/(vA+vB);
    }
}

extern "C" void ScaLBL_Color_BC_Z(int *list, int *Map, double *Phi, double *Den, double vA, double vB, int count, int Np)
{
    int idx,n,nm;
    // Fill the outlet with component b
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        Den[n] = vA;
        Den[Np+n] = vB;
        
        nm = Map[n];
        Phi[nm] = (vA-vB)/(vA+vB);
    }
}

extern "C" void ScaLBL_Color_BC_z_2(int *list, int *Map, double *Phi, double *DenA, double* DenB, double vA, double vB, int count, int Np)
{
    int idx,n,nm;
    // Fill the outlet with component b
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        DenA[n] = vA;
        DenB[n] = vB;
        //double valB = Den[Np+n]; // mass that reaches inlet is conserved
        
        nm = Map[n];
        Phi[nm] = (vA-vB)/(vA+vB);
    }
}

extern "C" void ScaLBL_Color_BC_Z_2(int *list, int *Map, double *Phi, double *DenA, double * DenB, double vA, double vB, int count, int Np)
{
    int idx,n,nm;
    // Fill the outlet with component b
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        DenA[n] = vA;
        DenB[n] = vB;
        
        nm = Map[n];
        Phi[nm] = (vA-vB)/(vA+vB);
    }
}
//*************************************************************************

//*************************************************************************
extern "C" void ScaLBL_D3Q19_ColorGradient(char *ID, double *phi, double *ColorGrad, int Nx, int Ny, int Nz)
{
    int n,N,i,j,k,nn;
    // distributions
    double f1,f2,f3,f4,f5,f6,f7,f8,f9;
    double f10,f11,f12,f13,f14,f15,f16,f17,f18;
    double nx,ny,nz;
    
    // non-conserved moments
    // additional variables needed for computations
    
    N = Nx*Ny*Nz;
    
    for (n=0; n<N; n++){
        
        //.......Back out the 3-D indices for node n..............
        k = n/(Nx*Ny);
        j = (n-Nx*Ny*k)/Nx;
        i = n-Nx*Ny*k-Nx*j;
        //........................................................................
        //........Get 1-D index for this thread....................
        //		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
        //........................................................................
        //					COMPUTE THE COLOR GRADIENT
        //........................................................................
        //.................Read Phase Indicator Values............................
        //........................................................................
        nn = n-1;							// neighbor index (get convention)
        if (i-1<0)		nn += Nx;			// periodic BC along the x-boundary
        f1 = phi[nn];						// get neighbor for phi - 1
        //........................................................................
        nn = n+1;							// neighbor index (get convention)
        if (!(i+1<Nx))	nn -= Nx;			// periodic BC along the x-boundary
        f2 = phi[nn];						// get neighbor for phi - 2
        //........................................................................
        nn = n-Nx;							// neighbor index (get convention)
        if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
        f3 = phi[nn];					// get neighbor for phi - 3
        //........................................................................
        nn = n+Nx;							// neighbor index (get convention)
        if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
        f4 = phi[nn];					// get neighbor for phi - 4
        //........................................................................
        nn = n-Nx*Ny;						// neighbor index (get convention)
        if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
        f5 = phi[nn];					// get neighbor for phi - 5
        //........................................................................
        nn = n+Nx*Ny;						// neighbor index (get convention)
        if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
        f6 = phi[nn];					// get neighbor for phi - 6
        //........................................................................
        nn = n-Nx-1;						// neighbor index (get convention)
        if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
        if (j-1<0)			nn += Nx*Ny;	// Perioidic BC along the y-boundary
        f7 = phi[nn];					// get neighbor for phi - 7
        //........................................................................
        nn = n+Nx+1;						// neighbor index (get convention)
        if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
        if (!(j+1<Ny))		nn -= Nx*Ny;	// Perioidic BC along the y-boundary
        f8 = phi[nn];					// get neighbor for phi - 8
        //........................................................................
        nn = n+Nx-1;						// neighbor index (get convention)
        if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
        if (!(j+1<Ny))		nn -= Nx*Ny;	// Perioidic BC along the y-boundary
        f9 = phi[nn];					// get neighbor for phi - 9
        //........................................................................
        nn = n-Nx+1;						// neighbor index (get convention)
        if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
        if (j-1<0)			nn += Nx*Ny;	// Perioidic BC along the y-boundary
        f10 = phi[nn];					// get neighbor for phi - 10
        //........................................................................
        nn = n-Nx*Ny-1;						// neighbor index (get convention)
        if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
        if (k-1<0)			nn += Nx*Ny*Nz;	// Perioidic BC along the z-boundary
        f11 = phi[nn];					// get neighbor for phi - 11
        //........................................................................
        nn = n+Nx*Ny+1;						// neighbor index (get convention)
        if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
        if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
        f12 = phi[nn];					// get neighbor for phi - 12
        //........................................................................
        nn = n+Nx*Ny-1;						// neighbor index (get convention)
        if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
        if (!(k+1<Nz))		nn -= Nx*Ny*Nz;	// Perioidic BC along the z-boundary
        f13 = phi[nn];					// get neighbor for phi - 13
        //........................................................................
        nn = n-Nx*Ny+1;						// neighbor index (get convention)
        if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
        if (k-1<0)			nn += Nx*Ny*Nz;	// Perioidic BC along the z-boundary
        f14 = phi[nn];					// get neighbor for phi - 14
        //........................................................................
        nn = n-Nx*Ny-Nx;					// neighbor index (get convention)
        if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
        if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
        f15 = phi[nn];					// get neighbor for phi - 15
        //........................................................................
        nn = n+Nx*Ny+Nx;					// neighbor index (get convention)
        if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
        if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
        f16 = phi[nn];					// get neighbor for phi - 16
        //........................................................................
        nn = n+Nx*Ny-Nx;					// neighbor index (get convention)
        if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
        if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
        f17 = phi[nn];					// get neighbor for phi - 17
        //........................................................................
        nn = n-Nx*Ny+Nx;					// neighbor index (get convention)
        if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
        if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
        f18 = phi[nn];					// get neighbor for phi - 18
        //............Compute the Color Gradient...................................
        nx = -(f1-f2+0.5*(f7-f8+f9-f10+f11-f12+f13-f14));
        ny = -(f3-f4+0.5*(f7-f8-f9+f10+f15-f16+f17-f18));
        nz = -(f5-f6+0.5*(f11-f12-f13+f14+f15-f16-f17+f18));
        //...........Normalize the Color Gradient.................................
        //	C = sqrt(nx*nx+ny*ny+nz*nz);
        //	nx = nx/C;
        //	ny = ny/C;
        //	nz = nz/C;
        //...Store the Color Gradient....................
        ColorGrad[n] = nx;
        ColorGrad[N+n] = ny;
        ColorGrad[2*N+n] = nz;
        //...............................................
    }
}
//*************************************************************************
extern "C" void ColorCollide( char *ID, double *disteven, double *distodd, double *ColorGrad,
                             double *Velocity, int Nx, int Ny, int Nz, double rlx_setA, double rlx_setB,
                             double alpha, double beta, double Fx, double Fy, double Fz, bool pBC)
{
    
    int n,N;
    // distributions
    double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
    double f10,f11,f12,f13,f14,f15,f16,f17,f18;
    
    // non-conserved moments
    double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
    // additional variables needed for computations
    double rho,jx,jy,jz,C,nx,ny,nz;
    
    N = Nx*Ny*Nz;
    char id;
    
    for (n=0; n<N; n++){
        
        id = ID[n];
        
        if (id > 0){
            
            // Retrieve the color gradient
            nx = ColorGrad[n];
            ny = ColorGrad[N+n];
            nz = ColorGrad[2*N+n];
            //...........Normalize the Color Gradient.................................
            C = sqrt(nx*nx+ny*ny+nz*nz);
            if (C==0.0) C=1.0;
            nx = nx/C;
            ny = ny/C;
            nz = nz/C;
            //......No color gradient at z-boundary if pressure BC are set.............
            //	if (pBC && k==0) nx = ny = nz = 0.f;
            //	if (pBC && k==Nz-1) nx = ny = nz = 0.f;
            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
            f2 = distodd[n];
            f4 = distodd[N+n];
            f6 = distodd[2*N+n];
            f8 = distodd[3*N+n];
            f10 = distodd[4*N+n];
            f12 = distodd[5*N+n];
            f14 = distodd[6*N+n];
            f16 = distodd[7*N+n];
            f18 = distodd[8*N+n];
            //........................................................................
            f0 = disteven[n];
            f1 = disteven[N+n];
            f3 = disteven[2*N+n];
            f5 = disteven[3*N+n];
            f7 = disteven[4*N+n];
            f9 = disteven[5*N+n];
            f11 = disteven[6*N+n];
            f13 = disteven[7*N+n];
            f15 = disteven[8*N+n];
            f17 = disteven[9*N+n];
            //........................................................................
            //					PERFORM RELAXATION PROCESS
            //........................................................................
            //....................compute the moments...............................................
            rho = f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;
            m1 = -30*f0-11*(f2+f1+f4+f3+f6+f5)+8*(f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18 +f17);
            m2 = 12*f0-4*(f2+f1 +f4+f3+f6 +f5)+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;
            jx = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
            m4 = 4*(-f1+f2)+f7-f8+f9-f10+f11-f12+f13-f14;
            jy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
            m6 = -4*(f3-f4)+f7-f8-f9+f10+f15-f16+f17-f18;
            jz = f5-f6+f11-f12-f13+f14+f15-f16-f17+f18;
            m8 = -4*(f5-f6)+f11-f12-f13+f14+f15-f16-f17+f18;
            m9 = 2*(f1+f2)-f3-f4-f5-f6+f7+f8+f9+f10+f11+f12+f13+f14-2*(f15+f16+f17+f18);
            m10 = -4*(f1+f2)+2*(f4+f3+f6+f5)+f8+f7+f10+f9+f12+f11+f14+f13-2*(f16+f15+f18+f17);
            m11 = f4+f3-f6-f5+f8+f7+f10+f9-f12-f11-f14-f13;
            m12 = -2*(f4+f3-f6-f5)+f8+f7+f10+f9-f12-f11-f14-f13;
            m13 = f8+f7-f10-f9;
            m14 = f16+f15-f18-f17;
            m15 = f12+f11-f14-f13;
            m16 = f7-f8+f9-f10-f11+f12-f13+f14;
            m17 = -f7+f8+f9-f10+f15-f16+f17-f18;
            m18 = f11-f12-f13+f14-f15+f16+f17-f18;
            //..........Toelke, Fruediger et. al. 2006...............
            if (C == 0.0)	nx = ny = nz = 1.0;
#ifdef STOKES
            m1 = m1 + rlx_setA*(- 11*rho -alpha*C - m1);
            m2 = m2 + rlx_setA*(3*rho - m2);
            m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
            m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
            m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
            m9 = m9 + rlx_setA*( 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
            m10 = m10 + rlx_setA*( - m10);
            m11 = m11 + rlx_setA*( 0.5*alpha*C*(ny*ny-nz*nz)- m11);
            m12 = m12 + rlx_setA*( - m12);
            m13 = m13 + rlx_setA*(  0.5*alpha*C*nx*ny - m13);
            m14 = m14 + rlx_setA*(  0.5*alpha*C*ny*nz - m14);
            m15 = m15 + rlx_setA*(  0.5*alpha*C*nx*nz - m15);
            m16 = m16 + rlx_setB*( - m16);
            m17 = m17 + rlx_setB*( - m17);
            m18 = m18 + rlx_setB*( - m18);
#else
            m1 = m1 + rlx_setA*((19*(jx*jx+jy*jy+jz*jz)/rho - 11*rho) -alpha*C - m1);
            m2 = m2 + rlx_setA*((3*rho - 5.5*(jx*jx+jy*jy+jz*jz)/rho)- m2);
            m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
            m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
            m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
            m9 = m9 + rlx_setA*(((2*jx*jx-jy*jy-jz*jz)/rho) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
            m10 = m10 + rlx_setA*( - m10);
            m11 = m11 + rlx_setA*(((jy*jy-jz*jz)/rho) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
            m12 = m12 + rlx_setA*( - m12);
            m13 = m13 + rlx_setA*( (jx*jy/rho) + 0.5*alpha*C*nx*ny - m13);
            m14 = m14 + rlx_setA*( (jy*jz/rho) + 0.5*alpha*C*ny*nz - m14);
            m15 = m15 + rlx_setA*( (jx*jz/rho) + 0.5*alpha*C*nx*nz - m15);
            m16 = m16 + rlx_setB*( - m16);
            m17 = m17 + rlx_setB*( - m17);
            m18 = m18 + rlx_setB*( - m18);
#endif
            //.................inverse transformation......................................................
            f0 = 0.05263157894736842*rho-0.012531328320802*m1+0.04761904761904762*m2;
            f1 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
            +0.1*(jx-m4)+0.0555555555555555555555555*(m9-m10);
            f2 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
            +0.1*(m4-jx)+0.0555555555555555555555555*(m9-m10);
            f3 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
            +0.1*(jy-m6)+0.02777777777777778*(m10-m9)+0.08333333333333333333333*(m11-m12);
            f4 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
            +0.1*(m6-jy)+0.02777777777777778*(m10-m9)+0.08333333333333333333333*(m11-m12);
            f5 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
            +0.1*(jz-m8)+0.02777777777777778*(m10-m9)+0.08333333333333333333333*(m12-m11);
            f6 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
            +0.1*(m8-jz)+0.02777777777777778*(m10-m9)+0.08333333333333333333333*(m12-m11);
            f7 = 0.05263157894736842*rho+0.003341687552213868*m1+0.003968253968253968*m2+0.1*(jx+jy)+0.025*(m4+m6)
            +0.02777777777777778*m9+0.01388888888888889*m10+0.08333333333333333333333*m11
            +0.04166666666666666*m12+0.25*m13+0.125*(m16-m17);
            f8 = 0.05263157894736842*rho+0.003341687552213868*m1+0.003968253968253968*m2-0.1*(jx+jy)-0.025*(m4+m6)
            +0.02777777777777778*m9+0.01388888888888889*m10+0.08333333333333333333333*m11
            +0.04166666666666666*m12+0.25*m13+0.125*(m17-m16);
            f9 = 0.05263157894736842*rho+0.003341687552213868*m1+0.003968253968253968*m2+0.1*(jx-jy)+0.025*(m4-m6)
            +0.02777777777777778*m9+0.01388888888888889*m10+0.08333333333333333333333*m11
            +0.04166666666666666*m12-0.25*m13+0.125*(m16+m17);
            f10 = 0.05263157894736842*rho+0.003341687552213868*m1+0.003968253968253968*m2+0.1*(jy-jx)+0.025*(m6-m4)
            +0.02777777777777778*m9+0.01388888888888889*m10+0.08333333333333333333333*m11
            +0.04166666666666666*m12-0.25*m13-0.125*(m16+m17);
            f11 = 0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2+0.1*(jx+jz)+0.025*(m4+m8)
            +0.02777777777777778*m9+0.01388888888888889*m10-0.08333333333333333333333*m11
            -0.04166666666666666*m12+0.25*m15+0.125*(m18-m16);
            f12 = 0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2-0.1*(jx+jz)-0.025*(m4+m8)
            +0.02777777777777778*m9+0.01388888888888889*m10-0.08333333333333333333333*m11
            -0.04166666666666666*m12+0.25*m15+0.125*(m16-m18);
            f13 = 0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2+0.1*(jx-jz)+0.025*(m4-m8)
            +0.02777777777777778*m9+0.01388888888888889*m10-0.08333333333333333333333*m11
            -0.04166666666666666*m12-0.25*m15-0.125*(m16+m18);
            f14 = 0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2+0.1*(jz-jx)+0.025*(m8-m4)
            +0.02777777777777778*m9+0.01388888888888889*m10-0.08333333333333333333333*m11
            -0.04166666666666666*m12-0.25*m15+0.125*(m16+m18);
            f15 = 0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2+0.1*(jy+jz)+0.025*(m6+m8)
            -0.0555555555555555555555555*m9-0.02777777777777778*m10+0.25*m14+0.125*(m17-m18);
            f16 =  0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2-0.1*(jy+jz)-0.025*(m6+m8)
            -0.0555555555555555555555555*m9-0.02777777777777778*m10+0.25*m14+0.125*(m18-m17);
            f17 = 0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2+0.1*(jy-jz)+0.025*(m6-m8)
            -0.0555555555555555555555555*m9-0.02777777777777778*m10-0.25*m14+0.125*(m17+m18);
            f18 = 0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2+0.1*(jz-jy)+0.025*(m8-m6)
            -0.0555555555555555555555555*m9-0.02777777777777778*m10-0.25*m14-0.125*(m17+m18);
            //.......................................................................................................
            // incorporate external force
            f1 += 0.16666666666666667*Fx;
            f2 -= 0.16666666666666667*Fx;
            f3 += 0.16666666666666667*Fy;
            f4 -= 0.16666666666666667*Fy;
            f5 += 0.16666666666666667*Fz;
            f6 -= 0.16666666666666667*Fz;
            f7 += 0.08333333333333333*(Fx+Fy);
            f8 -= 0.08333333333333333*(Fx+Fy);
            f9 += 0.08333333333333333*(Fx-Fy);
            f10 -= 0.08333333333333333*(Fx-Fy);
            f11 += 0.08333333333333333*(Fx+Fz);
            f12 -= 0.08333333333333333*(Fx+Fz);
            f13 += 0.08333333333333333*(Fx-Fz);
            f14 -= 0.08333333333333333*(Fx-Fz);
            f15 += 0.08333333333333333*(Fy+Fz);
            f16 -= 0.08333333333333333*(Fy+Fz);
            f17 += 0.08333333333333333*(Fy-Fz);
            f18 -= 0.08333333333333333*(Fy-Fz);
            //*********** WRITE UPDATED VALUES TO MEMORY ******************
            // Write the updated distributions
            //....EVEN.....................................
            disteven[n] = f0;
            disteven[N+n] = f2;
            disteven[2*N+n] = f4;
            disteven[3*N+n] = f6;
            disteven[4*N+n] = f8;
            disteven[5*N+n] = f10;
            disteven[6*N+n] = f12;
            disteven[7*N+n] = f14;
            disteven[8*N+n] = f16;
            disteven[9*N+n] = f18;
            //....ODD......................................
            distodd[n] = f1;
            distodd[N+n] = f3;
            distodd[2*N+n] = f5;
            distodd[3*N+n] = f7;
            distodd[4*N+n] = f9;
            distodd[5*N+n] = f11;
            distodd[6*N+n] = f13;
            distodd[7*N+n] = f15;
            distodd[8*N+n] = f17;
            
            //...Store the Velocity..........................
            Velocity[n] = jx;
            Velocity[N+n] = jy;
            Velocity[2*N+n] = jz;
            /*	Velocity[3*n] = jx;
             Velocity[3*n+1] = jy;
             Velocity[3*n+2] = jz;
             */	//...Store the Color Gradient....................
            //			ColorGrad[3*n] = nx*C;
            //			ColorGrad[3*n+1] = ny*C;
            //			ColorGrad[3*n+2] = nz*C;
            //...............................................
            //***************************************************************
        }	// check if n is in the solid
    } // loop over n
}

extern "C" void ScaLBL_D3Q19_ColorCollide( char *ID, double *disteven, double *distodd, double *phi, double *ColorGrad,
                                          double *Velocity, int Nx, int Ny, int Nz, double rlx_setA, double rlx_setB,
                                          double alpha, double beta, double Fx, double Fy, double Fz)
{
    
    int i,j,k,n,nn,N;
    // distributions
    double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
    double f10,f11,f12,f13,f14,f15,f16,f17,f18;
    
    // non-conserved moments
    double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
    // additional variables needed for computations
    double rho,jx,jy,jz,C,nx,ny,nz;
    
    N = Nx*Ny*Nz;
    char id;
    
    for (n=0; n<N; n++){
        
        id = ID[n];
        
        if (id > 0){
            
            //.......Back out the 3-D indices for node n..............
            k = n/(Nx*Ny);
            j = (n-Nx*Ny*k)/Nx;
            i = n-Nx*Ny*k-Nx*j;
            //........................................................................
            //........Get 1-D index for this thread....................
            //		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
            //........................................................................
            //					COMPUTE THE COLOR GRADIENT
            //........................................................................
            //.................Read Phase Indicator Values............................
            //........................................................................
            nn = n-1;							// neighbor index (get convention)
            if (i-1<0)		nn += Nx;			// periodic BC along the x-boundary
            f1 = phi[nn];						// get neighbor for phi - 1
            //........................................................................
            nn = n+1;							// neighbor index (get convention)
            if (!(i+1<Nx))	nn -= Nx;			// periodic BC along the x-boundary
            f2 = phi[nn];						// get neighbor for phi - 2
            //........................................................................
            nn = n-Nx;							// neighbor index (get convention)
            if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
            f3 = phi[nn];					// get neighbor for phi - 3
            //........................................................................
            nn = n+Nx;							// neighbor index (get convention)
            if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
            f4 = phi[nn];					// get neighbor for phi - 4
            //........................................................................
            nn = n-Nx*Ny;						// neighbor index (get convention)
            if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
            f5 = phi[nn];					// get neighbor for phi - 5
            //........................................................................
            nn = n+Nx*Ny;						// neighbor index (get convention)
            if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
            f6 = phi[nn];					// get neighbor for phi - 6
            //........................................................................
            nn = n-Nx-1;						// neighbor index (get convention)
            if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
            if (j-1<0)			nn += Nx*Ny;	// Perioidic BC along the y-boundary
            f7 = phi[nn];					// get neighbor for phi - 7
            //........................................................................
            nn = n+Nx+1;						// neighbor index (get convention)
            if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
            if (!(j+1<Ny))		nn -= Nx*Ny;	// Perioidic BC along the y-boundary
            f8 = phi[nn];					// get neighbor for phi - 8
            //........................................................................
            nn = n+Nx-1;						// neighbor index (get convention)
            if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
            if (!(j+1<Ny))		nn -= Nx*Ny;	// Perioidic BC along the y-boundary
            f9 = phi[nn];					// get neighbor for phi - 9
            //........................................................................
            nn = n-Nx+1;						// neighbor index (get convention)
            if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
            if (j-1<0)			nn += Nx*Ny;	// Perioidic BC along the y-boundary
            f10 = phi[nn];					// get neighbor for phi - 10
            //........................................................................
            nn = n-Nx*Ny-1;						// neighbor index (get convention)
            if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
            if (k-1<0)			nn += Nx*Ny*Nz;	// Perioidic BC along the z-boundary
            f11 = phi[nn];					// get neighbor for phi - 11
            //........................................................................
            nn = n+Nx*Ny+1;						// neighbor index (get convention)
            if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
            if (!(k+1<Nz))		nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
            f12 = phi[nn];					// get neighbor for phi - 12
            //........................................................................
            nn = n+Nx*Ny-1;						// neighbor index (get convention)
            if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
            if (!(k+1<Nz))		nn -= Nx*Ny*Nz;	// Perioidic BC along the z-boundary
            f13 = phi[nn];					// get neighbor for phi - 13
            //........................................................................
            nn = n-Nx*Ny+1;						// neighbor index (get convention)
            if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
            if (k-1<0)			nn += Nx*Ny*Nz;	// Perioidic BC along the z-boundary
            f14 = phi[nn];					// get neighbor for phi - 14
            //........................................................................
            nn = n-Nx*Ny-Nx;					// neighbor index (get convention)
            if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
            if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
            f15 = phi[nn];					// get neighbor for phi - 15
            //........................................................................
            nn = n+Nx*Ny+Nx;					// neighbor index (get convention)
            if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
            if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
            f16 = phi[nn];					// get neighbor for phi - 16
            //........................................................................
            nn = n+Nx*Ny-Nx;					// neighbor index (get convention)
            if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
            if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
            f17 = phi[nn];					// get neighbor for phi - 17
            //........................................................................
            nn = n-Nx*Ny+Nx;					// neighbor index (get convention)
            if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
            if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
            f18 = phi[nn];					// get neighbor for phi - 18
            //............Compute the Color Gradient...................................
            nx = -(f1-f2+0.5*(f7-f8+f9-f10+f11-f12+f13-f14));
            ny = -(f3-f4+0.5*(f7-f8-f9+f10+f15-f16+f17-f18));
            nz = -(f5-f6+0.5*(f11-f12-f13+f14+f15-f16-f17+f18));
            //...Store the Color Gradient....................
            ColorGrad[n] = nx;
            ColorGrad[N+n] = ny;
            ColorGrad[2*N+n] = nz;
            //...............................................
            //...........Normalize the Color Gradient.................................
            C = sqrt(nx*nx+ny*ny+nz*nz);
            if (C==0.0) C=1.0;
            nx = nx/C;
            ny = ny/C;
            nz = nz/C;
            //......No color gradient at z-boundary if pressure BC are set.............
            //	if (pBC && k==0) nx = ny = nz = 0.f;
            //	if (pBC && k==Nz-1) nx = ny = nz = 0.f;
            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
            f2 = distodd[n];
            f4 = distodd[N+n];
            f6 = distodd[2*N+n];
            f0 = disteven[n];
            f1 = disteven[N+n];
            f3 = disteven[2*N+n];
            f5 = disteven[3*N+n];
            //........................................................................
            //....................compute the moments...............................................
            rho = f0+f2+f1+f4+f3+f6+f5;
            m1 = -30*f0-11*(f2+f1+f4+f3+f6+f5);
            m2 = 12*f0-4*(f2+f1 +f4+f3+f6 +f5);
            jx = f1-f2;
            m4 = 4*(-f1+f2);
            jy = f3-f4;
            m6 = -4*(f3-f4);
            jz = f5-f6;
            m8 = -4*(f5-f6);
            m9 = 2*(f1+f2)-f3-f4-f5-f6;
            m10 = -4*(f1+f2)+2*(f4+f3+f6+f5);
            m11 = f4+f3-f6-f5;
            m12 = -2*(f4+f3-f6-f5);
            //........................................................................
            f8 = distodd[3*N+n];
            f10 = distodd[4*N+n];
            f7 = disteven[4*N+n];
            f9 = disteven[5*N+n];
            //........................................................................
            rho += f8+f7+f10+f9;
            m1 += 8*(f8+f7+f10+f9);
            m2 += f8+f7+f10+f9;
            jx += f7-f8+f9-f10;
            m4 += f7-f8+f9-f10;
            jy += f7-f8-f9+f10;
            m6 += f7-f8-f9+f10;
            m9 += f7+f8+f9+f10;
            m10 += f8+f7+f10+f9;
            m11 += f8+f7+f10+f9;
            m12 += f8+f7+f10+f9;
            m13 = f8+f7-f10-f9;
            m16 = f7-f8+f9-f10;
            m17 = -f7+f8+f9-f10;
            //........................................................................
            f11 = disteven[6*N+n];
            f13 = disteven[7*N+n];
            f12 = distodd[5*N+n];
            f14 = distodd[6*N+n];
            //........................................................................
            //........................................................................
            f15 = disteven[8*N+n];
            f17 = disteven[9*N+n];
            f16 = distodd[7*N+n];
            f18 = distodd[8*N+n];
            //........................................................................
            //....................compute the moments...............................................
            rho += f12+f11+f14+f13+f16+f15+f18+f17;
            m1 += 8*(f12+f11+f14+f13+f16+f15+f18+f17);
            m2 += f12+f11+f14+f13+f16+f15+f18+f17;
            jx += f11-f12+f13-f14;
            m4 += f11-f12+f13-f14;
            jy += f15-f16+f17-f18;
            m6 += f15-f16+f17-f18;
            jz += f11-f12-f13+f14+f15-f16-f17+f18;
            m8 += f11-f12-f13+f14+f15-f16-f17+f18;
            m9 += f11+f12+f13+f14-2*(f15+f16+f17+f18);
            m10 += f12+f11+f14+f13-2*(f16+f15+f18+f17);
            m11 += -f12-f11-f14-f13;
            m12 += -f12-f11-f14-f13;
            m14 = f16+f15-f18-f17;
            m15 = f12+f11-f14-f13;
            m16 += -f11+f12-f13+f14;
            m17 += f15-f16+f17-f18;
            m18 = f11-f12-f13+f14-f15+f16+f17-f18;
            //........................................................................
            
            /*				f2 = distodd[n];
             f4 = distodd[N+n];
             f6 = distodd[2*N+n];
             f8 = distodd[3*N+n];
             //........................................................................
             f0 = disteven[n];
             f1 = disteven[N+n];
             f3 = disteven[2*N+n];
             f5 = disteven[3*N+n];
             f7 = disteven[4*N+n];
             //........................................................................
             //........................................................................
             //....................compute the moments...............................................
             rho = f0+f2+f1+f4+f3+f6+f5+f8+f7;
             m1 = -30*f0-11*(f2+f1+f4+f3+f6+f5)+8*(f8+f7);
             m2 = 12*f0-4*(f2+f1 +f4+f3+f6 +f5)+f8+f7;
             jx = f1-f2+f7-f8;
             m4 = 4*(-f1+f2)+f7-f8;
             jy = f3-f4+f7-f8;
             m6 = -4*(f3-f4)+f7-f8;
             jz = f5-f6;
             m8 = -4*(f5-f6);
             m9 = 2*(f1+f2)-f3-f4-f5-f6+f7+f8;
             m10 = -4*(f1+f2)+2*(f4+f3+f6+f5)+f8+f7;
             m11 = f4+f3-f6-f5+f8+f7;
             m12 = -2*(f4+f3-f6-f5)+f8+f7;
             m13 = f8+f7;
             m16 = f7-f8;
             m17 = -f7+f8;
             //........................................................................
             f9 = disteven[5*N+n];
             f11 = disteven[6*N+n];
             f13 = disteven[7*N+n];
             f15 = disteven[8*N+n];
             f17 = disteven[9*N+n];
             f10 = distodd[4*N+n];
             f12 = distodd[5*N+n];
             f14 = distodd[6*N+n];
             f16 = distodd[7*N+n];
             f18 = distodd[8*N+n];
             //........................................................................
             rho += f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;
             m1 += 8*(f10+f9+f12+f11+f14+f13+f16+f15+f18 +f17);
             m2 += f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;
             jx += f9-f10+f11-f12+f13-f14;
             m4 += f9-f10+f11-f12+f13-f14;
             jy += -f9+f10+f15-f16+f17-f18;
             m6 += -f9+f10+f15-f16+f17-f18;
             jz += f11-f12-f13+f14+f15-f16-f17+f18;
             m8 += f11-f12-f13+f14+f15-f16-f17+f18;
             m9 += f9+f10+f11+f12+f13+f14-2*(f15+f16+f17+f18);
             m10 += f10+f9+f12+f11+f14+f13-2*(f16+f15+f18+f17);
             m11 += f10+f9-f12-f11-f14-f13;
             m12 += f10+f9-f12-f11-f14-f13;
             m13 += -f10-f9;
             m14 = f16+f15-f18-f17;
             m15 = f12+f11-f14-f13;
             m16 += f9-f10-f11+f12-f13+f14;
             m17 += f9-f10+f15-f16+f17-f18;
             m18 = f11-f12-f13+f14-f15+f16+f17-f18;
             */			//........................................................................
            //					PERFORM RELAXATION PROCESS
            //........................................................................
            //..........Toelke, Fruediger et. al. 2006...............
            if (C == 0.0)	nx = ny = nz = 0.0;
            m1 = m1 + rlx_setA*((19*(jx*jx+jy*jy+jz*jz)/rho - 11*rho) -alpha*C - m1);
            m2 = m2 + rlx_setA*((3*rho - 5.5*(jx*jx+jy*jy+jz*jz)/rho)- m2);
            m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
            m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
            m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
            m9 = m9 + rlx_setA*(((2*jx*jx-jy*jy-jz*jz)/rho) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
            m10 = m10 + rlx_setA*( - m10);
            m11 = m11 + rlx_setA*(((jy*jy-jz*jz)/rho) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
            m12 = m12 + rlx_setA*( - m12);
            m13 = m13 + rlx_setA*( (jx*jy/rho) + 0.5*alpha*C*nx*ny - m13);
            m14 = m14 + rlx_setA*( (jy*jz/rho) + 0.5*alpha*C*ny*nz - m14);
            m15 = m15 + rlx_setA*( (jx*jz/rho) + 0.5*alpha*C*nx*nz - m15);
            m16 = m16 + rlx_setB*( - m16);
            m17 = m17 + rlx_setB*( - m17);
            m18 = m18 + rlx_setB*( - m18);
            //.................inverse transformation......................................................
            f0 = 0.05263157894736842*rho-0.012531328320802*m1+0.04761904761904762*m2;
            f1 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
            +0.1*(jx-m4)+0.0555555555555555555555555*(m9-m10);
            f2 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
            +0.1*(m4-jx)+0.0555555555555555555555555*(m9-m10);
            f3 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
            +0.1*(jy-m6)+0.02777777777777778*(m10-m9)+0.08333333333333333333333*(m11-m12);
            f4 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
            +0.1*(m6-jy)+0.02777777777777778*(m10-m9)+0.08333333333333333333333*(m11-m12);
            f5 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
            +0.1*(jz-m8)+0.02777777777777778*(m10-m9)+0.08333333333333333333333*(m12-m11);
            f6 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
            +0.1*(m8-jz)+0.02777777777777778*(m10-m9)+0.08333333333333333333333*(m12-m11);
            f7 = 0.05263157894736842*rho+0.003341687552213868*m1+0.003968253968253968*m2+0.1*(jx+jy)+0.025*(m4+m6)
            +0.02777777777777778*m9+0.01388888888888889*m10+0.08333333333333333333333*m11
            +0.04166666666666666*m12+0.25*m13+0.125*(m16-m17);
            f8 = 0.05263157894736842*rho+0.003341687552213868*m1+0.003968253968253968*m2-0.1*(jx+jy)-0.025*(m4+m6)
            +0.02777777777777778*m9+0.01388888888888889*m10+0.08333333333333333333333*m11
            +0.04166666666666666*m12+0.25*m13+0.125*(m17-m16);
            f9 = 0.05263157894736842*rho+0.003341687552213868*m1+0.003968253968253968*m2+0.1*(jx-jy)+0.025*(m4-m6)
            +0.02777777777777778*m9+0.01388888888888889*m10+0.08333333333333333333333*m11
            +0.04166666666666666*m12-0.25*m13+0.125*(m16+m17);
            f10 = 0.05263157894736842*rho+0.003341687552213868*m1+0.003968253968253968*m2+0.1*(jy-jx)+0.025*(m6-m4)
            +0.02777777777777778*m9+0.01388888888888889*m10+0.08333333333333333333333*m11
            +0.04166666666666666*m12-0.25*m13-0.125*(m16+m17);
            f11 = 0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2+0.1*(jx+jz)+0.025*(m4+m8)
            +0.02777777777777778*m9+0.01388888888888889*m10-0.08333333333333333333333*m11
            -0.04166666666666666*m12+0.25*m15+0.125*(m18-m16);
            f12 = 0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2-0.1*(jx+jz)-0.025*(m4+m8)
            +0.02777777777777778*m9+0.01388888888888889*m10-0.08333333333333333333333*m11
            -0.04166666666666666*m12+0.25*m15+0.125*(m16-m18);
            f13 = 0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2+0.1*(jx-jz)+0.025*(m4-m8)
            +0.02777777777777778*m9+0.01388888888888889*m10-0.08333333333333333333333*m11
            -0.04166666666666666*m12-0.25*m15-0.125*(m16+m18);
            f14 = 0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2+0.1*(jz-jx)+0.025*(m8-m4)
            +0.02777777777777778*m9+0.01388888888888889*m10-0.08333333333333333333333*m11
            -0.04166666666666666*m12-0.25*m15+0.125*(m16+m18);
            f15 = 0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2+0.1*(jy+jz)+0.025*(m6+m8)
            -0.0555555555555555555555555*m9-0.02777777777777778*m10+0.25*m14+0.125*(m17-m18);
            f16 =  0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2-0.1*(jy+jz)-0.025*(m6+m8)
            -0.0555555555555555555555555*m9-0.02777777777777778*m10+0.25*m14+0.125*(m18-m17);
            f17 = 0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2+0.1*(jy-jz)+0.025*(m6-m8)
            -0.0555555555555555555555555*m9-0.02777777777777778*m10-0.25*m14+0.125*(m17+m18);
            f18 = 0.05263157894736842*rho+0.003341687552213868*m1
            +0.003968253968253968*m2+0.1*(jz-jy)+0.025*(m8-m6)
            -0.0555555555555555555555555*m9-0.02777777777777778*m10-0.25*m14-0.125*(m17+m18);
            //.......................................................................................................
            // incorporate external force
            f1 += 0.16666666666666667*Fx;
            f2 -= 0.16666666666666667*Fx;
            f3 += 0.16666666666666667*Fy;
            f4 -= 0.16666666666666667*Fy;
            f5 += 0.16666666666666667*Fz;
            f6 -= 0.16666666666666667*Fz;
            f7 += 0.08333333333333333*(Fx+Fy);
            f8 -= 0.08333333333333333*(Fx+Fy);
            f9 += 0.08333333333333333*(Fx-Fy);
            f10 -= 0.08333333333333333*(Fx-Fy);
            f11 += 0.08333333333333333*(Fx+Fz);
            f12 -= 0.08333333333333333*(Fx+Fz);
            f13 += 0.08333333333333333*(Fx-Fz);
            f14 -= 0.08333333333333333*(Fx-Fz);
            f15 += 0.08333333333333333*(Fy+Fz);
            f16 -= 0.08333333333333333*(Fy+Fz);
            f17 += 0.08333333333333333*(Fy-Fz);
            f18 -= 0.08333333333333333*(Fy-Fz);
            //*********** WRITE UPDATED VALUES TO MEMORY ******************
            // Write the updated distributions
            //....EVEN.....................................
            disteven[n] = f0;
            disteven[N+n] = f2;
            disteven[2*N+n] = f4;
            disteven[3*N+n] = f6;
            disteven[4*N+n] = f8;
            disteven[5*N+n] = f10;
            disteven[6*N+n] = f12;
            disteven[7*N+n] = f14;
            disteven[8*N+n] = f16;
            disteven[9*N+n] = f18;
            //....ODD......................................
            distodd[n] = f1;
            distodd[N+n] = f3;
            distodd[2*N+n] = f5;
            distodd[3*N+n] = f7;
            distodd[4*N+n] = f9;
            distodd[5*N+n] = f11;
            distodd[6*N+n] = f13;
            distodd[7*N+n] = f15;
            distodd[8*N+n] = f17;
            //...Store the Velocity..........................
            Velocity[n] = jx;
            Velocity[N+n] = jy;
            Velocity[2*N+n] = jz;
            //***************************************************************
        }	// check if n is in the solid
    } // loop over n
}

extern "C" void ScaLBL_D3Q7_ColorCollideMass(char *ID, double *A_even, double *A_odd, double *B_even, double *B_odd, 
                                             double *Den, double *Phi, double *ColorGrad, double *Velocity, double beta, int N, bool pBC)
{
    char id;
    
    int idx,n,q,Cqx,Cqy,Cqz;
    //	int sendLoc;
    
    double f0,f1,f2,f3,f4,f5,f6;
    double na,nb,nab;		// density values
    double ux,uy,uz;	// flow velocity
    double nx,ny,nz,C;	// color gradient components
    double a1,a2,b1,b2;
    double sp,delta;
    //double feq[6];		// equilibrium distributions
    // Set of Discrete velocities for the D3Q19 Model
    //int D3Q7[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    
    for (n=0; n<N; n++){
        id = ID[n];
        if (id != 0 ){
            
            //.....Load the Color gradient.........
            nx = ColorGrad[n];
            ny = ColorGrad[N+n];
            nz = ColorGrad[2*N+n];
            C = sqrt(nx*nx+ny*ny+nz*nz);
            if (C==0.0) C=1.0;
            nx = nx/C;
            ny = ny/C;
            nz = nz/C;
            //....Load the flow velocity...........
            ux = Velocity[n];
            uy = Velocity[N+n];
            uz = Velocity[2*N+n];
            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
            f2 = A_odd[n];
            f4 = A_odd[N+n];
            f6 = A_odd[2*N+n];
            f0 = A_even[n];
            f1 = A_even[N+n];
            f3 = A_even[2*N+n];
            f5 = A_even[3*N+n];
            na = f0+f1+f2+f3+f4+f5+f6;
            //........................................................................
            f2 = B_odd[n];
            f4 = B_odd[N+n];
            f6 = B_odd[2*N+n];
            f0 = B_even[n];
            f1 = B_even[N+n];
            f3 = B_even[2*N+n];
            f5 = B_even[3*N+n];
            nb = f0+f1+f2+f3+f4+f5+f6;
            nab = 1.0/(na+nb);
            //........................................................................
            //....Instantiate the density distributions
            // Generate Equilibrium Distributions and stream
            // Stationary value - distribution 0
            A_even[n] = 0.3333333333333333*na;
            B_even[n] = 0.3333333333333333*nb;
            // Non-Stationary equilibrium distributions
            //feq[0] = 0.1111111111111111*(1+4.5*ux);
            //feq[1] = 0.1111111111111111*(1-4.5*ux);
            //feq[2] = 0.1111111111111111*(1+4.5*uy);
            //feq[3] = 0.1111111111111111*(1-4.5*uy);
            //feq[4] = 0.1111111111111111*(1+4.5*uz);
            //feq[5] = 0.1111111111111111*(1-4.5*uz);
            
            //...............................................
            // q = 0,2,4
            // Cq = {1,0,0}, {0,1,0}, {0,0,1}
            delta = beta*na*nb*nab*0.1111111111111111*nx;
            if (!(na*nb*nab>0)) delta=0;
            a1 = na*(0.1111111111111111*(1+4.5*ux))+delta;
            b1 = nb*(0.1111111111111111*(1+4.5*ux))-delta;
            a2 = na*(0.1111111111111111*(1-4.5*ux))-delta;
            b2 = nb*(0.1111111111111111*(1-4.5*ux))+delta;
            
            A_odd[n] 	= a1;
            A_even[N+n] = a2;
            B_odd[n] 	= b1;
            B_even[N+n] = b2;
            //...............................................
            // q = 2
            // Cq = {0,1,0}
            delta = beta*na*nb*nab*0.1111111111111111*ny;
            if (!(na*nb*nab>0)) delta=0;
            a1 = na*(0.1111111111111111*(1+4.5*uy))+delta;
            b1 = nb*(0.1111111111111111*(1+4.5*uy))-delta;
            a2 = na*(0.1111111111111111*(1-4.5*uy))-delta;
            b2 = nb*(0.1111111111111111*(1-4.5*uy))+delta;
            
            A_odd[N+n] 	= a1;
            A_even[2*N+n] = a2;
            B_odd[N+n] 	= b1;
            B_even[2*N+n] = b2;
            //...............................................
            // q = 4
            // Cq = {0,0,1}
            delta = beta*na*nb*nab*0.1111111111111111*nz;
            if (!(na*nb*nab>0)) delta=0;
            a1 = na*(0.1111111111111111*(1+4.5*uz))+delta;
            b1 = nb*(0.1111111111111111*(1+4.5*uz))-delta;
            a2 = na*(0.1111111111111111*(1-4.5*uz))-delta;
            b2 = nb*(0.1111111111111111*(1-4.5*uz))+delta;
            
            A_odd[2*N+n] = a1;
            A_even[3*N+n] = a2;
            B_odd[2*N+n] = b1;
            B_even[3*N+n] = b2;
            //...............................................
            
            /*		// Construction and streaming for the components
             for (idx=0; idx<3; idx++){
             //...............................................
             // Distribution index
             q = 2*idx;
             // Associated discrete velocity
             Cqx = D3Q7[idx][0];
             Cqy = D3Q7[idx][1];
             Cqz = D3Q7[idx][2];
             // Generate the Equilibrium Distribution
             a1 = nBC*eq[q];
             b1 = nb*feq[q];
             a2 = nBC*eq[q+1];
             b2 = nb*feq[q+1];
             // Recolor the distributions
             if (C > 0.0){
             sp = nx*double(Cqx)+ny*double(Cqy)+nz*double(Cqz);
             //if (idx > 2)	sp = 0.7071067811865475*sp;
             //delta = sp*min( min(a1,a2), min(b1,b2) );
             delta = na*nb/(na+nb)*0.1111111111111111*sp;
             //if (a1>0 && b1>0){
             a1 += beta*delta;
             a2 -= beta*delta;
             b1 -= beta*delta;
             b2 += beta*delta;
             }
             // Save the re-colored distributions
             A_odd[N*idx+n] 		= a1;
             A_even[N*(idx+1)+n] = a2;
             B_odd[N*idx+n] 		= b1;
             B_even[N*(idx+1)+n] = b2;
             //...............................................
             }
             */
        }
    }
}

//*************************************************************************
extern "C" void DensityStreamD3Q7(char *ID, double *Den, double *Copy, double *Phi, double *ColorGrad, double *Velocity,
                                  double beta, int Nx, int Ny, int Nz, bool pBC, int S)
{
    char id;
    
    int idx;
    int in,jn,kn,n,nn,N;
    int q,Cqx,Cqy,Cqz;
    //	int sendLoc;
    
    double na,nb;		// density values
    double ux,uy,uz;	// flow velocity
    double nx,ny,nz,C;	// color gradient components
    double a1,a2,b1,b2;
    double sp,delta;
    double feq[6];		// equilibrium distributions
    // Set of Discrete velocities for the D3Q19 Model
    int D3Q7[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    N = Nx*Ny*Nz;
    
    for (n=0; n<N; n++){
        id = ID[n];
        // Local Density Values
        na = Copy[2*n];
        nb = Copy[2*n+1];
        if (id > 0 && na+nb > 0.0){
            //.......Back out the 3-D indices for node n..............
            int	k = n/(Nx*Ny);
            int j = (n-Nx*Ny*k)/Nx;
            int i = n-Nx*Ny*k-Nx*j;
            //.....Load the Color gradient.........
            nx = ColorGrad[n];
            ny = ColorGrad[N+n];
            nz = ColorGrad[2*N+n];
            C = sqrt(nx*nx+ny*ny+nz*nz);
            nx = nx/C;
            ny = ny/C;
            nz = nz/C;
            //....Load the flow velocity...........
            ux = Velocity[n];
            uy = Velocity[N+n];
            uz = Velocity[2*N+n];
            //....Instantiate the density distributions
            // Generate Equilibrium Distributions and stream
            // Stationary value - distribution 0
            //			Den[2*n] += 0.3333333333333333*na;
            //			Den[2*n+1] += 0.3333333333333333*nb;
            Den[2*n] += 0.3333333333333333*na;
            Den[2*n+1] += 0.3333333333333333*nb;
            // Non-Stationary equilibrium distributions
            feq[0] = 0.1111111111111111*(1+3*ux);
            feq[1] = 0.1111111111111111*(1-3*ux);
            feq[2] = 0.1111111111111111*(1+3*uy);
            feq[3] = 0.1111111111111111*(1-3*uy);
            feq[4] = 0.1111111111111111*(1+3*uz);
            feq[5] = 0.1111111111111111*(1-3*uz);
            // Construction and streaming for the components
            for (idx=0; idx<3; idx++){
                // Distribution index
                q = 2*idx;
                // Associated discrete velocity
                Cqx = D3Q7[idx][0];
                Cqy = D3Q7[idx][1];
                Cqz = D3Q7[idx][2];
                // Generate the Equilibrium Distribution
                a1 = na*feq[q];
                b1 = nb*feq[q];
                a2 = na*feq[q+1];
                b2 = nb*feq[q+1];
                // Recolor the distributions
                if (C > 0.0){
                    sp = nx*double(Cqx)+ny*double(Cqy)+nz*double(Cqz);
                    //if (idx > 2)	sp = 0.7071067811865475*sp;
                    //delta = sp*min( min(a1,a2), min(b1,b2) );
                    delta = na*nb/(na+nb)*0.1111111111111111*sp;
                    //if (a1>0 && b1>0){
                    a1 += beta*delta;
                    a2 -= beta*delta;
                    b1 -= beta*delta;
                    b2 += beta*delta;
                }
                
                // .......Get the neighbor node..............
                //nn = n + Stride[idx];
                in = i+Cqx;
                jn = j+Cqy;
                kn = k+Cqz;
                
                // Adjust for periodic BC, if necessary
                //				if (in<0) in+= Nx;
                //				if (jn<0) jn+= Ny;
                //				if (kn<0) kn+= Nz;
                //				if (!(in<Nx)) in-= Nx;
                //				if (!(jn<Ny)) jn-= Ny;
                //				if (!(kn<Nz)) kn-= Nz;
                // Perform streaming or bounce-back as needed
                id = ID[kn*Nx*Ny+jn*Nx+in];
                if (id == 0){							//.....Bounce-back Rule...........
                    //						Den[2*n] += a1;
                    //						Den[2*n+1] += b1;
                    Den[2*n] += a1;
                    Den[2*n+1] += b1;
                }
                else{
                    //......Push the "distribution" to neighboring node...........
                    // Index of the neighbor in the local process
                    //nn = (kn-zmin[rank]+1)*Nxp*Nyp + (jn-ymin[rank]+1)*Nxp + (in-xmin[rank]+1);
                    nn = kn*Nx*Ny+jn*Nx+in;
                    // Push to neighboring node
                    //						Den[2*nn] += a1;
                    //						Den[2*nn+1] += b1;
                    Den[2*nn] += a1;
                    Den[2*nn+1] += b1;
                }
                
                // .......Get the neighbor node..............
                q = 2*idx+1;
                in = i-Cqx;
                jn = j-Cqy;
                kn = k-Cqz;
                // Adjust for periodic BC, if necessary
                //				if (in<0) in+= Nx;
                //				if (jn<0) jn+= Ny;
                //				if (kn<0) kn+= Nz;
                //				if (!(in<Nx)) in-= Nx;
                //				if (!(jn<Ny)) jn-= Ny;
                //				if (!(kn<Nz)) kn-= Nz;
                // Perform streaming or bounce-back as needed
                id = ID[kn*Nx*Ny+jn*Nx+in];
                if (id == 0){
                    //.....Bounce-back Rule...........
                    //						Den[2*n] += a2;
                    //					Den[2*n+1] += b2;
                    Den[2*n] += a2;
                    Den[2*n+1] += b2;
                }
                else{
                    //......Push the "distribution" to neighboring node...........
                    // Index of the neighbor in the local process
                    //nn = (kn-zmin[rank]+1)*Nxp*Nyp + (jn-ymin[rank]+1)*Nxp + (in-xmin[rank]+1);
                    nn = kn*Nx*Ny+jn*Nx+in;
                    // Push to neighboring node
                    //					Den[2*nn] += a2;
                    //					Den[2*nn+1] += b2;
                    Den[2*nn] += a2;
                    Den[2*nn+1] += b2;
                }
            }
        }
    }
}

extern "C" void ScaLBL_ComputePhaseField(char *ID, double *Phi, double *Den, int N)
{
    int n;
    double Na,Nb;
    //...................................................................
    // Update Phi
    for (n=0; n<N; n++){
        
        if (ID[n] > 0 ){
            // Get the density value (Streaming already performed)
            Na = Den[n];
            Nb = Den[N+n];
            Phi[n] = (Na-Nb)/(Na+Nb);
        }
    }
    //...................................................................
}

extern "C" void ScaLBL_SetSlice_z(double *Phi, double value, int Nx, int Ny, int Nz, int Slice){
    int n;
    for (n=Slice*Nx*Ny; n<(Slice+1)*Nx*Ny; n++){
        Phi[n] = value;
    }
}







extern "C" void ScaLBL_D3Q19_AAeven_Color(int *Map, double *dist, double *Aq, double *Bq, double *Den, double *Phi,
                                          double *Vel, double* Press, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
                                          double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np) {
    
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
    
    double jxp,jyp,jzp;
    
    const double mrt_V1=0.05263157894736842;
    const double mrt_V2=0.012531328320802;
    const double mrt_V3=0.04761904761904762;
    const double mrt_V4=0.004594820384294068;
    const double mrt_V5=0.01587301587301587;
    const double mrt_V6=0.0555555555555555556;
    const double mrt_V7=0.02777777777777778;
    const double mrt_V8=0.08333333333333333333333;
    const double mrt_V9=0.003341687552213868;
    const double mrt_V10=0.003968253968253968;
    const double mrt_V11=0.01388888888888889;
    const double mrt_V12=0.04166666666666666;
    
    
    for (int n=start; n<finish; n++){
        
        nA = Den[n];
        nB = Den[Np + n];
        
        // compute phase indicator field
        phi=(nA-nB)/(nA+nB);
        
        // local density
        rho0=1.0 + 0.5*(1.0-phi)*(rhoB-rhoA);  // 1.0 for single phase
        // local relaxation time
        tau=tauA + 0.5*(1.0-phi)*(tauB-tauA); // tau = tauA for single phase
        rlx_setA = 1.0/tau;
        rlx_setB = 8.*(2.-rlx_setA)/(8.-rlx_setA); // 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
        
        // Get the 1D index based on regular data layout
        ijk = Map[n];
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
       
        /* Removed normalization */
        nx = -(m1-m2+(m7-m8+m9-m10+m11-m12+m13-m14));
        ny = -(m3-m4+(m7-m8-m9+m10+m15-m16+m17-m18));
        nz = -(m5-m6+(m11-m12-m13+m14+m15-m16-m17+m18));
        
        //...........Normalize the Color Gradient.................................
        C = sqrt(nx*nx+ny*ny+nz*nz);
//        if (C==0.0) C=1.0;
//        nx = nx/C;
//        ny = ny/C;
//        nz = nz/C;
        
        // q = 0
        fq = dist[n];
        rho = fq;
        m1  = -30.0*fq;
        m2  = 12.0*fq;
        
        // q = 1
        fq = dist[2*Np+n];
        rho += fq;
        m1 -= 11.0*fq;
        m2 -= 4.0*fq;
        jx = fq;
        m4 = -4.0*fq;
        m9 = 2.0*fq;
        m10 = -4.0*fq;
        
        // q = 2
        fq = dist[1*Np+n];
        rho += fq;
        m1 -= 11.0*(fq);
        m2 -= 4.0*(fq);
        jx -= fq;
        m4 += 4.0*(fq);
        m9 += 2.0*(fq);
        m10 -= 4.0*(fq);
        
        // q = 3
        fq = dist[4*Np+n];
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
        fq = dist[3*Np+n];
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
        fq = dist[6*Np+n];
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
        fq = dist[5*Np+n];
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
        fq = dist[8*Np+n];
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
        fq = dist[7*Np+n];
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
        fq = dist[10*Np+n];
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
        fq = dist[9*Np+n];
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
        fq = dist[12*Np+n];
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
        fq = dist[11*Np+n];
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
        fq = dist[14*Np+n];
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
        fq = dist[13*Np+n];
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
        fq = dist[16*Np+n];
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
        fq = dist[15*Np+n];
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
        fq = dist[18*Np+n];
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
        fq = dist[17*Np+n];
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
        
        
        // write the MOMENTUM
        Press[n] = 0.3333333333333333*rho;
        Vel[n] = jx + 0.5*Fx;
        Vel[Np+n] = jy + 0.5*Fy;
        Vel[2*Np+n] = jz + 0.5*Fz;
        
        
        jxp = jx + 0.5*Fx;
        jyp = jy + 0.5*Fy;
        jzp = jz + 0.5*Fz;
        
        
        // ENTER MOMENTUM INTO THE MOMENTS, NOT VELOCITY
        //..........Toelke, Fruediger et. al. 2006................................
        if (C == 0.0)    nx = ny = nz = 0.0;
        // C = nx = ny = nz = 0;
        m1 = m1 + rlx_setA*((19*(jxp*jxp+jyp*jyp+jzp*jzp)/rho0 - 11*rho) -alpha*C - m1);
        m2 = m2 + rlx_setA*((3*rho - 5.5*(jxp*jxp+jyp*jyp+jzp*jzp)/rho0)- m2);
        m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
        m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
        m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
        // pxx
        m9 = m9 + rlx_setA*(   ((2*jxp*jxp-jyp*jyp-jzp*jzp)/(3.0*rho0)) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
        // PIxx
        m10 = m10 + rlx_setA*( -0.5*((2*jxp*jxp-jyp*jyp-jzp*jzp)/(3.0*rho0))  - m10);
        // pww
        m11 = m11 + rlx_setA*(   ((jyp*jyp-jzp*jzp)/rho0) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
        // PIww
        m12 = m12 + rlx_setA*( -0.5*((jyp*jyp-jzp*jzp)/rho0)  - m12);
        m13 = m13 + rlx_setA*( (jxp*jyp/rho0) + 0.5*alpha*C*nx*ny - m13);
        m14 = m14 + rlx_setA*( (jyp*jzp/rho0) + 0.5*alpha*C*ny*nz - m14);
        m15 = m15 + rlx_setA*( (jxp*jzp/rho0) + 0.5*alpha*C*nx*nz - m15);
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
        
        
        
        
        
        //.......................................................................................................
        //.................inverse transformation......................................................
        
        fq = 0. - 0.012531328320802004*m1 + 0.047619047619047616*m2 + 0.05263157894736842*rho;
        dist[n] = fq;
        
        fq = 0. + 0.1*jx - 0.004594820384294068*m1 - 0.05555555555555555*m10 - 0.015873015873015872*m2 -
        0.1*m4 + 0.05555555555555555*m9 + 0.05263157894736842*rho + 0.16666666666666667*Fx;
        dist[n+Np] = fq;
        
        fq =  0. - 0.1*jx - 0.004594820384294068*m1 - 0.05555555555555555*m10 - 0.015873015873015872*m2 +
        0.1*m4 + 0.05555555555555555*m9 + 0.05263157894736842*rho - 0.16666666666666667*Fx;
        dist[n+2*Np] = fq;
        
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
        
        fq =  0. - 0.1*jy - 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 + 0.25*m14 -
        0.125*m17 + 0.125*m18 + 0.003968253968253968*m2 - 0.025*m6 - 0.025*m8 -
        0.05555555555555555*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fy+Fz);
        dist[n+16*Np] = fq;
        
        fq =   0. + 0.1*jy - 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 - 0.25*m14 +
        0.125*m17 + 0.125*m18 + 0.003968253968253968*m2 + 0.025*m6 - 0.025*m8 -
        0.05555555555555555*m9 + 0.05263157894736842*rho + 0.08333333333333333*(Fy-Fz);
        dist[n+17*Np] = fq;
        
        fq =   0. - 0.1*jy + 0.1*jz + 0.003341687552213868*m1 - 0.027777777777777776*m10 - 0.25*m14 -
        0.125*m17 - 0.125*m18 + 0.003968253968253968*m2 - 0.025*m6 + 0.025*m8 -
        0.05555555555555555*m9 + 0.05263157894736842*rho - 0.08333333333333333*(Fy-Fz);
        dist[n+18*Np] = fq;
        
        
        ux = jx/rho0;
        ux = jy/rho0;
        ux = jz/rho0;
        
       
        
         Aq[n] = nA/3.;
         Aq[n+1*Np] = ((beta*nA*nB*nx)/(nA + nB) + nA*(1 + 3*ux))/18.;
         Aq[n+2*Np] = (-((beta*nA*nB*nx)/(nA + nB)) + nA*(1 - 3*ux))/18.;
         Aq[n+3*Np] = ((beta*nA*nB*ny)/(nA + nB) + nA*(1 + 3*uy))/18.;
         Aq[n+4*Np] = (-((beta*nA*nB*ny)/(nA + nB)) + nA*(1 - 3*uy))/18.;
         Aq[n+5*Np] = ((beta*nA*nB*nz)/(nA + nB) + nA*(1 + 3*uz))/18.;
         Aq[n+6*Np] = (-((beta*nA*nB*nz)/(nA + nB)) + nA*(1 - 3*uz))/18.;
         Aq[n+7*Np] = ((beta*nA*nB*(nx + ny))/(nA + nB) + nA*(1 + 3*(ux + uy)))/36.;
         Aq[n+8*Np] = ((beta*nA*nB*(-nx - ny))/(nA + nB) + nA*(1 + 3*(-ux - uy)))/36.;
         Aq[n+9*Np] = ((beta*nA*nB*(nx - ny))/(nA + nB) + nA*(1 + 3*(ux - uy)))/36.;
         Aq[n+10*Np] = ((beta*nA*nB*(-nx + ny))/(nA + nB) + nA*(1 + 3*(-ux + uy)))/36.;
         Aq[n+11*Np] = ((beta*nA*nB*(nx + nz))/(nA + nB) + nA*(1 + 3*(ux + uz)))/36.;
         Aq[n+12*Np] = ((beta*nA*nB*(-nx - nz))/(nA + nB) + nA*(1 + 3*(-ux - uz)))/36.;
         Aq[n+13*Np] = ((beta*nA*nB*(nx - nz))/(nA + nB) + nA*(1 + 3*(ux - uz)))/36.;
         Aq[n+14*Np] = ((beta*nA*nB*(-nx + nz))/(nA + nB) + nA*(1 + 3*(-ux + uz)))/36.;
         Aq[n+15*Np] = ((beta*nA*nB*(ny + nz))/(nA + nB) + nA*(1 + 3*(uy + uz)))/36.;
         Aq[n+16*Np] = ((beta*nA*nB*(-ny - nz))/(nA + nB) + nA*(1 + 3*(-uy - uz)))/36.;
         Aq[n+17*Np] = ((beta*nA*nB*(ny - nz))/(nA + nB) + nA*(1 + 3*(uy - uz)))/36.;
         Aq[n+18*Np] = ((beta*nA*nB*(-ny + nz))/(nA + nB) + nA*(1 + 3*(-uy + uz)))/36.;

         
        
        Bq[n] = nB/3.;
        Bq[n+1*Np] = (-((beta*nA*nB*nx)/(nA + nB)) + nB*(1 + 3*ux))/18.;
        Bq[n+2*Np] = ((beta*nA*nB*nx)/(nA + nB) + nB*(1 - 3*ux))/18.;
        Bq[n+3*Np] = (-((beta*nA*nB*ny)/(nA + nB)) + nB*(1 + 3*uy))/18.;
        Bq[n+4*Np] = ((beta*nA*nB*ny)/(nA + nB) + nB*(1 - 3*uy))/18.;
        Bq[n+5*Np] = (-((beta*nA*nB*nz)/(nA + nB)) + nB*(1 + 3*uz))/18.;
        Bq[n+6*Np] = ((beta*nA*nB*nz)/(nA + nB) + nB*(1 - 3*uz))/18.;
        Bq[n+7*Np] = (-((beta*nA*nB*(nx + ny))/(nA + nB)) + nB*(1 + 3*(ux + uy)))/36.;
        Bq[n+8*Np] = (-((beta*nA*nB*(-nx - ny))/(nA + nB)) + nB*(1 + 3*(-ux - uy)))/36.;
        Bq[n+9*Np] = (-((beta*nA*nB*(nx - ny))/(nA + nB)) + nB*(1 + 3*(ux - uy)))/36.;
        Bq[n+10*Np] = (-((beta*nA*nB*(-nx + ny))/(nA + nB)) + nB*(1 + 3*(-ux + uy)))/36.;
        Bq[n+11*Np] = (-((beta*nA*nB*(nx + nz))/(nA + nB)) + nB*(1 + 3*(ux + uz)))/36.;
        Bq[n+12*Np] = (-((beta*nA*nB*(-nx - nz))/(nA + nB)) + nB*(1 + 3*(-ux - uz)))/36.;
        Bq[n+13*Np] = (-((beta*nA*nB*(nx - nz))/(nA + nB)) + nB*(1 + 3*(ux - uz)))/36.;
        Bq[n+14*Np] = (-((beta*nA*nB*(-nx + nz))/(nA + nB)) + nB*(1 + 3*(-ux + uz)))/36.;
        Bq[n+15*Np] = (-((beta*nA*nB*(ny + nz))/(nA + nB)) + nB*(1 + 3*(uy + uz)))/36.;
        Bq[n+16*Np] = (-((beta*nA*nB*(-ny - nz))/(nA + nB)) + nB*(1 + 3*(-uy - uz)))/36.;
        Bq[n+17*Np] = (-((beta*nA*nB*(ny - nz))/(nA + nB)) + nB*(1 + 3*(uy - uz)))/36.;
        Bq[n+18*Np] = (-((beta*nA*nB*(-ny + nz))/(nA + nB)) + nB*(1 + 3*(-uy + uz)))/36.;
        
    }
    
}









extern "C" void  ScaLBL_D3Q19_AAodd_Color(int *neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den,
                                         double *Phi, double *Vel, double* Press, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
                                         double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np){
    
    int n,nn,ijk,nread;
    int nr1,nr2,nr3,nr4,nr5,nr6;
    int nr7,nr8,nr9,nr10;
    int nr11,nr12,nr13,nr14;
    int nr15,nr16,nr17,nr18;
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
    
    double jxp,jyp,jzp;
    
    const double mrt_V1=0.05263157894736842;
    const double mrt_V2=0.012531328320802;
    const double mrt_V3=0.04761904761904762;
    const double mrt_V4=0.004594820384294068;
    const double mrt_V5=0.01587301587301587;
    const double mrt_V6=0.055555555555555556;
    const double mrt_V7=0.02777777777777778;
    const double mrt_V8=0.08333333333333333333333;
    const double mrt_V9=0.003341687552213868;
    const double mrt_V10=0.003968253968253968;
    const double mrt_V11=0.01388888888888889;
    const double mrt_V12=0.04166666666666666;
    
    

    
    for (int n=start; n<finish; n++){
        
        nA = Den[n];
        nB = Den[Np + n];
        
        // compute phase indicator field
        phi=(nA-nB)/(nA+nB);
        
        // local density
        rho0=1.0 + 0.5*(1.0-phi)*(rhoB-rhoA);  // 1.0 for single phase
        // local relaxation time
        tau=tauA + 0.5*(1.0-phi)*(tauB-tauA); // tau = tauA for single phase
        
        
        
        
        rlx_setA = 1.0/tau;
        rlx_setB = 8.*(2.-rlx_setA)/(8.-rlx_setA);
        
        
        
        
        // Get the 1D index based on regular data layout
        ijk = Map[n];
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
        
        //...........Normalize the Color Gradient.................................
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
        nr1 = neighborList[n];
        fq = dist[nr1];
        rho += fq;
        m1 -= 11.0*fq;
        m2 -= 4.0*fq;
        jx = fq;
        m4 = -4.0*fq;
        m9 = 2.0*fq;
        m10 = -4.0*fq;
        
        // q = 2
        nr2 = neighborList[n+Np];
        fq = dist[nr2];
        rho += fq;
        m1 -= 11.0*(fq);
        m2 -= 4.0*(fq);
        jx -= fq;
        m4 += 4.0*(fq);
        m9 += 2.0*(fq);
        m10 -= 4.0*(fq);
        
        // q = 3
        nr3 = neighborList[n+2*Np];
        fq = dist[nr3];
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
        nr4 = neighborList[n+3*Np];
        fq = dist[nr4];
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
        nr5 = neighborList[n+4*Np];
        fq = dist[nr5];
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
        nr6 = neighborList[n+5*Np];
        fq = dist[nr6];
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
        nr7 = neighborList[n+6*Np];
        fq = dist[nr7];
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
        nr8 = neighborList[n+7*Np];
        fq = dist[nr8];
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
        nr9 = neighborList[n+8*Np];
        fq = dist[nr9];
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
        nr10 = neighborList[n+9*Np];
        fq = dist[nr10];
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
        nr11 = neighborList[n+10*Np];
        fq = dist[nr11];
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
        nr12 = neighborList[n+11*Np];
        fq = dist[nr12];
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
        nr13 = neighborList[n+12*Np];
        fq = dist[nr13];
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
        nr14 = neighborList[n+13*Np];
        fq = dist[nr14];
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
        nr15 = neighborList[n+14*Np];
        fq = dist[nr15];
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
        nr16 = neighborList[n+15*Np];
        fq = dist[nr16];
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
        nr17 = neighborList[n+16*Np];
        fq = dist[nr17];
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
        nr18 = neighborList[n+17*Np];
        fq = dist[nr18];
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
        
        Press[n] = 0.3333333333333333*rho;
        Vel[n] = jx + 0.5*Fx;
        Vel[Np+n] = jy + 0.5*Fy;
        Vel[2*Np+n] = jz + 0.5*Fz;
        
        
        jxp = jx + 0.5*Fx;
        jyp = jy + 0.5*Fy;
        jzp = jz + 0.5*Fz;
        
        
        // ENTER MOMENTUM INTO THE MOMENTS, NOT VELOCITY
        //..........Toelke, Fruediger et. al. 2006................................
        if (C == 0.0)    nx = ny = nz = 0.0;
//        // C = nx = ny = nz = 0;
        m1 = m1 + rlx_setA*((19*(jxp*jxp+jyp*jyp+jzp*jzp)/rho0 - 11*rho) -alpha*C - m1);
        m2 = m2 + rlx_setA*((3*rho - 5.5*(jxp*jxp+jyp*jyp+jzp*jzp)/rho0)- m2);
        m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
        m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
        m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
        // pxx
        m9 = m9 + rlx_setA*(   ((2*jxp*jxp-jyp*jyp-jzp*jzp)/(3.0*rho0)) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
        // PIxx
        m10 = m10 + rlx_setA*( -0.5*((2*jxp*jxp-jyp*jyp-jzp*jzp)/(3.0*rho0))  - m10);
        // pww
        m11 = m11 + rlx_setA*(   ((jyp*jyp-jzp*jzp)/rho0) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
        // PIww
        m12 = m12 + rlx_setA*( -0.5*((jyp*jyp-jzp*jzp)/rho0)  - m12);
        m13 = m13 + rlx_setA*( (jxp*jyp/rho0) + 0.5*alpha*C*nx*ny - m13);
        m14 = m14 + rlx_setA*( (jyp*jzp/rho0) + 0.5*alpha*C*ny*nz - m14);
        m15 = m15 + rlx_setA*( (jxp*jzp/rho0) + 0.5*alpha*C*nx*nz - m15);
        m16 = m16 + rlx_setB*( - m16);
        m17 = m17 + rlx_setB*( - m17);
        m18 = m18 + rlx_setB*( - m18);
//
        /* Stokes flow */
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
        
        
        ux = jx / rho0;
        uy = jy / rho0;
        uz = jz / rho0;
        
        
         Aq[n] = nA/3.;
         Aq[nr2] = ((beta*nA*nB*nx)/(nA + nB) + nA*(1 + 3*ux))/18.;
         Aq[nr1] = (-((beta*nA*nB*nx)/(nA + nB)) + nA*(1 - 3*ux))/18.;
         Aq[nr4] = ((beta*nA*nB*ny)/(nA + nB) + nA*(1 + 3*uy))/18.;
         Aq[nr3] = (-((beta*nA*nB*ny)/(nA + nB)) + nA*(1 - 3*uy))/18.;
         Aq[nr6] = ((beta*nA*nB*nz)/(nA + nB) + nA*(1 + 3*uz))/18.;
         Aq[nr5] = (-((beta*nA*nB*nz)/(nA + nB)) + nA*(1 - 3*uz))/18.;
         Aq[nr8] = ((beta*nA*nB*(nx + ny))/(nA + nB) + nA*(1 + 3*(ux + uy)))/36.;
         Aq[nr7] = ((beta*nA*nB*(-nx - ny))/(nA + nB) + nA*(1 + 3*(-ux - uy)))/36.;
         Aq[nr10] = ((beta*nA*nB*(nx - ny))/(nA + nB) + nA*(1 + 3*(ux - uy)))/36.;
         Aq[nr9] = ((beta*nA*nB*(-nx + ny))/(nA + nB) + nA*(1 + 3*(-ux + uy)))/36.;
         Aq[nr12] = ((beta*nA*nB*(nx + nz))/(nA + nB) + nA*(1 + 3*(ux + uz)))/36.;
         Aq[nr11] = ((beta*nA*nB*(-nx - nz))/(nA + nB) + nA*(1 + 3*(-ux - uz)))/36.;
         Aq[nr14] = ((beta*nA*nB*(nx - nz))/(nA + nB) + nA*(1 + 3*(ux - uz)))/36.;
         Aq[nr13] = ((beta*nA*nB*(-nx + nz))/(nA + nB) + nA*(1 + 3*(-ux + uz)))/36.;
         Aq[nr16] = ((beta*nA*nB*(ny + nz))/(nA + nB) + nA*(1 + 3*(uy + uz)))/36.;
         Aq[nr15] = ((beta*nA*nB*(-ny - nz))/(nA + nB) + nA*(1 + 3*(-uy - uz)))/36.;
         Aq[nr18] = ((beta*nA*nB*(ny - nz))/(nA + nB) + nA*(1 + 3*(uy - uz)))/36.;
         Aq[nr17] = ((beta*nA*nB*(-ny + nz))/(nA + nB) + nA*(1 + 3*(-uy + uz)))/36.;

         
        
        Bq[n] = nB/3.;
        Bq[nr2] = (-((beta*nA*nB*nx)/(nA + nB)) + nB*(1 + 3*ux))/18.;
        Bq[nr1] = ((beta*nA*nB*nx)/(nA + nB) + nB*(1 - 3*ux))/18.;
        Bq[nr4] = (-((beta*nA*nB*ny)/(nA + nB)) + nB*(1 + 3*uy))/18.;
        Bq[nr3] = ((beta*nA*nB*ny)/(nA + nB) + nB*(1 - 3*uy))/18.;
        Bq[nr6] = (-((beta*nA*nB*nz)/(nA + nB)) + nB*(1 + 3*uz))/18.;
        Bq[nr5] = ((beta*nA*nB*nz)/(nA + nB) + nB*(1 - 3*uz))/18.;
        Bq[nr8] = (-((beta*nA*nB*(nx + ny))/(nA + nB)) + nB*(1 + 3*(ux + uy)))/36.;
        Bq[nr7] = (-((beta*nA*nB*(-nx - ny))/(nA + nB)) + nB*(1 + 3*(-ux - uy)))/36.;
        Bq[nr10] = (-((beta*nA*nB*(nx - ny))/(nA + nB)) + nB*(1 + 3*(ux - uy)))/36.;
        Bq[nr9] = (-((beta*nA*nB*(-nx + ny))/(nA + nB)) + nB*(1 + 3*(-ux + uy)))/36.;
        Bq[nr12] = (-((beta*nA*nB*(nx + nz))/(nA + nB)) + nB*(1 + 3*(ux + uz)))/36.;
        Bq[nr11] = (-((beta*nA*nB*(-nx - nz))/(nA + nB)) + nB*(1 + 3*(-ux - uz)))/36.;
        Bq[nr14] = (-((beta*nA*nB*(nx - nz))/(nA + nB)) + nB*(1 + 3*(ux - uz)))/36.;
        Bq[nr13] = (-((beta*nA*nB*(-nx + nz))/(nA + nB)) + nB*(1 + 3*(-ux + uz)))/36.;
        Bq[nr16] = (-((beta*nA*nB*(ny + nz))/(nA + nB)) + nB*(1 + 3*(uy + uz)))/36.;
        Bq[nr15] = (-((beta*nA*nB*(-ny - nz))/(nA + nB)) + nB*(1 + 3*(-uy - uz)))/36.;
        Bq[nr18] = (-((beta*nA*nB*(ny - nz))/(nA + nB)) + nB*(1 + 3*(uy - uz)))/36.;
        Bq[nr17] = (-((beta*nA*nB*(-ny + nz))/(nA + nB)) + nB*(1 + 3*(-uy + uz)))/36.;
    }
}




extern "C" void ScaLBL_D3Q7_AAodd_PhaseField(int *neighborList, int *Map, double *Aq, double *Bq, 
                                             double *Den, double *Phi, int start, int finish, int Np){
    
    int idx,n,nread;
    double fq,nA,nB;
    
    for (int n=start; n<finish; n++){
        
        //..........Compute the number density for component A............
        // q=0
        fq = Aq[n];
        nA = fq;
        nread = neighborList[n];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+1*Np];
        fq = Aq[nread];
        nA += fq;
        // Autogenerating LBM block...
        nread = neighborList[n+2*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+3*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+4*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+5*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+6*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+7*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+8*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+9*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+10*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+11*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+12*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+13*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+14*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+15*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+16*Np];
        fq = Aq[nread];
        nA += fq;
        nread = neighborList[n+17*Np];
        fq = Aq[nread];
        nA += fq;
        
        
        //..........Compute the number density for component B............
        // q=0
        fq = Bq[n];
        nB = fq;
        nread = neighborList[n];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+1*Np];
        fq = Bq[nread];
        nB += fq;
        // Autogenerating LBM block...
        nread = neighborList[n+2*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+3*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+4*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+5*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+6*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+7*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+8*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+9*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+10*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+11*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+12*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+13*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+14*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+15*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+16*Np];
        fq = Bq[nread];
        nB += fq;
        nread = neighborList[n+17*Np];
        fq = Bq[nread];
        nB += fq;
        
        // save the number densities
        Den[n] = nA;
        Den[Np+n] = nB;
        
        // save the phase indicator field
        idx = Map[n];
        Phi[idx] = (nA-nB)/(nA+nB);
    }
}

extern "C" void ScaLBL_D3Q7_AAeven_PhaseField(int *Map, double *Aq, double *Bq, double *Den, double *Phi, 
                                              int start, int finish, int Np){
    int idx,n,nread;
    double fq,nA,nB;
    for (int n=start; n<finish; n++){
        
        // compute number density for component A
        // q=0
        fq = Aq[n];
        nA = fq;
        fq = Aq[1*Np+n];
        nA += fq;
        // Autogenerating LBM block...
        fq = Aq[2*Np+n];
        nA += fq;
        fq = Aq[3*Np+n];
        nA += fq;
        fq = Aq[4*Np+n];
        nA += fq;
        fq = Aq[5*Np+n];
        nA += fq;
        fq = Aq[6*Np+n];
        nA += fq;
        fq = Aq[7*Np+n];
        nA += fq;
        fq = Aq[8*Np+n];
        nA += fq;
        fq = Aq[9*Np+n];
        nA += fq;
        fq = Aq[10*Np+n];
        nA += fq;
        fq = Aq[11*Np+n];
        nA += fq;
        fq = Aq[12*Np+n];
        nA += fq;
        fq = Aq[13*Np+n];
        nA += fq;
        fq = Aq[14*Np+n];
        nA += fq;
        fq = Aq[15*Np+n];
        nA += fq;
        fq = Aq[16*Np+n];
        nA += fq;
        fq = Aq[17*Np+n];
        nA += fq;
        fq = Aq[18*Np+n];
        nA += fq;    //done...
        
        
        
        // compute number density for component B
        // q=0
        fq = Bq[n];
        nB = fq;
        fq = Bq[1*Np+n];
        nB += fq;
        // Autogenerating LBM block...
        fq = Bq[2*Np+n];
        nB += fq;
        fq = Bq[3*Np+n];
        nB += fq;
        fq = Bq[4*Np+n];
        nB += fq;
        fq = Bq[5*Np+n];
        nB += fq;
        fq = Bq[6*Np+n];
        nB += fq;
        fq = Bq[7*Np+n];
        nB += fq;
        fq = Bq[8*Np+n];
        nB += fq;
        fq = Bq[9*Np+n];
        nB += fq;
        fq = Bq[10*Np+n];
        nB += fq;
        fq = Bq[11*Np+n];
        nB += fq;
        fq = Bq[12*Np+n];
        nB += fq;
        fq = Bq[13*Np+n];
        nB += fq;
        fq = Bq[14*Np+n];
        nB += fq;
        fq = Bq[15*Np+n];
        nB += fq;
        fq = Bq[16*Np+n];
        nB += fq;
        fq = Bq[17*Np+n];
        nB += fq;
        fq = Bq[18*Np+n];
        nB += fq;    //done...
        
        // save the number densities
        Den[n] = nA;
        Den[Np+n] = nB;
        
        // save the phase indicator field
        idx = Map[n];
        Phi[idx] = (nA-nB)/(nA+nB);
    }
}




extern "C" void ScaLBL_D3Q19_Gradient(int *Map, double *phi, double *ColorGrad, int start, int finish, int Np, int Nx, int Ny, int Nz){
    int idx,n,N,i,j,k,nn;
    // distributions
    double f1,f2,f3,f4,f5,f6,f7,f8,f9;
    double f10,f11,f12,f13,f14,f15,f16,f17,f18;
    double nx,ny,nz;
    
    for (idx=0; idx<Np; idx++){
        
        // Get the 1D index based on regular data layout
        n = Map[idx];
        
        //.......Back out the 3D indices for node n..............
        k = n/(Nx*Ny);
        j = (n-Nx*Ny*k)/Nx;
        i = n-Nx*Ny*k-Nx*j;
        //........................................................................
        //........Get 1-D index for this thread....................
        //		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
        //........................................................................
        //					COMPUTE THE COLOR GRADIENT
        //........................................................................
        //.................Read Phase Indicator Values............................
        //........................................................................
        nn = n-1;							// neighbor index (get convention)
        if (i-1<0)		nn += Nx;			// periodic BC along the x-boundary
        f1 = phi[nn];						// get neighbor for phi - 1
        //........................................................................
        nn = n+1;							// neighbor index (get convention)
        if (!(i+1<Nx))	nn -= Nx;			// periodic BC along the x-boundary
        f2 = phi[nn];						// get neighbor for phi - 2
        //........................................................................
        nn = n-Nx;							// neighbor index (get convention)
        if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
        f3 = phi[nn];					// get neighbor for phi - 3
        //........................................................................
        nn = n+Nx;							// neighbor index (get convention)
        if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
        f4 = phi[nn];					// get neighbor for phi - 4
        //........................................................................
        nn = n-Nx*Ny;						// neighbor index (get convention)
        if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
        f5 = phi[nn];					// get neighbor for phi - 5
        //........................................................................
        nn = n+Nx*Ny;						// neighbor index (get convention)
        if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
        f6 = phi[nn];					// get neighbor for phi - 6
        //........................................................................
        nn = n-Nx-1;						// neighbor index (get convention)
        if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
        if (j-1<0)			nn += Nx*Ny;	// Perioidic BC along the y-boundary
        f7 = phi[nn];					// get neighbor for phi - 7
        //........................................................................
        nn = n+Nx+1;						// neighbor index (get convention)
        if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
        if (!(j+1<Ny))		nn -= Nx*Ny;	// Perioidic BC along the y-boundary
        f8 = phi[nn];					// get neighbor for phi - 8
        //........................................................................
        nn = n+Nx-1;						// neighbor index (get convention)
        if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
        if (!(j+1<Ny))		nn -= Nx*Ny;	// Perioidic BC along the y-boundary
        f9 = phi[nn];					// get neighbor for phi - 9
        //........................................................................
        nn = n-Nx+1;						// neighbor index (get convention)
        if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
        if (j-1<0)			nn += Nx*Ny;	// Perioidic BC along the y-boundary
        f10 = phi[nn];					// get neighbor for phi - 10
        //........................................................................
        nn = n-Nx*Ny-1;						// neighbor index (get convention)
        if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
        if (k-1<0)			nn += Nx*Ny*Nz;	// Perioidic BC along the z-boundary
        f11 = phi[nn];					// get neighbor for phi - 11
        //........................................................................
        nn = n+Nx*Ny+1;						// neighbor index (get convention)
        if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
        if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
        f12 = phi[nn];					// get neighbor for phi - 12
        //........................................................................
        nn = n+Nx*Ny-1;						// neighbor index (get convention)
        if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
        if (!(k+1<Nz))		nn -= Nx*Ny*Nz;	// Perioidic BC along the z-boundary
        f13 = phi[nn];					// get neighbor for phi - 13
        //........................................................................
        nn = n-Nx*Ny+1;						// neighbor index (get convention)
        if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
        if (k-1<0)			nn += Nx*Ny*Nz;	// Perioidic BC along the z-boundary
        f14 = phi[nn];					// get neighbor for phi - 14
        //........................................................................
        nn = n-Nx*Ny-Nx;					// neighbor index (get convention)
        if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
        if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
        f15 = phi[nn];					// get neighbor for phi - 15
        //........................................................................
        nn = n+Nx*Ny+Nx;					// neighbor index (get convention)
        if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
        if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
        f16 = phi[nn];					// get neighbor for phi - 16
        //........................................................................
        nn = n+Nx*Ny-Nx;					// neighbor index (get convention)
        if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
        if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
        f17 = phi[nn];					// get neighbor for phi - 17
        //........................................................................
        nn = n-Nx*Ny+Nx;					// neighbor index (get convention)
        if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
        if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
        f18 = phi[nn];					// get neighbor for phi - 18
        //............Compute the Color Gradient...................................
        nx = -(f1-f2+0.5*(f7-f8+f9-f10+f11-f12+f13-f14));
        ny = -(f3-f4+0.5*(f7-f8-f9+f10+f15-f16+f17-f18));
        nz = -(f5-f6+0.5*(f11-f12-f13+f14+f15-f16-f17+f18));
        //...............................................
        //...Store the Color Gradient....................
        ColorGrad[idx] = nx;
        ColorGrad[Np+idx] = ny;
        ColorGrad[2*Np+idx] = nz;
        //...............................................
    }
}

extern "C" void ScaLBL_PhaseField_Init(int *Map, double *Phi, double *DenA, double * DenB, double *Aq, double *Bq, int start, int finish, int Np){
    int idx,n;
    double phi,nA,nB;
    
    for (idx=start; idx<finish; idx++){
        
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
        DenA[idx] = nA;
        DenB[idx] = nB;
         
       
        
        Aq[idx] = 0.3333333333333333*nA;
        Aq[1*Np+idx] = 0.055555555555555555*nA;
        Aq[2*Np+idx] = 0.055555555555555555*nA;
        Aq[3*Np+idx] = 0.055555555555555555*nA;
        Aq[4*Np+idx] = 0.055555555555555555*nA;
        Aq[5*Np+idx] = 0.055555555555555555*nA;
        Aq[6*Np+idx] = 0.055555555555555555*nA;
        
        Aq[7*Np+idx] = 0.0277777777777778*nA;
        Aq[8*Np+idx] = 0.0277777777777778*nA;
        Aq[9*Np+idx] = 0.0277777777777778*nA;
        Aq[10*Np+idx] = 0.0277777777777778*nA;
        Aq[11*Np+idx] = 0.0277777777777778*nA;
        Aq[12*Np+idx] = 0.0277777777777778*nA;
        Aq[13*Np+idx] = 0.0277777777777778*nA;
        Aq[14*Np+idx] = 0.0277777777777778*nA;
        Aq[15*Np+idx] = 0.0277777777777778*nA;
        Aq[16*Np+idx] = 0.0277777777777778*nA;
        Aq[17*Np+idx] = 0.0277777777777778*nA;
        Aq[18*Np+idx] = 0.0277777777777778*nA;
        
        
        
        Bq[idx] = 0.3333333333333333*nB;
        Bq[1*Np+idx] = 0.055555555555555555*nB;
        Bq[2*Np+idx] = 0.055555555555555555*nB;
        Bq[3*Np+idx] = 0.055555555555555555*nB;
        Bq[4*Np+idx] = 0.055555555555555555*nB;
        Bq[5*Np+idx] = 0.055555555555555555*nB;
        Bq[6*Np+idx] = 0.055555555555555555*nB;
        
        Bq[7*Np+idx] = 0.0277777777777778*nB;
        Bq[8*Np+idx] = 0.0277777777777778*nB;
        Bq[9*Np+idx] = 0.0277777777777778*nB;
        Bq[10*Np+idx] = 0.0277777777777778*nB;
        Bq[11*Np+idx] = 0.0277777777777778*nB;
        Bq[12*Np+idx] = 0.0277777777777778*nB;
        Bq[13*Np+idx] = 0.0277777777777778*nB;
        Bq[14*Np+idx] = 0.0277777777777778*nB;
        Bq[15*Np+idx] = 0.0277777777777778*nB;
        Bq[16*Np+idx] = 0.0277777777777778*nB;
        Bq[17*Np+idx] = 0.0277777777777778*nB;
        Bq[18*Np+idx] = 0.0277777777777778*nB;
        
        
        
    }
   
}

extern "C" void ScaLBL_PhaseField_Init_D3Q7(int *Map, double *Phi, double *DenA, double * DenB, double *Aq, double *Bq, int start, int finish, int Np,int N) {
    int idx,n;
    double phi,nA,nB;
    
    for (idx=start; idx<finish; idx++) {
        
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



extern "C" void ScaLBL_PhaseField_Init_D3Q7_2(int *Map, double *Phi, double *DenA, double * DenB, double *laggedDenA, double * laggedDenB, double *analysisDenA, double * analysisDenB, double *Aq, double *Bq, int start, int finish, int Np,int N) {
    int idx,n;
    double phi,nA,nB;
    
    for (idx=start; idx<finish; idx++) {
        
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
        
        laggedDenA[n] = nA;
        laggedDenB[n] = nB;
        
        analysisDenA[n] = nA;
        analysisDenB[n] = nB;
         
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

extern "C" void ScaLBL_Color_InitDistancePacked(int *Map, double *Phi, double *Distance,
                                                double das, double dbs, double beta, double xp, int start, int finish, int Np)
{
    double d;
    int idx,n;
    for (idx=start; idx<finish; idx++){
        n = Map[idx];
        d = Distance[n];
        if (d == -2) Phi[n] = -1; // Set to w fluid
        Phi[n] = (2.0*(exp(-2.0*beta*(d+xp)))/(1.0+exp(-2.0*beta*(d+xp))) - 1.0);
    }
}

//
//\char *ID, double *Den, double *Phi, double das, double dbs, int Nx, int Ny, int Nz)
//{
//    int n,N;
//    
//    N = Nx*Ny*Nz;
//    
//    for (n=0; n<N; n++){
//        
//        if ( ID[n] == 1){
//            Den[n] = 1.0;
//            Den[N+n] = 0.0;
//            Phi[n] = 1.0;
//        }
//        else if ( ID[n] == 2){
//            Den[n] = 0.0;
//            Den[N+n] = 1.0;
//            Phi[n] = -1.0;
//        }

extern "C" void ScaLBL_Color_InitDistanceFull(char *ID, double *Phi, double* Distance, double beta, double xp, int N)
{
    double d;
    double val;
    int idx,n;
    //beta = 1.0;
    for (idx=0; idx<N; idx++){
        d = Distance[idx];
        Phi[idx] = 1.0;
        // Phi[idx] *= (2.0*(exp(-2.0*beta*(d)))/(1.0+exp(-2.0*beta*(d)))-1.0);
        if (ID[idx] == 0) Phi[idx] = -0.2;
        if (ID[idx] == 2) Phi[idx] = -1.0;
        if (ID[idx] == 1) Phi[idx] = 1.0;
        
        //if (d < 0) Phi[idx]*=-0.2;
        
        if (Phi[idx] > 1) Phi[idx] = 1.0;
        if (Phi[idx] < -1) Phi[idx] = -1.0;
        
        // Flip contact angle from > 90 to < 90
        //        Phi[idx]*=-1;
        
        
    }
}





// Adjustable phase affinity
extern "C" double AdjustPhaseAffinity(int *Map, double *Phi, int Np, double contactAngle, char* id){
    int idx,n;
    double phi,nA,nB;
    
    double currentPhaseAffinity;
    
    
    contactAngle += 0.01;
    
    currentPhaseAffinity = contactAngle;
    
    for (idx=0; idx<Np; idx++){
        /* if Solid phase */
        if (id[idx] == 0) {
            
            Phi[idx] = -0.8;  // shouldn't do this for all solid, it's slow.
        }
    }
    
    return currentPhaseAffinity;
    
}



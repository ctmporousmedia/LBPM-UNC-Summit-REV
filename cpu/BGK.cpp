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

#include <iostream>
using namespace std;

extern "C" void ScaLBL_D3Q19_AAeven_BGK(double *dist, int start, int finish, int Np, double tau, double Fx, double Fy, double Fz,double*Velocity){
    // conserved momemnts
    double rho,ux,uy,uz,uu;
    // non-conserved moments
    double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
    
    f0 = 0;
    f1= 0; f2= 0; f3= 0; f4= 0; f5= 0; f6= 0; f7= 0; f8= 0; f9= 0; f10= 0; f11= 0; f12= 0; f13= 0; f14= 0; f15= 0; f16= 0; f17= 0; f18 = 0;
    // cout << "wp=" << wp << " Fx=" << Fx << " Fy=" << Fy << " Fz=" << Fz << endl;
    
    double wp = 1./tau;
    
    double feq0, feq1, feq2, feq3, feq4, feq5, feq6, feq7, feq8, feq9 ,feq10;
       double feq11, feq12, feq13, feq14, feq15, feq16, feq17, feq18;
    
    for (int n=start; n<finish; n++) {
        /* READ OPPOSITE */
        
        f0 = dist[n];
        // q = 1
        f1 = dist[2*Np+n];
        // q = 2
        f2 = dist[1*Np+n];
        f3 = dist[4*Np+n];
        f4 = dist[3*Np+n];
        f5 = dist[6*Np+n];
        f6 = dist[5*Np+n];
        f7 = dist[8*Np+n];
        f8 = dist[7*Np+n];
        f9 = dist[10*Np+n];
        f10 = dist[9*Np+n];
        f11 = dist[12*Np+n];
        f12 = dist[11*Np+n];
        f13 = dist[14*Np+n];
        f14 = dist[13*Np+n];
        f15 = dist[16*Np+n];
        f16 = dist[15*Np+n];
        f17 = dist[18*Np+n];
        f18 = dist[17*Np+n];
        
        rho = f0+f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14+f15+f16+f17+f18;
        ux = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
        uy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
        uz = f5-f6+f11-f12-f13+f14+f15-f16-f17+f18;
        
//        ux += 0.5 * Fx;
//        uy += 0.5 * Fy;
//        uz += 0.5 * Fz;
        
        // write the velocity
        Velocity[n] = ux + 0.5*Fx*rho;
        Velocity[Np+n] = uy + 0.5*Fy*rho;
        Velocity[2*Np+n] = uz + 0.5*Fz*rho;
        
        ux /= rho;
        uy /= rho;
        uz /= rho;
        
        uu = 1.5*(ux*ux+uy*uy+uz*uz);
       
        
//        feq0 =  0.3333333333333333*rho*(1.0-uu);
//        feq1 =  0.05555555555555555*rho*(1 + 3.0*ux + 4.5*ux*ux - uu);
//        feq2 =  0.05555555555555555*rho*(1 - 3.0*ux + 4.5*ux*ux - uu);
//        feq3 =  0.05555555555555555*rho*(1 + 3.0*uy + 4.5*uy*uy - uu);
//        feq4 =  0.05555555555555555*rho*(1 - 3.0*uy + 4.5*uy*uy - uu);
//        feq5 =  0.05555555555555555*rho*(1 + 3.0*uz + 4.5*uz*uz - uu);
//        feq6 =  0.05555555555555555*rho*(1 - 3.0*uz + 4.5*uz*uz - uu);
//        feq7 =  0.02777777777777778*rho*(1 + 3.0*(ux+uy)  + 4.5*(ux+uy)*(ux+uy) - uu);
//        feq8 =  0.02777777777777778*rho*(1 + 3.0*(-ux-uy) + 4.5*(-ux-uy)*(-ux-uy) - uu);
//        feq9 =  0.0277777777777778*rho*(1  + 3.0*(ux-uy)  + 4.5*(ux-uy)*(ux-uy) - uu);
//        feq10 = 0.0277777777777778*rho*(1  + 3.0*(-ux+uy) + 4.5*(-ux+uy)*(-ux+uy) - uu);
//        feq11 = 0.0277777777777778*rho*(1  + 3.0*(ux+uz)  + 4.5*(ux+uz)*(ux+uz) - uu);
//        feq12 = 0.0277777777777778*rho*(1  + 3.0*(-ux-uz) + 4.5*(-ux-uz)*(-ux-uz) - uu);
//        feq13 = 0.0277777777777778*rho*(1  + 3.0*(ux-uz)  + 4.5*(ux-uz)*(ux-uz) - uu);
//        feq14 = 0.0277777777777778*rho*(1  + 3.0*(-ux+uz) + 4.5*(-ux+uz)*(-ux+uz) - uu);
//        feq15 = 0.0277777777777778*rho*(1  + 3.0*(uy+uz)  + 4.5*(uy+uz)*(uy+uz) - uu);
//        feq16 = 0.0277777777777778*rho*(1  + 3.0*(-uy-uz) + 4.5*(-uy-uz)*(-uy-uz) - uu);
//        feq17 = 0.0277777777777778*rho*(1  + 3.0*(uy-uz)  + 4.5*(uy-uz)*(uy-uz) - uu);
//        feq18 = 0.0277777777777778*rho*(1  + 3.0*(-uy+uz) + 4.5*(-uy+uz)*(-uy+uz) - uu);
        
        feq0 =  0.3333333333333333*rho*1.0; // -uu);
        feq1 =  0.05555555555555555*rho*(1 + 3.0*ux); //  + 4.5*ux*ux - uu);
        feq2 =  0.05555555555555555*rho*(1 - 3.0*ux); //+ 4.5*ux*ux - uu);
        feq3 =  0.05555555555555555*rho*(1 + 3.0*uy); //+ 4.5*uy*uy - uu);
        feq4 =  0.05555555555555555*rho*(1 - 3.0*uy); //+ 4.5*uy*uy - uu);
        feq5 =  0.05555555555555555*rho*(1 + 3.0*uz); //+ 4.5*uz*uz - uu);
        feq6 =  0.05555555555555555*rho*(1 - 3.0*uz); //+ 4.5*uz*uz - uu);
        feq7 =  0.02777777777777778*rho*(1 + 3.0*(ux+uy) ); //+ 4.5*(ux+uy)*(ux+uy) - uu);
        feq8 =  0.02777777777777778*rho*(1 + 3.0*(-ux-uy)); //+ 4.5*(-ux-uy)*(-ux-uy) - uu);
        feq9 =  0.0277777777777778*rho*(1  + 3.0*(ux-uy) ); //+ 4.5*(ux-uy)*(ux-uy) - uu);
        feq10 = 0.0277777777777778*rho*(1  + 3.0*(-ux+uy)); //+ 4.5*(-ux+uy)*(-ux+uy) - uu);
        feq11 = 0.0277777777777778*rho*(1  + 3.0*(ux+uz) ); //+ 4.5*(ux+uz)*(ux+uz) - uu);
        feq12 = 0.0277777777777778*rho*(1  + 3.0*(-ux-uz)); //+ 4.5*(-ux-uz)*(-ux-uz) - uu);
        feq13 = 0.0277777777777778*rho*(1  + 3.0*(ux-uz) ); //+ 4.5*(ux-uz)*(ux-uz) - uu);
        feq14 = 0.0277777777777778*rho*(1  + 3.0*(-ux+uz)); //+ 4.5*(-ux+uz)*(-ux+uz) - uu);
        feq15 = 0.0277777777777778*rho*(1  + 3.0*(uy+uz) ); //+ 4.5*(uy+uz)*(uy+uz) - uu);
        feq16 = 0.0277777777777778*rho*(1  + 3.0*(-uy-uz)); //+ 4.5*(-uy-uz)*(-uy-uz) - uu);
        feq17 = 0.0277777777777778*rho*(1  + 3.0*(uy-uz) ); //+ 4.5*(uy-uz)*(uy-uz) - uu);
        feq18 = 0.0277777777777778*rho*(1  + 3.0*(-uy+uz)); //+ 4.5*(-uy+uz)*(-uy+uz) - uu);
        
        

        
        // q=0
        dist[n] = f0*(1.0-wp)+wp*feq0;
        
        // q = 1
        dist[1*Np+n] = f1*(1.0-wp) + wp*feq1 + (1.-0.0*wp)*0.166666666666667*Fx;
        
        // q = 2
        dist[2*Np+n] = f2*(1.0-wp) + wp*feq2 - (1.-0.0*wp)*0.166666666666667*Fx;
        
        // q = 3
        dist[3*Np+n] = f3*(1.0-wp) + wp*feq3 + (1.-0.0*wp)*0.166666666666667*Fy;
        
        // q = 4
        dist[4*Np+n] = f4*(1.0-wp) + wp*feq4 - (1.-0.0*wp)*0.166666666666667*Fy;
        
        // q = 5
        dist[5*Np+n] = f5*(1.0-wp) + wp*feq5 + (1.-0.0*wp)*0.166666666666667*Fz;
        
        // q = 6
        dist[6*Np+n] = f6*(1.0-wp) + wp*feq6 - (1.-0.0*wp)*0.166666666666667*Fz;
        
        // q = 7
        dist[7*Np+n] = f7*(1.0-wp) + wp*feq7    + (1.-0.0*wp)*0.083333333333333333333*(Fx+Fy);
        
        // q = 8
        dist[8*Np+n] = f8*(1.0-wp) + wp*feq8  + (1.-0.0*wp)*0.083333333333333333333*(-Fx-Fy);

        // q = 9
        dist[9*Np+n] = f9*(1.0-wp) + wp*feq9    + (1.-0.0*wp)*0.083333333333333333333*(Fx-Fy);

        // q = 10
        dist[10*Np+n] = f10*(1.0-wp) + wp*feq10 + (1.-0.0*wp)*0.083333333333333333333*(-Fx+Fy);

        // q = 11
        dist[11*Np+n] = f11*(1.0-wp) + wp*feq11   + (1.-0.0*wp)*0.083333333333333333333*(Fx+Fz);

        // q = 12
        dist[12*Np+n] = f12*(1.0-wp) + wp*feq12  + (1.-0.0*wp)*0.083333333333333333333*(-Fx-Fz);

        // q = 13
        dist[13*Np+n] = f13*(1.0-wp) + wp*feq13   + (1.-0.0*wp)*0.083333333333333333333*(Fx-Fz);

        // q= 14
        dist[14*Np+n] = f14*(1.0-wp) + wp*feq14  + (1.-0.0*wp)*0.083333333333333333333*(-Fx+Fz);

        // q = 15
        dist[15*Np+n] = f15*(1.0-wp) + wp*feq15    + (1.-0.0*wp)*0.083333333333333333333*(Fy+Fz);

        // q = 16
        dist[16*Np+n] = f16*(1.0-wp) + wp*feq16  + (1.-0.0*wp)*0.083333333333333333333*(-Fy-Fz);

        // q = 17
        dist[17*Np+n] = f17*(1.0-wp) + wp*feq17    + (1.-0.0*wp)*0.083333333333333333333*(Fy-Fz);

        // q = 18
        dist[18*Np+n] = f18*(1.0-wp) + wp*feq18  + (1.-0.0*wp)*0.083333333333333333333*(-Fy+Fz);
        
        
    }
}

extern "C" void ScaLBL_D3Q19_AAodd_BGK(int *neighborList, double *dist, int start, int finish, int Np, double tau, double Fx, double Fy, double Fz,double*Velocity){
    // conserved momemnts
    double rho,ux,uy,uz,uu;
    // non-conserved moments
    double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
    int nr1,nr2,nr3,nr4,nr5,nr6,nr7,nr8,nr9,nr10,nr11,nr12,nr13,nr14,nr15,nr16,nr17,nr18;
    
    double feq0, feq1, feq2, feq3, feq4, feq5, feq6, feq7, feq8, feq9 ,feq10;
    double feq11, feq12, feq13, feq14, feq15, feq16, feq17, feq18;
    
    f0= 0; f1= 0; f2= 0; f3= 0; f4= 0; f5= 0; f6= 0; f7= 0; f8= 0; f9= 0; f10= 0; f11= 0; f12= 0; f13= 0; f14= 0; f15= 0; f16= 0; f17= 0; f18= 0;
     // cout << "Np=" << Np << " wp=" << wp << " Fx=" << Fx << " Fy=" << Fy << " Fz=" << Fz << " start=" << start << " finish=" << finish <<  endl;
    
    double wp = 1./tau;
    
    
    for (int n=start; n<finish; n++){
        
        // q=0
        f0 = dist[n];
        
        // q = 1
        nr1 = neighborList[n];
        f1 = dist[nr1];
        
        // q = 2
        nr2 = neighborList[n+Np];
        f2 = dist[nr2];
        
        // q = 3
        nr3 = neighborList[n+2*Np];
        f3 = dist[nr3];
        
        // q = 4
        nr4 = neighborList[n+3*Np];
        f4 = dist[nr4];
        
        // q = 5
        nr5 = neighborList[n+4*Np];
        f5 = dist[nr5];
        
        // q = 6
        nr6 = neighborList[n+5*Np];
        f6 = dist[nr6];
        
        // q = 7
        nr7 = neighborList[n+6*Np];
        f7 = dist[nr7];

        // q = 8
        nr8 = neighborList[n+7*Np];
        f8 = dist[nr8];

        // q = 9
        nr9 = neighborList[n+8*Np];
        f9 = dist[nr9];

        // q = 10
        nr10 = neighborList[n+9*Np];
        f10 = dist[nr10];

        // q=11
        nr11 = neighborList[n+10*Np];
        f11 = dist[nr11];

        // q=12
        nr12 = neighborList[n+11*Np];
        f12 = dist[nr12];

        // q=13
        nr13 = neighborList[n+12*Np];
        f13 = dist[nr13];

        // q=14
        nr14 = neighborList[n+13*Np];
        f14 = dist[nr14];

        // q=15
        nr15 = neighborList[n+14*Np];
        f15 = dist[nr15];

        // q=16
        nr16 = neighborList[n+15*Np];
        f16 = dist[nr16];

        // q=17
        nr17 = neighborList[n+16*Np];
        f17 = dist[nr17];

        // q=18
        nr18 = neighborList[n+17*Np];
        f18 = dist[nr18];
        
        rho = f0+f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14+f15+f16+f17+f18;
        ux = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
        uy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
        uz = f5-f6+f11-f12-f13+f14+f15-f16-f17+f18;
        
//        ux +=  Fx;
//        uy +=  Fy;
//        uz +=  Fz;
        ux /= rho;
        uy /= rho;
        uz /= rho;

        
//        uu = 1.5*(ux*ux+uy*uy+uz*uz);
//        feq0 =  0.3333333333333333*rho*(1.0-uu);
//        feq1 =  0.05555555555555555*rho*(1 + 3.0*ux + 4.5*ux*ux - uu);
//        feq2 =  0.05555555555555555*rho*(1 - 3.0*ux + 4.5*ux*ux - uu);
//        feq3 =  0.05555555555555555*rho*(1 + 3.0*uy + 4.5*uy*uy - uu);
//        feq4 =  0.05555555555555555*rho*(1 - 3.0*uy + 4.5*uy*uy - uu);
//        feq5 =  0.05555555555555555*rho*(1 + 3.0*uz + 4.5*uz*uz - uu);
//        feq6 =  0.05555555555555555*rho*(1 - 3.0*uz + 4.5*uz*uz - uu);
//        feq7 =  0.02777777777777778*rho*(1 + 3.0*(ux+uy)  + 4.5*(ux+uy)*(ux+uy) - uu);
//        feq8 =  0.02777777777777778*rho*(1 + 3.0*(-ux-uy) + 4.5*(-ux-uy)*(-ux-uy) - uu);
//        feq9 =  0.0277777777777778*rho*(1  + 3.0*(ux-uy)  + 4.5*(ux-uy)*(ux-uy) - uu);
//        feq10 = 0.0277777777777778*rho*(1  + 3.0*(-ux+uy) + 4.5*(-ux+uy)*(-ux+uy) - uu);
//        feq11 = 0.0277777777777778*rho*(1  + 3.0*(ux+uz)  + 4.5*(ux+uz)*(ux+uz) - uu);
//        feq12 = 0.0277777777777778*rho*(1  + 3.0*(-ux-uz) + 4.5*(-ux-uz)*(-ux-uz) - uu);
//        feq13 = 0.0277777777777778*rho*(1  + 3.0*(ux-uz)  + 4.5*(ux-uz)*(ux-uz) - uu);
//        feq14 = 0.0277777777777778*rho*(1  + 3.0*(-ux+uz) + 4.5*(-ux+uz)*(-ux+uz) - uu);
//        feq15 = 0.0277777777777778*rho*(1  + 3.0*(uy+uz)  + 4.5*(uy+uz)*(uy+uz) - uu);
//        feq16 = 0.0277777777777778*rho*(1  + 3.0*(-uy-uz) + 4.5*(-uy-uz)*(-uy-uz) - uu);
//        feq17 = 0.0277777777777778*rho*(1  + 3.0*(uy-uz)  + 4.5*(uy-uz)*(uy-uz) - uu);
//        feq18 = 0.0277777777777778*rho*(1  + 3.0*(-uy+uz) + 4.5*(-uy+uz)*(-uy+uz) - uu);
        
        feq0 =  0.3333333333333333*rho*1.0; // -uu);
        feq1 =  0.05555555555555555*rho*(1 + 3.0*ux); //  + 4.5*ux*ux - uu);
        feq2 =  0.05555555555555555*rho*(1 - 3.0*ux); //+ 4.5*ux*ux - uu);
        feq3 =  0.05555555555555555*rho*(1 + 3.0*uy); //+ 4.5*uy*uy - uu);
        feq4 =  0.05555555555555555*rho*(1 - 3.0*uy); //+ 4.5*uy*uy - uu);
        feq5 =  0.05555555555555555*rho*(1 + 3.0*uz); //+ 4.5*uz*uz - uu);
        feq6 =  0.05555555555555555*rho*(1 - 3.0*uz); //+ 4.5*uz*uz - uu);
        feq7 =  0.02777777777777778*rho*(1 + 3.0*(ux+uy) ); //+ 4.5*(ux+uy)*(ux+uy) - uu);
        feq8 =  0.02777777777777778*rho*(1 + 3.0*(-ux-uy)); //+ 4.5*(-ux-uy)*(-ux-uy) - uu);
        feq9 =  0.0277777777777778*rho*(1  + 3.0*(ux-uy) ); //+ 4.5*(ux-uy)*(ux-uy) - uu);
        feq10 = 0.0277777777777778*rho*(1  + 3.0*(-ux+uy)); //+ 4.5*(-ux+uy)*(-ux+uy) - uu);
        feq11 = 0.0277777777777778*rho*(1  + 3.0*(ux+uz) ); //+ 4.5*(ux+uz)*(ux+uz) - uu);
        feq12 = 0.0277777777777778*rho*(1  + 3.0*(-ux-uz)); //+ 4.5*(-ux-uz)*(-ux-uz) - uu);
        feq13 = 0.0277777777777778*rho*(1  + 3.0*(ux-uz) ); //+ 4.5*(ux-uz)*(ux-uz) - uu);
        feq14 = 0.0277777777777778*rho*(1  + 3.0*(-ux+uz)); //+ 4.5*(-ux+uz)*(-ux+uz) - uu);
        feq15 = 0.0277777777777778*rho*(1  + 3.0*(uy+uz) ); //+ 4.5*(uy+uz)*(uy+uz) - uu);
        feq16 = 0.0277777777777778*rho*(1  + 3.0*(-uy-uz)); //+ 4.5*(-uy-uz)*(-uy-uz) - uu);
        feq17 = 0.0277777777777778*rho*(1  + 3.0*(uy-uz) ); //+ 4.5*(uy-uz)*(uy-uz) - uu);
        feq18 = 0.0277777777777778*rho*(1  + 3.0*(-uy+uz)); //+ 4.5*(-uy+uz)*(-uy+uz) - uu);
        

        // q=0
        dist[n] = f0*(1.0-wp)+wp*feq0;
        
        // q = 1
        dist[nr2] = f1*(1.0-wp) + wp*feq1        + (1.-0.0*wp)*0.166666666666667*Fx;
        
        // q = 2
        dist[nr1] = f2*(1.0-wp) + wp*feq2        - (1.-0.0*wp)*0.166666666666667*Fx;
        
        // q = 3
        dist[nr4] = f3*(1.0-wp) + wp*feq3        + (1.-0.0*wp)*0.166666666666667*Fy;
        
        // q = 4
        dist[nr3] = f4*(1.0-wp) + wp*feq4       - (1.-0.0*wp)*0.166666666666667*Fy;
        
        // q = 5
        dist[nr6] = f5*(1.0-wp) + wp*feq5       + (1.-0.0*wp)*0.166666666666667*Fz;
        
        // q = 6
        dist[nr5] = f6*(1.0-wp) + wp*feq6       - (1.-0.0*wp)*0.166666666666667*Fz;
        
        // q = 7
        dist[nr8] = f7*(1.0-wp) + wp*feq7       + (1.-0.0*wp)*0.083333333333333333333*(Fx+Fy);

        // q = 8
        dist[nr7] = f8*(1.0-wp) + wp*feq8       + (1.-0.0*wp)*0.083333333333333333333*(-Fx-Fy);

        // q = 9
        dist[nr10] = f9*(1.0-wp) + wp*feq9      + (1.-0.0*wp)*0.083333333333333333333*(Fx-Fy);

        // q = 10
        dist[nr9] = f10*(1.0-wp) + wp*feq10     + (1.-0.0*wp)*0.083333333333333333333*(-Fx+Fy);

        // q = 11
        dist[nr12] = f11*(1.0-wp) + wp*feq11    + (1.-0.0*wp)*0.083333333333333333333*(Fx+Fz);

        // q = 12
        dist[nr11] = f12*(1.0-wp) + wp*feq12    + (1.-0.0*wp)*0.083333333333333333333*(-Fx-Fz);

        // q = 13
        dist[nr14] = f13*(1.0-wp) + wp*feq13    + (1.-0.0*wp)*0.083333333333333333333*(Fx-Fz);

        // q= 14
        dist[nr13] = f14*(1.0-wp) + wp*feq14    + (1.-0.0*wp)*0.083333333333333333333*(-Fx+Fz);

        // q = 15
        dist[nr16] = f15*(1.0-wp) + wp*feq15    + (1.-0.0*wp)*0.083333333333333333333*(Fy+Fz);

        // q = 16
        dist[nr15] = f16*(1.0-wp) + wp*feq16    + (1.-0.0*wp)*0.083333333333333333333*(-Fy-Fz);

        // q = 17
        dist[nr18] = f17*(1.0-wp) + wp*feq17    + (1.-0.0*wp)*0.083333333333333333333*(Fy-Fz);

        // q = 18
        dist[nr17] = f18*(1.0-wp) + wp*feq18    + (1.-0.0*wp)*0.083333333333333333333*(-Fy+Fz);
      
        
    }
}




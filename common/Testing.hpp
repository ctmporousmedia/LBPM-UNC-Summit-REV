//
//  Testing.hpp
//  
//
//  Created by Christopher Fowler on 10/29/19.
//

#ifndef Testing_hpp
#define Testing_hpp


#include <iostream>
#include <fstream>
#include <assert.h>
#include <stdio.h>
#include <vector>

void save_3Ddouble_serial(double* field, int Nx, int Ny, int Nz,int depth, std::string s)
{
    assert(depth >= 0);
    char Filename[40];
    char cstr[s.size() + 1];
    s.copy(cstr, s.size() + 1);
    cstr[s.size()] = '\0';
    sprintf(Filename,"%s.txt",cstr);
    FILE *f2 = fopen(Filename, "w");
    if (f2 == NULL) { printf("Error opening file!\n");  exit(1);  }

    double var = 0;
    for (int z =depth; z < Nz-depth; z++) {
        for (int y =depth; y < Ny-depth; y++) {
            for (int x=depth; x < Nx-depth; x++) {
                int n = x + Nx*y + Nx*Ny*z;
                var = field[n];
                fprintf(f2, "%d %d %d %.0f\n",x, y, z, var);
                
            } // i // iproc //    printf("\n");
        } // j // jproc //  printf("\n");
    }  // k // kproc
    fclose(f2);
    
    // printf("LINE=%d\n",__LINE__);
    //    delete[] rankArrayNoPM;
}

void save_3Dint_serial(int* field, int Nx, int Ny, int Nz,int depth, std::string s)
{
    assert(depth >= 0);
    char Filename[40];
    char cstr[s.size() + 1];
    s.copy(cstr, s.size() + 1);
    cstr[s.size()] = '\0';
    sprintf(Filename,"%s.txt",cstr);
    FILE *f2 = fopen(Filename, "w");
    if (f2 == NULL) { printf("Error opening file!\n");  exit(1);  }

    int var = 0;
    for (int z =depth; z < Nz-depth; z++) {
        for (int y =depth; y < Ny-depth; y++) {
            for (int x=depth; x < Nx-depth; x++) {
                int n = x + Nx*y + Nx*Ny*z;
                var = field[n];
                fprintf(f2, "%d %d %d %d\n",x, y, z, var);
                
            } // i // iproc //    printf("\n");
        } // j // jproc //  printf("\n");
    }  // k // kproc
    fclose(f2);
    
    // printf("LINE=%d\n",__LINE__);
    //    delete[] rankArrayNoPM;
}

void save_4Dint_serial(int* field, int Nx, int Ny, int Nz, int slice, int depth, std::string s)
{
    assert(depth >= 0);
    char Filename[40];
    char cstr[s.size() + 1];
    s.copy(cstr, s.size() + 1);
    cstr[s.size()] = '\0';
    sprintf(Filename,"%s.txt",cstr);
    FILE *f2 = fopen(Filename, "w");
    if (f2 == NULL) { printf("Error opening file!\n");  exit(1);  }

    int var = 0;
    
    for (int z =depth; z < Nz-depth; z++) {
        for (int y =depth; y < Ny-depth; y++) {
            for (int x=depth; x < Nx-depth; x++) {
                int n = x + Nx*y + Nx*Ny*z + slice*Nx*Ny*Nz;
                var = field[n];
                fprintf(f2, "%d %d %d %d\n",x, y, z, var);
                
            } // i // iproc //    printf("\n");
        } // j // jproc //  printf("\n");
    }  // k // kproc
    fclose(f2);
    
    // printf("LINE=%d\n",__LINE__);
    //    delete[] rankArrayNoPM;
}

void save_4Dsize_t_serial(size_t* field, int Nx, int Ny, int Nz, int slice, int depth, std::string s)
{
    assert(depth >= 0);
    char Filename[40];
    char cstr[s.size() + 1];
    s.copy(cstr, s.size() + 1);
    cstr[s.size()] = '\0';
    sprintf(Filename,"%s.txt",cstr);
    FILE *f2 = fopen(Filename, "w");
    if (f2 == NULL) { printf("Error opening file!\n");  exit(1);  }

    size_t var = 0;
    
    for (int z =depth; z < Nz-depth; z++) {
        for (int y =depth; y < Ny-depth; y++) {
            for (int x=depth; x < Nx-depth; x++) {
                int n = x + Nx*y + Nx*Ny*z + slice*Nx*Ny*Nz;
                var = field[n];
                fprintf(f2, "%d %d %d %d\n",x, y, z, var);
                
            } // i // iproc //    printf("\n");
        } // j // jproc //  printf("\n");
    }  // k // kproc
    fclose(f2);
    
    // printf("LINE=%d\n",__LINE__);
    //    delete[] rankArrayNoPM;
}

void save_3Ddouble_parallel(int rank, int* field, int Nx, int Ny, int Nz, int depth, std::string s)
{
    assert(depth >= 0);
    char cstr[s.size() + 1];
    s.copy(cstr, s.size() + 1);
    cstr[s.size()] = '\0';
    char LocalRankString[8];
    char LocalRankFilename[40];
    sprintf(LocalRankString,".%05d",rank);
    sprintf(LocalRankFilename,"%s%s.txt",cstr,LocalRankString);
    FILE *f2 = fopen(LocalRankFilename, "w");

    double var = 0;
    for (int z =depth; z < Nz-depth; z++) {
        for (int y =depth; y < Ny-depth; y++) {
            for (int x=depth; x < Nx-depth; x++) {
                int n = x + Nx*y + Nx*Ny*z;
                var = field[n];
                fprintf(f2, "%d %d %d %.0f\n",x, y, z, var);
                
            } // i // iproc //    printf("\n");
        } // j // jproc //  printf("\n");
    }  // k // kproc
    fclose(f2);
    
    // printf("LINE=%d\n",__LINE__);
    //    delete[] rankArrayNoPM;
}

void save_3Dint_parallel(int rank, int* field, int Nx, int Ny, int Nz, int depth, std::string s)
{
    
    assert(depth >= 0);
    char cstr[s.size() + 1];
    s.copy(cstr, s.size() + 1);
    cstr[s.size()] = '\0';
    char LocalRankString[8];
    char LocalRankFilename[40];
    sprintf(LocalRankString,".%05d",rank);
    sprintf(LocalRankFilename,"%s%s.txt",cstr,LocalRankString);
    FILE *f2 = fopen(LocalRankFilename, "w");

    int var = 0;
    for (int z =depth; z < Nz-depth; z++) {
        for (int y =depth; y < Ny-depth; y++) {
            for (int x=depth; x < Nx-depth; x++) {
                int  n = x + Nx*y + Nx*Ny*z;
                var = field[n];
                fprintf(f2, "%d %d %d %d\n",x, y, z, var);
                
            } // i // iproc //    printf("\n");
        } // j // jproc //  printf("\n");
    }  // k // kproc
    fclose(f2);
    
    // printf("LINE=%d\n",__LINE__);
    //    delete[] rankArrayNoPM;
}

void CALCULATE_FACE_POROSITY(int rank, int* field, int Nx, int Ny, int Nz, int depth, int edges, std::string s)
{
    
    // FACE RUNS FROM 1 TO 6
    assert(depth >= 0);
    
    char cstr[s.size() + 1];
    s.copy(cstr, s.size() + 1);
    cstr[s.size()] = '\0';
    char LocalRankString[8];
    char LocalRankFilename[40];
    sprintf(LocalRankString,".%05d",rank);
    sprintf(LocalRankFilename,"%s%s.txt",cstr,LocalRankString);
    
    
    int x_face_sum = 0;
    int X_face_sum = 0;
    int y_face_sum = 0;
    int Y_face_sum = 0;
    int z_face_sum = 0;
    int Z_face_sum = 0;
    
    
    // CREATE FACE ELEMENTS ONLY ARRAY
    int * bdry_array; bdry_array = new int[Nx*Ny*Nz]; for(int i=0; i <Nx*Ny*Nz; i++) {bdry_array[i] = 1; }
    for (int z =depth+1; z < Nz-depth-1; z++) {
        for (int y =depth+1; y < Ny-depth-1; y++) {
            for (int x=depth+1; x < Nx-depth-1; x++) {
                int  n = x + Nx*y + Nx*Ny*z;
                bdry_array[n] = 0;
            } // i // iproc //    printf("\n");
        } // j // jproc //  printf("\n");
    }  // k // kproc
    for (int z =depth; z < Nz-depth; z++) {
        for (int y =depth; y < Ny-depth; y++) {
            for (int x=depth; x < Nx-depth; x++) {
                int  n = x + Nx*y + Nx*Ny*z;
                if (edges == 0) {
                    if (x ==depth && y ==depth) bdry_array[n] = 0;
                    if (x ==(Nx-depth-1) && y ==depth ) bdry_array[n] = 0;
                    if (x ==depth && y ==(Ny-depth-1)) bdry_array[n] = 0;
                    if (x ==(Nx-depth-1) && y ==(Ny-depth-1)) bdry_array[n] = 0;
                    
                    if ( y ==depth && z ==depth) bdry_array[n] = 0;
                    if (y ==(Ny-depth-1) && z ==depth) bdry_array[n] = 0;
                    if ( y ==depth && z ==(Nz-depth-1)) bdry_array[n] = 0;
                    if(  y ==(Ny-depth-1) && z ==(Nz-depth-1)) bdry_array[n] = 0;
                    
                    if (x ==depth && z ==depth ) bdry_array[n] = 0;
                    if (x ==depth &&  z ==(Nz-depth-1)) bdry_array[n] = 0;
                    if (x ==(Nx-depth-1) && z ==depth) bdry_array[n] = 0;
                    if (x ==(Nx-depth-1) && z ==(Nz-depth-1)) bdry_array[n] = 0;
                }
                
                
            } // i // iproc //    printf("\n");
        } // j // jproc //  printf("\n");
    }  // k // kproc
    
  //  save_3Dint_parallel(rank,bdry_array,Nx,Ny,Nz,0,"bdry_array");
    
    for (int z =depth; z < Nz-depth; z++) {
        for (int y =depth; y < Ny-depth; y++) {
            for (int x=depth; x < Nx-depth; x++) {
                int  n = x + Nx*y + Nx*Ny*z;
                if (bdry_array[n] == 1) {
                    
                    if (x == (Nx-depth-1) && field[n] <= 0) {
                        X_face_sum++;
                    }
                    if (x == (depth) && field[n] <= 0) {
                        x_face_sum++;
                    }
                    if (y == (Ny-depth-1) && field[n] <= 0) {
                        Y_face_sum++;
                    }
                    if (y == (depth) && field[n] <= 0) {
                        y_face_sum++;
                    }
                    if (z == (Nz-depth-1) && field[n] <= 0) {
                        Z_face_sum++;
                    }
                    if (z == (depth) && field[n] <= 0) {
                        z_face_sum++;
                    }
                    
                    
                }
            } // i // iproc //    printf("\n");
        } // j // jproc //  printf("\n");
    }

printf("(%d) %s face porosity: X(1)=%d x(2)=%d Y(3)=%d y(4)=%d Z(5)=%d z(6)=%d\n", rank,LocalRankFilename, X_face_sum,x_face_sum,Y_face_sum,y_face_sum,Z_face_sum,z_face_sum);

// printf("LINE=%d\n",__LINE__);
//    delete[] rankArrayNoPM;
delete[] bdry_array;
}


void save_ptrArrayInt_parallel(int rank, std::vector<int> countVec, int** field, int count, std::string s)
{
    char cstr[s.size() + 1];
    s.copy(cstr, s.size() + 1);
    cstr[s.size()] = '\0';
    char LocalRankString[8];
    char LocalRankFilename[40];
    sprintf(LocalRankString,".%05d",rank);
    sprintf(LocalRankFilename,"%s%s.txt",cstr,LocalRankString);
    FILE *f2 = fopen(LocalRankFilename, "w");
    if (f2 == NULL) { printf("Error opening file!\n");  exit(1);  }
    for (int i =0; i < count; i++) {
        for (int j =0; j < countVec.at(i); j++) {
            int ranpm = field[i][j];
            fprintf(f2, "%d %d 1 %d\n",i, j, ranpm);
        } // i // jproc //  printf("\n");
    }  // j // kproc
    fclose(f2);
}



#endif /* Testing_hpp */

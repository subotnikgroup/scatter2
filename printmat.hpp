///////////////////////////////////////////////////////////////////////////////
//
// scatter2: 2D scattering calculation Library
// Copyright (C) 2020 Joseph Subotnik's group, The University of Pennsylvania
//
// This file is part of scatter2.
//
// scatter2 is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// scatter2 is distributed in the hope that it will be usefull, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with scatter2. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////


/* C-Style print of armadillo matrices */

#include <stdio.h>

template<typename T>
void fprint_mat(const T& mat, FILE* fp=stdout, const char* name="", bool longprint=true){

    char buffer[32];

    if (strlen(name) > 0){
        fprintf(fp, "%s\n", name);
    }

    const char* formatter1 = longprint ? "% .8g" : " % .4g";
    const char* formatter2 = longprint ? " %15s" : "  %11s";

    for(int i = 0; i < mat.n_rows;i++){
        for(int j = 0; j < mat.n_cols; j++){
            memset(buffer, 0, 32);
            sprintf(buffer, formatter1, mat(i,j));
            fprintf(fp, formatter2, buffer);
        }
        fprintf(fp, "\n");
    }
}



// For armadillo matrix only
template<typename T>
void fprint_cmat(const T& mat, FILE* fp=stdout, const char* name="", bool longprint=true){

    char buffer[64];

    if (strlen(name) > 0){
        fprintf(fp, "%s\n", name);
    }

    const char* formatter1 = longprint ? "% .8g%+.8gj" : "% .4g%+.4gj";
    const char* formatter2 = longprint ? " %31s" : " %19s";

    for(int i = 0; i < mat.n_rows;i++){
        for(int j = 0; j < mat.n_cols; j++){
            memset(buffer, 0, 64);
            sprintf(buffer, formatter1, mat(i,j).real(), mat(i,j).imag());
            fprintf(fp, formatter2, buffer);
        }
        fprintf(fp, "\n");
    }
}

#define print_cmat(mat) fprint_cmat((mat), stdout, #mat, false)

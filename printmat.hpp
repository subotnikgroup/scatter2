
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

#define ARMA_USE_SUPERLU
#ifndef ARMA_SUPERLU_INCLUDE_DIR
#define ARMA_SUPERLU_INCLUDE_DIR /opt/superlu/gnu/5.2.1/include/ 
#endif

#include <armadillo>
#include <array>
#include <string>
#include <numeric>
#include "external/cxxopts.hpp"
#include "printmat.hpp"

#include "potential.h"
#include "angularpot.hpp"
#include "linearpot.hpp"

using cx_mat = arma::cx_mat;
using mat = arma::mat;
using cx_vec = arma::cx_vec;
using vec = arma::vec;
using cx_double = arma::cx_double;

const std::complex<double> I(0, 1);

#define STR(x) std::to_string(x)

#define arange(a, b) (arma::linspace<vec>(a, b, ((b)-(a)+1)))
#define to_complex(X) (arma::conv_to<cx_mat>::from(X))
#define RSUBMAT(m,rbegin,rend,cbegin,cend) (m.submat((rbegin),(cbegin),(rend),(cend)))

// for armadillo matrices only; calculate submat for two sections
// includes end.
template<typename T>
void submat_2seg(const T& m, T& m2, int begin1, int end1, int begin2, int end2){
    m2 = T(end1-begin1+end2-begin2+2, end1-begin1+end2-begin2+2);
    
    int n1 = end1-begin1;
    int n2 = end2-begin2;

    RSUBMAT(m2, 0, n1, 0, n1) = RSUBMAT(m, begin1, end1, begin1, end1);
    RSUBMAT(m2, 0, n1, n1+1, n1+n2+1) = RSUBMAT(m, begin1, end1, begin2, end2);
    RSUBMAT(m2, n1+1, n1+n2+1, 0, n1) = RSUBMAT(m, begin2, end2, begin1, end1);
    RSUBMAT(m2, n1+1, n1+n2+1, n1+1, n1+n2+1) = RSUBMAT(m, begin2, end2, begin2, end2);

}

template<typename T>
void subvec_2seg(const T& m, T& m2, int begin1, int end1, int begin2, int end2){
    m2 = T(end1-begin1+end2-begin2+2, 1);

    int n1 = end1 - begin1;
    int n2 = end2 - begin2;

    m2.submat(0, 0, n1, 0) = m.submat(begin1, 0, end1, 0);
    m2.submat(n1+1, 0, n1+n2+1, 0) = m.submat(begin2, 0, end2, 0);
}

// kron only at row direction.
template<typename T>
T kron_row(const T& A, const T& B){
    T m(A.n_rows * B.n_rows, A.n_cols);
    for (int i = 0; i < A.n_cols; i++){
        m.col(i) = arma::kron(A.col(i), B.col(i));
    }
    return m;
}

// !! Use this in a single line!
// Bound safe setting lhs.submat(...) = rhs
#define SET_SUBMAT(lhs,rbegin,rend,cbegin,cend,rhs) \
    if(rend>=rbegin && cend>=cbegin)RSUBMAT(lhs,rbegin,rend,cbegin,cend)=rhs;

// Bound safe getting m.submat(...) 
#define GET_SUBMAT(m,rbegin,rend,cbegin,cend) \
    ((rend>=rbegin && cend>=cbegin) ? RSUBMAT(m,rbegin,rend,cbegin,cend) : std::remove_reference<decltype(m)>::type() )



class ScatterCalc{
public:

    /* Inputs */

    int Nx;
    int Ny = 0;
    int Nel = 2;
    int n_ext;
    int n_stencil;  // (len(stencil)+1)/2
    int n_boundstate;   // maximum number of bound states
    int in_boundstate;  // incident bound state
    int in_el_inc;      // incident electronic state
    int n_eig_out = 20;      // eigenstate output
    double L;
    double dx = 0;
    double m;
    double E0_tot;     // initial overall energy
    enum StartPos {DOWN, UP} start_pos;

    Potential* potential;

    /* Caches */

    vec half_kinetic_mat; // -1/2m/d^2/dx^2 (nonnegative part); do not use directly
//    cx_mat H;
    arma::sp_cx_mat H;
    arma::sp_cx_mat lhs;
    cx_mat rhs;
    
    std::vector<int> offsetA, offsetB;

    /* Results */

    std::vector<vec> trans;
    std::vector<vec> refl;
    std::vector<cx_vec> trans_raw;
    std::vector<cx_vec> refl_raw;
    double trans_tot, refl_tot;
    double overall;
    double kin;
    std::vector<int> n_bsA, n_bsB;      // Final number of bound states
    std::vector<vec> kA, kB;            // wavenumber of different bound states
    std::vector<mat> psi_bsA, psi_bsB;  // bound state wavefunctions
    std::vector<cx_mat> psi;            // Ny x Nx, inside elements

    ScatterCalc(){

    }

    void setup(){

        printf("Using %d-stencil and %d state represents incident waves\n", 2*n_stencil-1, n_ext);
        printf("Total grid points: %d x %d\n", Nx, Nx);

        // geometry; The grid spans from [-L/2,L/2]
        if (dx == 0) dx = L / (Nx - 1);
        if (Ny == 0) Ny = Nx;

        printf("dx = %f; x = [%f, %f], y=[%f, %f]\n", dx, -L/2, (Nx-1)*dx - L/2, -L/2, (Ny-1)*dx - L/2);

        // Calculate Laplacian matrix
        mat half_lap;
        get_half_laplacian(half_lap);
        half_kinetic_mat = -half_lap.t() /2/m/dx/dx;

        // Calculating boundary conditions

        Nel = potential->dim();
        printf("Model dimension is %d\n", Nel);
        std::vector<mat> V_edgeA(Nel, mat(Nx, 1)), V_edgeB(Nel, mat(Nx, 1));
        std::vector<double> _V_edge_buffer(Nel);
        for (int i = 0; i < Nx; i++){
            potential->get_PE_edgeA(i*dx - L/2, _V_edge_buffer);
            for (int n = 0; n < Nel; n++){
                V_edgeA[n](i) = _V_edge_buffer[n];
            }
            potential->get_PE_edgeB(i*dx - L/2, _V_edge_buffer);
            for (int n = 0; n < Nel; n++){
                V_edgeB[n](i) = _V_edge_buffer[n];
            }
        }

        psi_bsA.resize(Nel); psi_bsB.resize(Nel);
        std::vector<vec> E_bsA(Nel), E_bsB(Nel);
        kA.resize(Nel); kB.resize(Nel);
        n_bsA.resize(Nel); n_bsB.resize(Nel);
        
        for (int n = 0; n < Nel; n++){
            compute_boundstate(V_edgeA[n], psi_bsA[n], E_bsA[n]);
            compute_boundstate(V_edgeB[n], psi_bsB[n], E_bsB[n]);
        }

        for (int n = 0; n < Nel; n++){
            printf("Bound states |%d,A> (%d):\n", n, E_bsA[n].n_elem);
            fprint_mat(E_bsA[n]);
        }
        for (int n = 0; n < Nel; n++){
            printf("Bound states |%d,B> (%d):\n", n, E_bsB[n].n_elem);
            fprint_mat(E_bsB[n]);
        }

#ifdef _DEBUG
        psi_bsA[0].print("psi0_bound=");
        psi_bsA[1].print("psi1_bound=");
#endif

        auto _compare_zero = [](const vec& lhs, const vec& rhs){return (lhs(0)<rhs(0)?true:false);};
        int n_min_EA = std::min_element(E_bsA.begin(), E_bsA.end(), _compare_zero) - E_bsA.begin();
        int n_min_EB = std::min_element(E_bsB.begin(), E_bsB.end(), _compare_zero) - E_bsB.begin();

        // reduce bound states
        if (E0_tot < std::max(E_bsA[n_min_EA](0), E_bsB[n_min_EB](0))){
            throw std::runtime_error("Energy not enough for ground state");
        }

#define FIND_FIRST_LARGER(a, v) (std::upper_bound((a).begin(), (a).end(), (v)) - (a).begin())

        for (int n = 0; n < Nel; n++){
            n_bsA[n] = FIND_FIRST_LARGER(E_bsA[n], E0_tot);
            n_bsB[n] = FIND_FIRST_LARGER(E_bsB[n], E0_tot);
            E_bsA[n] = GET_SUBMAT(E_bsA[n], 0, n_bsA[n]-1, 0, 0);
            E_bsB[n] = GET_SUBMAT(E_bsB[n], 0, n_bsB[n]-1, 0, 0);
            psi_bsA[n] = GET_SUBMAT(psi_bsA[n], 0, psi_bsA[n].n_rows-1, 0, n_bsA[n]-1);
            psi_bsB[n] = GET_SUBMAT(psi_bsB[n], 0, psi_bsB[n].n_rows-1, 0, n_bsB[n]-1);
            kA[n] = arma::sqrt(2*m*(E0_tot - E_bsA[n]));
            kB[n] = arma::sqrt(2*m*(E0_tot - E_bsB[n]));
        }
#undef FIND_FIRST_LARGER

        printf("Number of bound states taken in A:");
        for (int n = 0; n < Nel; n++){
            printf(" %d", n_bsA[n]);
        }
        printf("\nNumber of bound states taken in B:");
        for (int n = 0; n < Nel; n++){
            printf(" %d", n_bsB[n]);
        }
        printf("\nIncident from boundstate |0,%d,%s>\n", in_boundstate, start_pos==DOWN?"A":"B");
        if (in_boundstate >= (start_pos==DOWN?n_bsA:n_bsB)[start_pos==DOWN?n_min_EA:n_min_EB]){
            throw std::runtime_error("Incident boundstate too large");
        }
        printf("Outcoming wavenumbers in A:\n");
        for (int n = 0; n < Nel; n++){
            fprint_mat(kA[n], stdout, ("state " + STR(n)).c_str());
        }
        printf("Outcoming wavenumbers in B:\n");
        for (int n = 0; n < Nel; n++){
            fprint_mat(kB[n], stdout, ("state " + STR(n)).c_str());
        }

        // The model is divided by A-side (down) and B-side (up)
        // Belows are matrices of Hamiltonian components
        // the j'th column/element represent j'th bound state

        std::vector<arma::cx_rowvec> edgeA_H_inA(Nel), edgeA_H_outA(Nel), edgeB_H_inB(Nel), edgeB_H_outB(Nel);
        std::vector<cx_mat> X_T_edgeA(Nel), X_T_edgeB(Nel);
        std::vector<cx_mat> X_T_inA(Nel), X_T_outA(Nel), X_T_inB(Nel), X_T_outB(Nel);
        std::vector<cx_mat> edgeA(Nel), edgeB(Nel);

        for (int n = 0; n < Nel; n++){
            edgeA[n] = to_complex(psi_bsA[n]);
            edgeB[n] = to_complex(psi_bsB[n]);
        }

        std::vector<cx_mat> inA_trans(Nel), outA_trans(Nel), inB_trans(Nel), outB_trans(Nel);
        std::vector<cx_mat> inA(Nel), outA(Nel), inB(Nel), outB(Nel);

        // To modify scatter coordinates, modify here.

        if (potential->is_angular()){

        // perpendicular components of kets; j'th column is j'th bound state
        for (int n = 0; n < Nel; n++){
            inA_trans[n] = arma::exp(I*dx*arange(-n_ext+1, 0) * kA[n].t());
            outA_trans[n] = arma::exp(-I*dx*arange(-n_ext+1, 0) * kA[n].t());
            inB_trans[n] = arma::exp(I*dx*arange(-n_ext+1, 0) * kB[n].t());
            outB_trans[n] = arma::exp(-I*dx*arange(-n_ext+1, 0) * kB[n].t());
        }

        for (int n = 0; n < Nel; n++){
            inA[n] = kron_row(edgeA[n], inA_trans[n]);
            outA[n] = kron_row(edgeA[n], outA_trans[n]);
            inB[n] = kron_row(inB_trans[n], edgeB[n]);  // since B goes x direction
            outB[n] = kron_row(outB_trans[n], edgeB[n]);
        }

        for (int n = 0; n < Nel; n++){
            edgeA_H_inA[n] = prod_Kmat_1d({-n_ext+1, 0}, inA_trans[n]) + E_bsA[n].t();
            edgeA_H_outA[n] = prod_Kmat_1d({-n_ext+1, 0}, outA_trans[n]) + E_bsA[n].t();
            edgeB_H_inB[n] = prod_Kmat_1d({-n_ext+1, 0}, inB_trans[n]) + E_bsB[n].t();
            edgeB_H_outB[n] = prod_Kmat_1d({-n_ext+1, 0}, outB_trans[n]) + E_bsB[n].t();
        }

        for (int n = 0; n < Nel; n++){
            times_Kmat({1,Nx,0,0}, {1,Nx,1,Ny}, edgeA[n], X_T_edgeA[n]);
            times_Kmat({0,0,1,Ny}, {1,Nx,1,Ny}, edgeB[n], X_T_edgeB[n]);
        }

        for (int n = 0; n < Nel; n++){
            times_Kmat({1,Nx,-n_ext+1,0}, {1,Nx,1,Ny}, inA[n], X_T_inA[n]);
            times_Kmat({1,Nx,-n_ext+1,0}, {1,Nx,1,Ny}, outA[n], X_T_outA[n]);
            times_Kmat({-n_ext+1,0,1,Ny}, {1,Nx,1,Ny}, inB[n], X_T_inB[n]);
            times_Kmat({-n_ext+1,0,1,Ny}, {1,Nx,1,Ny}, outB[n], X_T_outB[n]);
        }

        }
        else{   // straight

        printf("The potential is straight\n");

        for (int n = 0; n < Nel; n++){
            inA_trans[n] = arma::exp(I*dx*arange(-n_ext+1, 0) * kA[n].t());
            outA_trans[n] = arma::exp(-I*dx*arange(-n_ext+1, 0) * kA[n].t());
            inB_trans[n] = arma::exp(-I*dx*arange(0, n_ext-1) * kB[n].t());
            outB_trans[n] = arma::exp(I*dx*arange(0, n_ext-1) * kB[n].t());
        }

        for (int n = 0; n < Nel; n++){
            inA[n] = kron_row(edgeA[n], inA_trans[n]);
            outA[n] = kron_row(edgeA[n], outA_trans[n]);
            inB[n] = kron_row(edgeB[n], inB_trans[n]);
            outB[n] = kron_row(edgeB[n], outB_trans[n]);
        }

        for (int n = 0; n < Nel; n++){
            edgeA_H_inA[n] = prod_Kmat_1d({-n_ext+1, 0}, inA_trans[n]) + E_bsA[n].t();
            edgeA_H_outA[n] = prod_Kmat_1d({-n_ext+1, 0}, outA_trans[n]) + E_bsA[n].t();
            edgeB_H_inB[n] = prod_Kmat_1d({0, n_ext-1}, inB_trans[n]) + E_bsB[n].t();
            edgeB_H_outB[n] = prod_Kmat_1d({0, n_ext-1}, outB_trans[n]) + E_bsB[n].t();
        }

        for (int n = 0; n < Nel; n++){
            times_Kmat({1,Nx,0,0}, {1,Nx,1,Ny}, edgeA[n], X_T_edgeA[n]);
            times_Kmat({1,Nx,Ny+1,Ny+1}, {1,Nx,1,Ny}, edgeB[n], X_T_edgeB[n]);
        }

        for (int n = 0; n < Nel; n++){
            times_Kmat({1,Nx,-n_ext+1,0}, {1,Nx,1,Ny}, inA[n], X_T_inA[n]);
            times_Kmat({1,Nx,-n_ext+1,0}, {1,Nx,1,Ny}, outA[n], X_T_outA[n]);
            times_Kmat({1,Nx,Ny+1,Ny+n_ext}, {1,Nx,1,Ny}, inB[n], X_T_inB[n]);
            times_Kmat({1,Nx,Ny+1,Ny+n_ext}, {1,Nx,1,Ny}, outB[n], X_T_outB[n]);
        }

        }

        int N = Nx*Ny;
        int n_bsA_tot = std::accumulate(n_bsA.begin(), n_bsA.end(), 0);
        int n_bsB_tot = std::accumulate(n_bsB.begin(), n_bsB.end(), 0);
        offsetA.resize(Nel); offsetB.resize(Nel);
        offsetA[0] = Nel*N;
        offsetB[0] = Nel*N + n_bsA_tot;
        for (int n = 1; n < Nel; n++){
            offsetA[n] = offsetA[n-1] + n_bsA[n-1];
            offsetB[n] = offsetB[n-1] + n_bsB[n-1];
        }

        // The Hamiltonian is 2N + 2(bs[0]+bs[1]): state0, state1, Aout, Bout
        H = arma::sp_cx_mat(N*Nel + n_bsA_tot + n_bsB_tot, N*Nel + n_bsA_tot + n_bsB_tot);
        rhs = arma::zeros<cx_mat>(N*Nel + n_bsA_tot + n_bsB_tot, 1);

        for (int n = 0; n < Nel; n++){
            SET_SUBMAT(H, n*N, (n+1)*N-1, offsetA[n], offsetA[n]+n_bsA[n]-1, X_T_outA[n])
            SET_SUBMAT(H, offsetA[n], offsetA[n]+n_bsA[n]-1, n*N, (n+1)*N-1, X_T_edgeA[n].t())
            SET_SUBMAT(H, offsetA[n], offsetA[n]+n_bsA[n]-1, offsetA[n], offsetA[n]+n_bsA[n]-1, arma::diagmat(edgeA_H_outA[n]))

            SET_SUBMAT(H, n*N, (n+1)*N-1, offsetB[n], offsetB[n]+n_bsB[n]-1, X_T_outB[n])
            SET_SUBMAT(H, offsetB[n], offsetB[n]+n_bsB[n]-1, n*N, (n+1)*N-1, X_T_edgeB[n].t())
            SET_SUBMAT(H, offsetB[n], offsetB[n]+n_bsB[n]-1, offsetB[n], offsetB[n]+n_bsB[n]-1, arma::diagmat(edgeB_H_outB[n]))
        }

        if (start_pos == DOWN){

            // incoming from |0A, 0>

            if (in_el_inc < 0) in_el_inc = n_min_EA;
            kin = kA[in_el_inc](in_boundstate);
            rhs(offsetA[in_el_inc] + in_boundstate) = E0_tot - edgeA_H_inA[in_el_inc](in_boundstate);
            SET_SUBMAT(rhs, in_el_inc*N, (in_el_inc+1)*N-1, 0, 0, -X_T_inA[in_el_inc].col(in_boundstate))
        }
        else if (start_pos == UP){

            if (in_el_inc < 0) in_el_inc = n_min_EB;
            kin = kB[in_el_inc](in_boundstate);
            rhs(offsetB[in_el_inc] + in_boundstate) = E0_tot - edgeB_H_inB[in_el_inc](in_boundstate);
            SET_SUBMAT(rhs, in_el_inc*N, (in_el_inc+1)*N-1, 0, 0, -X_T_inB[in_el_inc].col(in_boundstate))
        }

        printf("Incident from |%d,%s>\n", in_el_inc, start_pos==DOWN?"A":"B");

        mat K;
        get_Kmat({1,Nx,1,Ny}, {1,Nx,1,Ny}, K);
        for (int n = 0; n < Nel; n++){
            H.submat(n*N, n*N, (n+1)*N-1, (n+1)*N-1) = to_complex(K);
        }

#define LOC(x_, y_, s_) ((x_*Ny + y_) + s_ * (N))

        // Potential energy
#pragma omp parallel for
        for (int i = 0; i < Nx; i++){
            std::vector<std::complex<double> > _V_buffer(Nel*Nel, 0.0);
            for (int j = 0; j < Ny; j++){
                potential->get_PE(i*dx - L/2, j*dx - L/2, _V_buffer);

                for (int n1 = 0; n1 < Nel; n1++){
                    for (int n2 = 0; n2 < Nel; n2++){
                        H(LOC(i, j, n1), LOC(i, j, n2)) += _V_buffer[n1*Nel + n2];
                    }
                }
            }
        }
#undef LOC

#ifdef _DEBUG
        print_cmat(H);
        print_cmat(rhs);
#endif

        lhs = arma::sp_cx_mat(H);
        lhs.diag() -= E0_tot;
    }

    // Calculate and return bound states with given potential
    void compute_boundstate(const mat& V_edge, mat& psi_b, vec& E_b){

        mat H_edge(V_edge.size(), V_edge.size(), arma::fill::zeros);
        get_Kmat_1d({1,(int)V_edge.size()}, {1,(int)V_edge.size()}, H_edge);
        H_edge += arma::diagmat(V_edge);

        arma::eig_sym(E_b, psi_b, H_edge);
        E_b = E_b.subvec(0, std::min(n_boundstate, (int)E_b.n_elem)-1);
        psi_b = psi_b.submat(0, 0, psi_b.n_rows-1, std::min(n_boundstate, (int)psi_b.n_cols)-1);
    }

    void compute_eigenstates(cx_mat& m_psi, vec& E){
        int N = Nx*Ny;
        arma::cx_mat H_sys = arma::conv_to<arma::cx_mat>::from(arma::sp_cx_mat(H.submat(0,0,N*Nel-1,N*Nel-1)));

        printf("Calculating eigenstates...\n");

        arma::eig_sym(E, m_psi, H_sys);
        for (int n = 0; n < n_eig_out; n++){
            for (int j = 0; j < Nel; j++){
                psi.push_back(arma::reshape(m_psi.submat(j*N, n, (j+1)*N-1, n), Ny, Nx));
            }
        }
    }

    void compute_exponential(cx_mat& exp_m){
        
        cx_mat U;
        vec E;
        compute_eigenstates(U, E);
        exp_m = U * arma::diagmat(arma::exp(-I*E)) * U.t();
    }

    void print_potential(FILE* fp){
        std::vector<cx_mat> V(Nel*Nel);
        potential->generate(arange(0,Nx-1)*dx - L/2, arange(0,Ny-1)*dx - L/2, V);
        for (int n1 = 0; n1 < Nel; n1++){
            for (int n2 = 0; n2 < Nel; n2++){
                fprint_cmat(V[n1*Nel+n2], fp, "", false);
            }
        }
    }

    void get_half_laplacian(mat& out)const{
        switch(n_stencil){
            case 2:out = {-2, 1};break;
            case 3:out = {-5/2.0, 4/3.0, -1/12.0}; break;
            case 4:out = {-49/18.0, 3/2.0, -3/20.0, 1/90.0}; break;
            case 5:out = {-205/72.0, 8/5.0, -1/5.0, 8/315.0, -1/560.0}; break;
            default: throw std::runtime_error("Bad stencil number");break;
        }
    }

    // Calculate <xr,ys|T|psi> = \sum_{pq}{<xr,ys|T|xp,yq><psi_{pq}>}
    // !! Convention of size: xbegin, xend, ybegin, yend
    void times_Kmat(std::array<int,4> size_in, std::array<int,4> size_out, const cx_mat& mat_in, cx_mat& mat_out){
        if (mat_in.n_elem == 0) return;
        arma::mat K;
        get_Kmat(size_out, size_in, K);
        mat_out = K * mat_in;
    }

    // Calculate <x0|T|psi> in 1D
    cx_mat prod_Kmat_1d(std::array<int,2> size_in, const cx_mat& mat_in){
        arma::mat K;
        get_Kmat_1d({0, 0}, size_in, K);
        return K * mat_in;
    }

    // Calculate <xp,yp|T|xq,yq>
    // !! Convention of size: xbegin, xend (include), ybegin, yend (include)
    void get_Kmat(std::array<int, 4> size_bra, std::array<int, 4> size_ket, mat& out){
        
        int nx_bra = size_bra[1]-size_bra[0]+1,
            ny_bra = size_bra[3]-size_bra[2]+1,
            nx_ket = size_ket[1]-size_ket[0]+1,
            ny_ket = size_ket[3]-size_ket[2]+1;

        out = mat(nx_bra*ny_bra, nx_ket*ny_ket, arma::fill::zeros);

        for (int px = size_bra[0]; px <= size_bra[1]; px++){
            for (int py = size_bra[2]; py <= size_bra[3]; py++){
                for (int qx = std::max(px-n_stencil+1, size_ket[0]); qx <= std::min(px+n_stencil-1, size_ket[1]); qx++){
                    for (int qy = std::max(py-n_stencil+1, size_ket[2]); qy <= std::min(py+n_stencil-1, size_ket[3]); qy++){                        
                        out((px-size_bra[0])*ny_bra + py-size_bra[2], (qx-size_ket[0])*ny_ket + qy-size_ket[2]) = 
                        (py==qy? half_kinetic_mat(abs(px-qx)):0.0) + (px==qx?half_kinetic_mat(abs(py-qy)):0.0);
                    }
                }
            }
        }
    }

    // Calculate <xp|T|xq>
    void get_Kmat_1d(std::array<int, 2> size_bra, std::array<int, 2> size_ket, mat& out){

        out = mat(size_bra[1]-size_bra[0]+1, size_ket[1]-size_ket[0]+1, arma::fill::zeros);
        for (int p = size_bra[0]; p <= size_bra[1]; p++){
            for (int q = std::max(p-n_stencil+1, size_ket[0]); q <= std::min(p+n_stencil-1, size_ket[1]); q++){
                out(p-size_bra[0], q-size_ket[0]) = half_kinetic_mat(abs(p-q));
            }
        }
    }

    void run(){
        printf("Start running\n");

        cx_vec ret = arma::spsolve(lhs, rhs);
        
        int N = Nx*Ny;

        refl.resize(Nel);
        trans.resize(Nel);
        refl_raw.resize(Nel);
        trans_raw.resize(Nel);
        printf("kin=%.10f\n", kin);

        if (start_pos == DOWN){
            for (int n = 0; n < Nel; n++){
                refl[n] = kA[n] % arma::pow(arma::abs(GET_SUBMAT(ret, offsetA[n], offsetA[n]+n_bsA[n]-1, 0, 0)), 2) / kin;
                trans[n] = kB[n] % arma::pow(arma::abs(GET_SUBMAT(ret, offsetB[n], offsetB[n]+n_bsB[n]-1, 0, 0)), 2) / kin;
                refl_raw[n] = GET_SUBMAT(ret, offsetA[n], offsetA[n]+n_bsA[n]-1, 0, 0);
                trans_raw[n] = GET_SUBMAT(ret, offsetB[n], offsetB[n]+n_bsB[n]-1, 0, 0);
            }
        }
        else{
            for (int n = 0; n < Nel; n++){
                refl[n] = kB[n] % arma::pow(arma::abs(GET_SUBMAT(ret, offsetB[n], offsetB[n]+n_bsB[n]-1, 0, 0)), 2) / kin;
                trans[n] = kA[n] % arma::pow(arma::abs(GET_SUBMAT(ret, offsetA[n], offsetA[n]+n_bsA[n]-1, 0, 0)), 2) / kin;
                refl_raw[n] = GET_SUBMAT(ret, offsetB[n], offsetB[n]+n_bsB[n]-1, 0, 0);
                trans_raw[n] = GET_SUBMAT(ret, offsetA[n], offsetA[n]+n_bsA[n]-1, 0, 0);
            }
        }

        for (int n = 0; n < Nel; n++){
            psi.push_back(arma::reshape(GET_SUBMAT(ret, n*N, (n+1)*N-1, 0, 0), Ny, Nx));
        }

        refl_tot = 0.0; trans_tot = 0.0;
        for (int n = 0; n < Nel; n++){
            refl_tot += arma::accu(refl[n]);
            trans_tot += arma::accu(trans[n]);
        }
        
        overall = refl_tot + trans_tot;

        printf("Reflection: %.10f\nTransmission: %.10f\nSum: %.10f\n", refl_tot, trans_tot, overall);
    }

};


Potential* get_potential(const char* name){
#define preprocfac 
preprocfac name;
}

int main(int argc, char** argv){

    ScatterCalc calc;

    cxxopts::Options desc("Scattering", "Scattering calculation");
    desc.add_options()
        ("help", "produce help message")
        ("N", "N", cxxopts::value<int>(calc.Nx))
        ("Ny", "Ny", cxxopts::value<int>(calc.Ny))
        ("dx", "dx", cxxopts::value<double>(calc.dx))
        ("n_ext", "n_ext", cxxopts::value<int>(calc.n_ext)->default_value("3"))
        ("n_stencil", "n_stencil", cxxopts::value<int>(calc.n_stencil)->default_value("3"))
        ("n_bs", "n_boundstate", cxxopts::value<int>(calc.n_boundstate)->default_value("1"))
        ("in_bs", "Incident boundstate", cxxopts::value<int>(calc.in_boundstate)->default_value("0"))
        ("in_el", "Incident electronic state", cxxopts::value<int>(calc.in_el_inc)->default_value("-1"))
        ("L", "L", cxxopts::value<double>(calc.L))
        ("m", "m", cxxopts::value<double>(calc.m)->default_value("1.0"))
        ("E", "Overall energy", cxxopts::value<double>(calc.E0_tot))
        ("startpos", "up/down", cxxopts::value<std::string>())
        ("potential", "angular", cxxopts::value<std::string>())
        ("param", "PE params", cxxopts::value<std::vector<double> >())
        ("eig", "Compute eigenstates", cxxopts::value<bool>())
        ("pe", "Output potential energy", cxxopts::value<bool>())
        ("output", "output prefix", cxxopts::value<std::string>()->default_value("out"))
        ;

    auto vm = desc.parse(argc, argv);
    if (vm.count("help")) {
        std::cout << desc.help() << std::endl;
        return 0;
    }

    std::string startpos = vm["startpos"].as<std::string>();
    if (startpos == "up"){
        calc.start_pos = ScatterCalc::UP;
    }
    else if (startpos == "down"){
        calc.start_pos = ScatterCalc::DOWN;
    }
    std::string outputprefix = vm["output"].as<std::string>();
    std::string potential = vm["potential"].as<std::string>();
    calc.potential = get_potential(potential.c_str()); 

    if (vm.count("param")){
        calc.potential->init(vm["param"].as<std::vector<double> >());
    }
    calc.setup();
    if (vm.count("pe")){
    FILE* f_pe = fopen(("pe-" + outputprefix + ".txt").c_str(), "w");
    calc.print_potential(f_pe);
    fclose(f_pe);
    }
    if (vm.count("eig")){
        cx_mat psi;
        vec E_eig;
        calc.compute_eigenstates(psi, E_eig);
        E_eig.subvec(0, calc.n_eig_out-1).print();
    }
    else{
    calc.run();

    if (outputprefix != "silent"){
    FILE* f_prob = fopen(("prob-" + outputprefix + ".txt").c_str(), "w");
    fprintf(f_prob, "refl=%.10g\ntrans=%.10g\noverall=%.10g\n",
     calc.refl_tot, calc.trans_tot, calc.overall);
        for (size_t n = 0; n < calc.refl.size(); n++){
            fprint_mat(calc.refl[n], f_prob, ("refl"+STR(n)).c_str());
        }
        for (size_t n = 0; n < calc.trans.size(); n++){
            fprint_mat(calc.trans[n], f_prob, ("trans"+STR(n)).c_str());
        }
        for (size_t n = 0; n < calc.refl_raw.size(); n++){
            fprint_cmat(calc.refl_raw[n], f_prob, ("reflraw"+STR(n)).c_str());
        }
        for (size_t n = 0; n < calc.trans_raw.size(); n++){
            fprint_cmat(calc.trans_raw[n], f_prob, ("transraw"+STR(n)).c_str());
        }
    fclose(f_prob);
    }
    }

    if (outputprefix != "silent"){
    FILE* f_matrix = fopen(("state-" + outputprefix + ".txt").c_str(), "w");
    for(size_t n = 0; n < calc.psi.size(); n++){
        fprint_cmat(calc.psi[n], f_matrix);
    }
    fclose(f_matrix);

    FILE* f_bs = fopen(("bound-" + outputprefix + ".txt").c_str(), "w");
    for(size_t n = 0; n < calc.psi_bsA.size(); n++){
        fprint_mat(calc.psi_bsA[n], f_bs, ("boundstates at |" + STR(n) + ",A>").c_str());
    }
    for(size_t n = 0; n < calc.psi_bsB.size(); n++){
        fprint_mat(calc.psi_bsB[n], f_bs, ("boundstates at |" + STR(n) + ",B>").c_str());
    }
    fclose(f_bs);
    }
    return 0;
}

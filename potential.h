
#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <armadillo>

class Potential{
public:
    virtual void init(const std::vector<double>& param) = 0;
    virtual void get_PE(double x, double y, std::vector<std::complex<double> >& ) = 0;  // v(i*n+j) =>V_{ij}
    virtual void get_PE_edgeA(double x, std::vector<double>&) = 0;   // v(i) => V_{ii}
    virtual void get_PE_edgeB(double x, std::vector<double>& V_edge) {
        if (symmetry() == 1){
            get_PE_edgeA(x, V_edge);
        }
        else if (symmetry() == -1){
            get_PE_edgeA(x, V_edge);
            std::reverse(V_edge.begin(), V_edge.end());
        }
        else{
            throw std::runtime_error("Undefined symmetry");
        }
    }
    virtual int symmetry()const = 0;
    virtual int dim()const{
        return 2;
    }
    virtual bool is_angular()const{
        return true;
    }

    // !Matrix coord is [y,x]
    void generate(const arma::mat& xrange, const arma::mat& yrange, std::vector<arma::cx_mat>& V){
        std::vector<std::complex<double> > _V_local(V.size(), 0.0);
        for (size_t n = 0; n < V.size(); n++){
            V[n] = arma::cx_mat(xrange.n_elem, yrange.n_elem, arma::fill::zeros);
        }

        for (int i = 0; i < xrange.n_elem; i++){
            for (int j = 0; j < yrange.n_elem; j++){
                get_PE(xrange(i), yrange(j), _V_local);
                for (size_t n = 0; n < V.size(); n++){
                    V[n](i, j) = _V_local[n];
                }
            }
        }
    }

};

#endif
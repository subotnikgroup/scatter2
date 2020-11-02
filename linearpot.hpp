#include "potential.h"

class BarrierPotential : public Potential {
public:
    double width = 0.1;
    double height = 1.0;

    void init(const std::vector<double>& param){
        width = param[0];
        height = param[1];
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        V[0] = (y>=0 && y < width) ? height : 0;
        V[3] = 0.1;
        V[2] = V[1] = 0.0; 
    }

    void get_PE_edgeA(double x, std::vector<double>& V ){
        V[0] = 0;
        V[1] = 0.1;
    }

    int symmetry()const{
        return 1;
    }

    bool is_angular()const{
        return false;
    }
    static const char* get_name(){
        return "barrier";
    }

};

class Tully1Potential : public Potential {
public:
    double A = 0.01;
    double B = 1.6;
    double C = 0.005;
    double D = 1.0;

    void init(const std::vector<double>& param){
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        double E1 = y > 0 ? A * (1 - exp(-B*y)) : -A*(1 - exp(B*y));
        double V12 = C*exp(-D*y*y);
        
            V[0] = E1 + A;
            V[3] = -E1 + A;
            V[1] = V[2] = V12;
        
    }

    void get_PE_edgeA(double x, std::vector<double>& V ){
        V[0] = 0;
        V[1] = 2*A;
    }

    int symmetry()const{
        return -1;
    }

    bool is_angular()const{
        return false;
    }

    static const char* get_name(){
        return "tully1";
    }

};


class LinearComplexVectorWPotential : public Potential {
public:

    double L = 1.0;
    double A = 1.0;
    double eps = 5.0;
    double C = 1.0;
    double W = 1.0;
    double ltheta = arma::datum::pi/2;
    const double forbidden = 1000.0;

    void init(const std::vector<double>& param){
        L = param[0];
        A = param[1];
        eps = param[2];
        C = param[3];
        W = param[4];
        ltheta = param[5];
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        
        if (x < -L/2 || x > L/2){
            V[0] = forbidden;
            V[3] = forbidden;
            V[1] = V[2] = 0.0;
            return;
        }

        V[0] = tanh(y*eps) * A + A;
        V[3] = tanh(-y*eps) * A + A;
        V[1] = C * exp(-eps*eps*y*y) * std::exp(std::complex<double>(0, W*(cos(ltheta)*x + sin(ltheta)*y)));
        V[2] = std::conj(V[1]);

    }

    void get_PE_edgeA(double x, std::vector<double>& V){

        if (x < -L/2 || x > L/2){
            V[0] = forbidden;
            V[1] = forbidden;
            return;
        }
        V[0] = 0;
        V[1] = 2*A;
    }

    int symmetry()const{
        return -1;
    }
    bool is_angular()const{
        return false;
    }
    static const char* get_name(){
        return "linel";
    }

};


class LinearComplexRadialPotential : public Potential {
public:

    double L = 1.0;
    double A = 1.0;
    double eps = 5.0;
    double C = 1.0;
    double W = 1.0;
    double R = 1.0;
    const double forbidden = 1000.0;

    void init(const std::vector<double>& param){
        L = param[0];
        A = param[1];
        eps = param[2];
        C = param[3];
        W = param[4];
        R = param[5];
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        
        if (x < -L/2 || x > L/2){
            V[0] = forbidden;
            V[3] = forbidden;
            V[1] = V[2] = 0.0;
            return;
        }

        double r = std::sqrt((x+R)*(x+R) + y*y);

        V[0] = tanh(y*eps) * A + A;
        V[3] = tanh(-y*eps) * A + A;
        V[1] = C * exp(-eps*eps*y*y) * std::exp(std::complex<double>(0, W*r));
        V[2] = std::conj(V[1]);

    }

    void get_PE_edgeA(double x, std::vector<double>& V){

        if (x < -L/2 || x > L/2){
            V[0] = forbidden;
            V[1] = forbidden;
            return;
        }
        V[0] = 0;
        V[1] = 2*A;
    }

    int symmetry()const{
        return -1;
    }
    bool is_angular()const{
        return false;
    }
    static const char* get_name(){
        return "line";
    }

};

class LinearComplexPotentialEx : public Potential {
public:

    double L = 1.0;
    double A = 1.0;
    double eps = 5.0;
    double C = 1.0;
    double W = 1.0;
    int mode = 0;
    double k;
    double alpha = 1.0;

    const double forbidden = 1000.0;

    void init(const std::vector<double>& param){
        L = param[0];
        A = param[1];
        eps = param[2];
        C = param[3];
        W = param[4];
        mode = (int)param[5];
        k = param[6];

        if (mode == 1 && param.size() > 7){
            alpha = param[7];
        }
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        
        if (x < -L/2 || x > L/2){
            V[0] = forbidden;
            V[3] = forbidden;
            V[1] = V[2] = 0.0;
            return;
        }

        if (mode == 0){
            V[0] = tanh(y*eps) * A + A;
            V[3] = tanh(-y*eps) * A + A;
            V[1] = C * exp(-eps*eps*y*y) * std::exp(std::complex<double>(0, W*x));        
        }
        else if (mode == 1){    // rectangle trapezoid
            double shift = alpha / (eps + k*x) - alpha / eps;

            V[0] = tanh((y-shift)*(eps+k*x))*A + A;
            V[3] = -tanh((y-shift)*(eps+k*x))*A + A;
            V[1] = C * exp(-(eps+k*x)*(eps+k*x)*(y-shift)*(y-shift)) * std::exp(std::complex<double>(0, W*x)); 
        }
        else if (mode == 2){
            V[0] = tanh((y-k*x)*eps)*A + A;
            V[3] = -tanh((y-k*x)*eps)*A + A;
            V[1] = C * exp(-eps*eps*(y-k*x)*(y-k*x)) * std::exp(std::complex<double>(0, W*x)); 
        }
        else if (mode == 3){
            V[0] = tanh(y*(eps+k*x))*A + A;
            V[3] = -tanh(y*(eps+k*x))*A + A;
            V[1] = C * exp(-(eps+k*x)*(eps+k*x)*y*y) * std::exp(std::complex<double>(0, W*x)); 
        }

        V[2] = std::conj(V[1]);

    }

    void get_PE_edgeA(double x, std::vector<double>& V){

        if (x < -L/2 || x > L/2){
            V[0] = forbidden;
            V[1] = forbidden;
            return;
        }
        V[0] = 0;
        V[1] = 2*A;
    }

    int symmetry()const{
        return -1;
    }
    bool is_angular()const{
        return false;
    }
    static const char* get_name(){
        return "lineex";
    }

};


class LinearComplexRadial4StatePotential : public Potential {
public:

    double L = 1.0;
    double A = 1.0;
    double eps = 5.0;
    double C = 1.0;
    double W = 1.0;
    double R = 1.0;
    double T = 1.0;
    const double forbidden = 1000.0;

    void init(const std::vector<double>& param){
        L = param[0];
        A = param[1];
        eps = param[2];
        C = param[3];
        W = param[4];
        R = param[5];
        T = param[6];
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        
#define LOC(a,b) (a*4+b)
        if (x < -L/2 || x > L/2){
            std::fill(V.begin(), V.end(), 0.0);
            V[LOC(0,0)] = V[LOC(1,1)] = V[LOC(2,2)] = V[LOC(3,3)] = forbidden;
            return;
        }

        double r = std::sqrt((x+R)*(x+R) + y*y);
        std::complex<double> s = C * exp(-eps*eps*y*y) * std::exp(std::complex<double>(0, W*r));

        V[LOC(1,1)] = V[LOC(0,0)] = tanh(y*eps) * A + A;
        V[LOC(3,3)] = V[LOC(2,2)] = tanh(-y*eps) * A + A;

        V[LOC(2,0)] = V[LOC(1,3)] = V[LOC(0,2)] = V[LOC(3,1)] = T * exp(-eps*eps*y*y);
        V[LOC(0,3)] = s;
        V[LOC(3,0)] = std::conj(s);
        V[LOC(2,1)] = -s;
        V[LOC(1,2)] = -std::conj(s);

        V[LOC(0,1)] = V[LOC(1,0)] = V[LOC(2,3)] = V[LOC(3,2)] = 0.0;
#undef LOC

    }

    void get_PE_edgeA(double x, std::vector<double>& V){

        if (x < -L/2 || x > L/2){
            V[0] = forbidden;
            V[1] = forbidden;
            V[2] = forbidden;
            V[3] = forbidden;
            return;
        }
        V[0] = V[1] = 0;
        V[2] = V[3] = 2*A;
    }

    int symmetry()const{
        return -1;
    }
    bool is_angular()const{
        return false;
    }
    int dim()const{
        return 4;
    }
    static const char* get_name(){
        return "line4";
    }

};


class ThreeStateStraightPotential : public Potential {
public:
    double A = 0.02;
    double eps = 5;
    double omega_a = 0.005;
    double omega_b = 0.005;
    double M = 1000;
    double C = 0.005;

    double L = 2;
    double r0 = 4;
    double lift = 0;

    std::complex<double> V12 = {0.5, 0.5};
    std::complex<double> V13 = {0.5, 0.1};
    std::complex<double> V23 = {0.1, 0.5};

    const double forbidden = 1000;

    void init(const std::vector<double>& params){
        A = params[0]; 
        eps = params[1];
        omega_a = params[2];
        omega_b = params[3];
        M = params[4];
        C = params[5];
        L = params[6];
        r0 = params[7];

        if (params.size() > 8)
            lift = params[8];

        if (params.size() > 9){
            V12.real(params[9]);
            V12.imag(params[10]);
            V13.real(params[11]);
            V13.imag(params[12]);
            V23.real(params[13]);
            V23.imag(params[14]);
        }
    }

   void get_PE(double x, double y, std::vector<std::complex<double> >& V){
#define LOC(a,b) (a*3+b)
        if (x < -L/2 || x > L/2){
            std::fill(V.begin(), V.end(), 0.0);
            V[LOC(0,0)] = V[LOC(1,1)] = V[LOC(2,2)] = forbidden;
            return;
        }

        double Et_a = 0.5 * M * omega_a*omega_a * x*x;
        double Et_b = 0.5 * M * omega_b*omega_b * x*x;
        double v = C * exp(-std::pow(y*eps/r0, 2));

        V[LOC(0,0)] = Et_a + A * tanh(y*eps/r0) + A;
        V[LOC(1,1)] = Et_a - A * tanh(y*eps/r0) + A;
        V[LOC(2,2)] = Et_b + 2*A + lift;

        V[LOC(0,1)] = v*V12;
        V[LOC(0,2)] = v*V13;
        V[LOC(1,2)] = v*V23;
        V[LOC(1,0)] = v*std::conj(V12);
        V[LOC(2,0)] = v*std::conj(V13);
        V[LOC(2,1)] = v*std::conj(V23);

#undef LOC
   }

    void get_PE_edgeA(double x, std::vector<double>& V){
        if (x < -L/2 || x > L/2){
            V[0] = V[1] = V[2] = forbidden;
        }
        else{
            double Et_a = 0.5 * M * omega_a*omega_a * x*x;
            double Et_b = 0.5 * M * omega_b*omega_b * x*x;
            V[0] = Et_a;
            V[1] = 2*A + Et_a;
            V[2] = 2*A + lift + Et_b;
        }
    }
    void get_PE_edgeB(double x, std::vector<double>& V){
        if (x < -L/2 || x > L/2){
            V[0] = V[1] = V[2] = forbidden;
        }
        else{
            double Et_a = 0.5 * M * omega_a*omega_a * x*x;
            double Et_b = 0.5 * M * omega_b*omega_b * x*x;
            V[0] = 2*A + Et_a;
            V[1] = Et_a;
            V[2] = 2*A + lift + Et_b;
        }
    }
    bool is_angular()const {
        return false;
    }
    int dim()const{
        return 3;
    }
    int symmetry()const{
        return 0;
    }
    static const char* get_name(){
        return "3statel";
    }

};


class ThreeStateStraight2Potential : public Potential {
public:
    double A = 0.02;
    double eps = 5;
    double omega_a = 0.005;
    double omega_b = 0.005;
    double M = 1000;
    double C = 0.005;

    double L = 2;
    double r0 = 4;
    double lift = 0;

    std::complex<double> V12 = {0.5, 0.1};
    std::complex<double> V13 = {0.1, 0.5};
    std::complex<double> V23 = {0.1, 0.1};

    const double forbidden = 1000;

    void init(const std::vector<double>& params){
        A = params[0]; 
        eps = params[1];
        omega_a = params[2];
        omega_b = params[3];
        M = params[4];
        C = params[5];
        L = params[6];
        r0 = params[7];

        if (params.size() > 8)
            lift = params[8];

        if (params.size() > 9){
            V12.real(params[9]);
            V12.imag(params[10]);
            V13.real(params[11]);
            V13.imag(params[12]);
            V23.real(params[13]);
            V23.imag(params[14]);
        }
    }

   void get_PE(double x, double y, std::vector<std::complex<double> >& V){
#define LOC(a,b) (a*3+b)
        if (x < -L/2 || x > L/2){
            std::fill(V.begin(), V.end(), 0.0);
            V[LOC(0,0)] = V[LOC(1,1)] = V[LOC(2,2)] = forbidden;
            return;
        }

        double Et_a = 0.5 * M * omega_a*omega_a * x*x;
        double Et_b = 0.5 * M * omega_b*omega_b * x*x;
        double v = C * exp(-std::pow(y*eps/r0, 2));

        V[LOC(0,0)] = Et_a + A * tanh(y*eps/r0) + A;
        V[LOC(1,1)] = Et_a - A * tanh(y*eps/r0) + A;
        V[LOC(2,2)] = Et_b - A * tanh(y*eps/r0) + A + lift;

        V[LOC(0,1)] = v*V12;
        V[LOC(0,2)] = v*V13;
        V[LOC(1,2)] = v*V23;
        V[LOC(1,0)] = v*std::conj(V12);
        V[LOC(2,0)] = v*std::conj(V13);
        V[LOC(2,1)] = v*std::conj(V23);

#undef LOC
   }

    void get_PE_edgeA(double x, std::vector<double>& V){
        if (x < -L/2 || x > L/2){
            V[0] = V[1] = V[2] = forbidden;
        }
        else{
            double Et_a = 0.5 * M * omega_a*omega_a * x*x;
            double Et_b = 0.5 * M * omega_b*omega_b * x*x;
            V[0] = Et_a;
            V[1] = 2*A + Et_a;
            V[2] = 2*A + lift + Et_b;
        }
    }
    void get_PE_edgeB(double x, std::vector<double>& V){
        if (x < -L/2 || x > L/2){
            V[0] = V[1] = V[2] = forbidden;
        }
        else{
            double Et_a = 0.5 * M * omega_a*omega_a * x*x;
            double Et_b = 0.5 * M * omega_b*omega_b * x*x;
            V[0] = 2*A + Et_a;
            V[1] = Et_a;
            V[2] = lift + Et_b;
        }
    }
    bool is_angular()const {
        return false;
    }
    int dim()const{
        return 3;
    }
    int symmetry()const{
        return 0;
    }
    static const char* get_name(){
        return "3statel2";
    }

};


class ThreeStateStraight3Potential : public Potential {
public:
    double A = 0.02;
    double eps = 1;
    double omega_a = 0.005;
    double omega_b = 0.005;
    double M = 1000;
    double C = 0.005;
    double x0 = 0.5;

    double L = 2;
    double lift = 0;

    std::complex<double> V12 = {0.5, 0.1};
    std::complex<double> V13 = {0.1, 0.5};
    std::complex<double> V23 = {0.1, 0.1};

    const double forbidden = 1000;

    void init(const std::vector<double>& params){
        A = params[0]; 
        eps = params[1];
        omega_a = params[2];
        omega_b = params[3];
        M = params[4];
        C = params[5];
        L = params[6];
        x0 = params[7];

        if (params.size() > 8)
            lift = params[8];

        if (params.size() > 9){
            V12.real(params[9]);
            V12.imag(params[10]);
            V13.real(params[11]);
            V13.imag(params[12]);
            V23.real(params[13]);
            V23.imag(params[14]);
        }
    }

   void get_PE(double x, double y, std::vector<std::complex<double> >& V){
#define LOC(a,b) (a*3+b)
        if (x < -L/2 || x > L/2){
            std::fill(V.begin(), V.end(), 0.0);
            V[LOC(0,0)] = V[LOC(1,1)] = V[LOC(2,2)] = forbidden;
            return;
        }

        double Et_1 = 0.5 * M * omega_a*omega_a * x*x;
        double Et_2 = 0.5 * M * omega_b*omega_b * (x+x0)*(x+x0);
        double Et_3 = 0.5 * M * omega_b*omega_b * (x-x0)*(x-x0);
        double v = C * exp(-std::pow(y*eps, 2));

        V[LOC(0,0)] = Et_1 + A * tanh(y*eps) + A;
        V[LOC(1,1)] = Et_2 - A * tanh(y*eps) + A + lift;
        V[LOC(2,2)] = Et_3 - A * tanh(y*eps) + A + lift;

        V[LOC(0,1)] = v*V12;
        V[LOC(0,2)] = v*V13;
        V[LOC(1,2)] = v*V23;
        V[LOC(1,0)] = v*std::conj(V12);
        V[LOC(2,0)] = v*std::conj(V13);
        V[LOC(2,1)] = v*std::conj(V23);

#undef LOC
   }

    void get_PE_edgeA(double x, std::vector<double>& V){
        if (x < -L/2 || x > L/2){
            V[0] = V[1] = V[2] = forbidden;
        }
        else{
            double Et_1 = 0.5 * M * omega_a*omega_a * x*x;
            double Et_2 = 0.5 * M * omega_b*omega_b * (x+x0)*(x+x0);
            double Et_3 = 0.5 * M * omega_b*omega_b * (x-x0)*(x-x0);
            V[0] = Et_1;
            V[1] = 2*A + lift + Et_2;
            V[2] = 2*A + lift + Et_3;
        }
    }
    void get_PE_edgeB(double x, std::vector<double>& V){
        if (x < -L/2 || x > L/2){
            V[0] = V[1] = V[2] = forbidden;
        }
        else{
            double Et_1 = 0.5 * M * omega_a*omega_a * x*x;
            double Et_2 = 0.5 * M * omega_b*omega_b * (x+x0)*(x+x0);
            double Et_3 = 0.5 * M * omega_b*omega_b * (x-x0)*(x-x0);
            V[0] = 2*A + Et_1;
            V[1] = lift + Et_2;
            V[2] = lift + Et_3;
        }
    }
    bool is_angular()const {
        return false;
    }
    int dim()const{
        return 3;
    }
    int symmetry()const{
        return 0;
    }
    static const char* get_name(){
        return "3statel3";
    }

};


class LinearTwoStateBifurPotential1 : public Potential {
public:

    double A = 0.02;
    double C = 0.005;
    double W = 1.0;
    double eps = 2.5;
    double eps_c = 1.0;
    double omega = 0.01;
    double r = 2;
    double M = 1000;
    
    int control = 0;
    int phase_vary = 0;
    double lambda = 0;
    double yshift = 0;
    

    void init(const std::vector<double>& param){
        A = param[0];
        C = param[1];
        W = param[2];
        eps = param[3];
        eps_c = param[4];
        omega = param[5];
        r = param[6];
        M = param[7];

        if (param.size() > 8){
            control = int(param[8]);
            if (control >= 10 && control < 20){
                phase_vary = 1;
                control -= 10;
            }
            else if (control >= 20){
                phase_vary = 2;
                control -= 20;
            }
            if (param.size() > 9){
                yshift = param[9];
            }
        }
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        
        V[0] = 0.5*M*omega*omega*x*x + A * (exp(eps*(y - yshift)) - 1);
        V[3] = 0.5*M*omega*omega*(y < r ?
            std::pow(sqrt(x*x + (y-r)*(y-r)) - r, 2) :
            (x < -1e-10 ? (x+r)*(x+r) : (x > 1e-10 ? (x-r)*(x-r) : 1000)));

        double c;
        if (control == 0) c = x;
        else if (control == 1) c = std::abs(x);
        else if (control == 2) c = 1.0;

        if (phase_vary == 0){
            V[1] = C * exp(-eps_c*eps_c*y*y) * std::exp(std::complex<double>(0, W*x)) * c;
        }
        else if (phase_vary == 1){
            V[1] = exp(-eps_c*eps_c*y*y) * std::complex<double>(C* c, W);
        }
        else if (phase_vary == 2){
            V[1] = exp(-eps_c*eps_c*y*y) * std::sqrt(C*C*c*c + W*W);
        }
        V[2] = std::conj(V[1]);

    }

    void get_PE_edgeA(double x, std::vector<double>& V){

        V[0] = 0.5*M*omega*omega*x*x - A;
        V[1] = 100;
    }
    void get_PE_edgeB(double x, std::vector<double>& V){
        V[0] = 100;
        V[1] = 0.5*M*omega*omega*(x < -1e-10 ? (x+r)*(x+r) : (x > 1e-10 ? (x-r)*(x-r) : 1000));
    }

    int symmetry()const{
        return 0;
    }
    bool is_angular()const{
        return false;
    }
    static const char* get_name(){
        return "2statelx";
    }

};

class LinearTwoStatePseudoBifurPotential1 : public Potential {
public:

    double A = 0.02;
    double C = 0.005;
    double W = 1.0;
    double eps = 2.5;
    double eps_c = 1.0;
    double omega = 0.01;
    double r = 2;
    double M = 1000;
    
    int control = 0;
    int phase_vary = 0;
    double lambda = 0;
    double sigma = arma::datum::pi / 3;
    double slope = 0.1;
    

    void init(const std::vector<double>& param){
        A = param[0];
        C = param[1];
        W = param[2];
        eps = param[3];
        eps_c = param[4];
        omega = param[5];
        r = param[6];
        M = param[7];

        if (param.size() > 8){
            control = int(param[8]);
            if (control >= 10){
                phase_vary = 1;
                control -= 10;
            }
        }
        if (param.size() > 9){
            sigma = param[9];
            slope = param[10];
        }
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){

        double m_sig = atan2(x, r-y);
        
        V[0] = 0.5*M*omega*omega*x*x + A * (exp(eps*y) - 1);
        V[3] = 0.5*M*omega*omega*(y < r ?
            std::pow(sqrt(x*x + (y-r)*(y-r)) - r, 2) :
            (x < 0 ? (x+r)*(x+r) : 1000.0/(0.5*M*omega*omega)));
        if (x > 0 && y < r && m_sig > sigma){
            V[3] += slope * (m_sig - sigma);
        }

        double c;
        if (control == 0) c = x;
        else if (control == 1) c = std::abs(x);
        else if (control == 2) c = 1.0;

        if (!phase_vary){
            V[1] = C * exp(-eps*eps*y*y) * std::exp(std::complex<double>(0, W*x)) * c;
        }
        else{
            V[1] = exp(-eps*eps*y*y) * std::complex<double>(C* c, W);
        }
        V[2] = std::conj(V[1]);

    }

    void get_PE_edgeA(double x, std::vector<double>& V){

        V[0] = 0.5*M*omega*omega*x*x - A;
        V[1] = 1000;
    }
    void get_PE_edgeB(double x, std::vector<double>& V){
        V[0] = 1000;
        V[1] = 0.5*M*omega*omega*(x < 0 ? (x+r)*(x+r) : 100);
    }

    int symmetry()const{
        return 0;
    }
    bool is_angular()const{
        return false;
    }
    static const char* get_name(){
        return "2statelxp";
    }

};


class LinearPlanarCIPotential : public Potential {
public:

    int control = 0;
    int set_ci = 0; // yes
    double A = 0.02;
    double C = 0.01;
    double W = 0.5;
    double eps = 2.5;
    double epsc = 1.0;
    double omega1 = 0.005;
    double omega2 = 0.005;
    const double M = 1000;

    void init(const std::vector<double>& param) {
        control = param[0];
        set_ci = control / 10;
        control = control % 10;
        A = param[1];
        C = param[2];
        W = param[3];
        eps = param[4];
        epsc = param[5];
        omega1 = param[6];
        omega2 = param[7];
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        V[0] = 0.5*M*omega1*omega1*x*x + (control == 0 ? A*(tanh(eps*y)+1) : A*exp(eps*y));
        V[3] = 0.5*M*omega2*omega2*x*x + (control == 0 ? A*(tanh(-eps*y)+1) : A*exp(-eps*y));
        V[1] = C * exp(-eps*eps*y*y) * exp(std::complex<double>(0, W*x)) * (set_ci == 0 ? x : (set_ci == 2 ? std::abs(x) : 1.0));
        V[2] = std::conj(V[1]);
    }

    void get_PE_edgeA(double x, std::vector<double>& V){
        V[0] = 0.5*M*omega1*omega1*x*x;
        V[1] = 0.5*M*omega2*omega2*x*x + (control == 0 ? 2*A : 1);   
    }

    void get_PE_edgeB(double x, std::vector<double>& V){
        V[0] = 0.5*M*omega1*omega1*x*x + (control == 0 ? 2*A : 1);
        V[1] = 0.5*M*omega2*omega2*x*x;
    }

    int symmetry()const {
        return 0;
    }
    bool is_angular()const {
        return false;
    }
    static const char* get_name(){
        return "lpci";
    }
};

class LinearSTCrossingPotential : public Potential {
public:

    double L = 1.0;
    double A = 1.0;
    double eps = 5.0;
    double W = 1.0;
    double C1 = 1.0;
    double C2 = 1.0;
    double a = 1.0;
    double xm = 0.0, xp = 0.0;
    double L2 = 1.0;

    int control = 0;

    const double forbidden = 1000.0;

    void init(const std::vector<double>& param){
        L = param[0];
        A = param[1];
        eps = param[2];
        W = param[3];
        C1 = param[4];
        C2 = param[5];

        // switch control
        if (param.size() > 6){
            control = param[6];
        }

        if (control == 1 || control == 2){
            a = param[7];
            xm = param[8];
            xp = param[9];
        }
        if (control == 2){
            L2 = param[10];
        }
        else{
            L2 = L;
        }

    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        
#define LOC(a,b) (a*4+b)

        if (control == 2){
        double currentL = (L + L2)/2.0 + std::erf(eps*(y-xm)) * (L2 - L)/2.0;

        if (x < -currentL/2 || x > currentL/2){
            std::fill(V.begin(), V.end(), 0.0);
            V[LOC(0,0)] = V[LOC(1,1)] = V[LOC(2,2)] = V[LOC(3,3)] = forbidden;
            return;
        }
 
        }
        else{

        if (x < -L/2 || x > L/2){
            std::fill(V.begin(), V.end(), 0.0);
            V[LOC(0,0)] = V[LOC(1,1)] = V[LOC(2,2)] = V[LOC(3,3)] = forbidden;
            return;
        }
        }

        double theta;
        if (control == 1 || control == 2){
            theta = a * arma::datum::pi/2 * (std::erf(eps*(y - xm)) + 1.0)
            + (1.0 - a) * arma::datum::pi/2 * (std::erf(eps*(y - xp)) + 1.0);
        }
        else{
            theta = arma::datum::pi/2 * (std::erf(eps*y) + 1.0);
        }

        V[LOC(0,0)] = -cos(theta) * A + A;
        V[LOC(3,3)] = V[LOC(2,2)] = V[LOC(1,1)] =  cos(theta) * A + A;
        V[LOC(0,1)] = V[LOC(3,0)] = A* C1 * sin(theta) * std::exp(std::complex<double>(0, W*x));
        V[LOC(0,3)] = V[LOC(1,0)] = A *C1 * sin(theta) * std::exp(std::complex<double>(0, -W*x));
        V[LOC(0,2)] = V[LOC(2,0)] = A *C2 * sin(theta); 

        V[LOC(1,2)] = V[LOC(1,3)] = V[LOC(2,1)] = V[LOC(2,3)] = V[LOC(3,1)] = V[LOC(3,2)] = 0.0;
#undef LOC

    }

    void get_PE_edgeA(double x, std::vector<double>& V){

        if (x < -L/2 || x > L/2){
            V[0] = forbidden;
            V[1] = forbidden;
            V[2] = forbidden;
            V[3] = forbidden;
            return;
        }
        V[0] = 0;
        V[1] = V[2] = V[3] = 2*A;
    }
    void get_PE_edgeB(double x, std::vector<double>& V){

        if (x < -L2/2 || x > L2/2){
            V[0] = forbidden;
            V[1] = forbidden;
            V[2] = forbidden;
            V[3] = forbidden;
            return;
        }
        V[0] = 2*A;
        V[1] = V[2] = V[3] = 0;
    }

    int symmetry()const{
        return 0;
    }
    bool is_angular()const{
        return false;
    }
    int dim()const{
        return 4;
    }
    static const char* get_name(){
        return "stc";
    }

};



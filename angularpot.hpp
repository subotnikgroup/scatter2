
#include "potential.h"

class RightAnglePotential : public Potential {
public:

double B=5, D=4, W=0, C=0.15;

const double theta1=arma::datum::pi/2, theta2=0;

    void init(const std::vector<double>& param){
        B = param[0];
        D = param[1];
        W = param[2];
        C = param[3];
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){

#define calc_H_base(x_, y_) (tanh((y_) - B) - tanh((y_) + B) + tanh(x_) + D)
        double x1 = x*cos(theta1)+y*sin(theta1);
        double y1 = -x*sin(theta1)+y*cos(theta1);

        double x2 = x*cos(theta2)+y*sin(theta2);
        double y2 = -x*sin(theta2)+y*cos(theta2);

        V[0] = calc_H_base(x1, y1);
        V[3] = calc_H_base(x2, y2);
        V[1] = C * std::exp(std::complex<double>(0.0, W*std::sqrt(x*x + y*y)));
        V[2] = std::conj(V[1]);

#undef calc_H_base

    }
    void get_PE_edgeA(double x, std::vector<double>& V){
        V[0] = tanh(x-B)  - tanh(x+B) + D - 1;
        V[1] = tanh(x-B) - tanh(x+B) + D + 1;
    }

    int symmetry()const{
        return -1;
    }
    static const char* get_name(){
        return "rightangle";
    }
};


class AngularPotential : public Potential {
public:

    double r = 0.5;
    double R = 2.0;
    double A = 1.0;
    double eps = 5.0;
    double C = 1.0;
    double W = 1.0;
    const double theta1 = 0;
    const double theta2 = arma::datum::pi/2;
    const double forbidden = 1000.0;

    void init(const std::vector<double>& param){
        r = param[0];
        R = param[1];
        A = param[2];
        eps = param[3];
        C = param[4];
        W = param[5];
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        double x1 = x + R/2, y1 = y + R/2;
        double m_r = sqrt(x1*x1 + y1*y1);
        double m_theta = std::atan2(y1, x1);
        
        if (m_r < r || m_r >= R){
            V[0] = forbidden;
            V[3] = forbidden;
            V[1] = V[2] = 0.0;
            return;
        }

        V[0] = tanh((m_theta - (theta1+theta2)/2)*eps) * A + A;
        V[3] = tanh(-(m_theta - (theta1+theta2)/2)*eps) * A + A;
        V[1] = C * exp(-eps*eps*std::pow(m_theta - (theta1+theta2)/2, 2)) * std::exp(std::complex<double>(0, W*m_r));
        V[2] = std::conj(V[1]);

    }

    void get_PE_edgeA(double x, std::vector<double>& V){

        if (x + R/2 < r || x + R/2 >= R){
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
    static const char* get_name(){
        return "angle";
    }
};

class AngularVectorWPotential : public Potential {
public:

    double r = 0.5;
    double R = 2.0;
    double A = 1.0;
    double eps = 5.0;
    double C = 1.0;
    double W = 1.0;
    double ltheta = arma::datum::pi/2;

    const double theta1 = 0;
    const double theta2 = arma::datum::pi/2;
    const double forbidden = 1000.0;

    void init(const std::vector<double>& param){
        r = param[0];
        R = param[1];
        A = param[2];
        eps = param[3];
        C = param[4];
        W = param[5];
        ltheta = param[6];
    }
    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        double x1 = x + R/2, y1 = y + R/2;
        double m_r = sqrt(x1*x1 + y1*y1);
        double m_theta = std::atan2(y1, x1);

        if (m_r < r || m_r >= R){
            V[0] = forbidden;
            V[3] = forbidden;
            V[1] = V[2] = 0.0;
            return;
        }
       
        V[0] = tanh((m_theta - (theta1+theta2)/2)*eps) * A + A;
        V[3] = tanh(-(m_theta - (theta1+theta2)/2)*eps) * A + A;
        V[1] = C * exp(-eps*eps*std::pow(m_theta - (theta1+theta2)/2, 2)) * 
            std::exp(std::complex<double>(0, W*(cos(ltheta)*x1 + sin(ltheta)*y1)));
        V[2] = std::conj(V[1]);

    }

    void get_PE_edgeA(double x, std::vector<double>& V){

        if (x + R/2 < r || x + R/2 >= R){
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
    static const char* get_name(){
        return "anglel";
    }
};

class AngularPotentialEqualThick : public Potential {
public:

    double r = 0.5;
    double R = 2.0;
    double A = 1.0;
    double eps = 5.0;
    double C = 1.0;
    double W = 1.0;
    const double theta1 = 0;
    const double theta2 = arma::datum::pi/2;
    const double forbidden = 1000.0;

    void init(const std::vector<double>& param){
        r = param[0];
        R = param[1];
        A = param[2];
        eps = param[3];
        C = param[4];
        W = param[5];
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        double x1 = x + R/2, y1 = y + R/2;
        double m_r = sqrt(x1*x1 + y1*y1);
        double m_theta = std::atan2(y1, x1);
        
        if (m_r < r || m_r >= R){
            V[0] = forbidden;
            V[3] = forbidden;
            V[1] = V[2] = 0.0;
            return;
        }

        V[0] = tanh(m_r * std::sin(m_theta - (theta1+theta2)/2)*eps) * A + A;
        V[3] = tanh(-m_r * std::sin(m_theta - (theta1+theta2)/2)*eps) * A + A;
        V[1] = C * exp(-eps*eps*std::pow(m_r * std::sin(m_theta - (theta1+theta2)/2), 2)) * std::exp(std::complex<double>(0, W*m_r));
        V[2] = std::conj(V[1]);

    }

    void get_PE_edgeA(double x, std::vector<double>& V){

        if (x + R/2 < r || x + R/2 >= R){
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
    static const char* get_name(){
        return "angleet";
    }
};


class AngularPotential2StateMag : public Potential {
public:

    double r = 4.0;
    double R = 8.0;
    double A = 0.02;
    double eps = 5.0;
    double C = 0.005;
    double W = 1.0;

    double C2 = 0.1;
    double eps2 = 8.0;
    double theta_c2_off = 0; // theta_c2 = theta_c2_off + theta_c
    double dtheta = 0.78539816; // pi/4

    const double theta_c = arma::datum::pi/4;
    const double forbidden = 1000.0;

    double theta_c2;

    void init(const std::vector<double>& param){
        r = param[0];
        R = param[1];
        A = param[2];
        eps = param[3];
        C = param[4];
        W = param[5];

        C2 = param[6];
        eps2 = param[7];
        theta_c2_off = param[8];
        dtheta = param[9];

        theta_c2 = theta_c2_off + theta_c;
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        double x1 = x + R/2, y1 = y + R/2;
        double m_r = sqrt(x1*x1 + y1*y1);
        double m_theta = std::atan2(y1, x1);
        
#define LOC(a,b) (a*4+b)
        if (m_r < r || m_r >= R){
            std::fill(V.begin(), V.end(), 0.0);
            V[LOC(0,0)] = V[LOC(1,1)] = V[LOC(2,2)] = V[LOC(3,3)] = forbidden;
            return;
        }
        V[LOC(1,1)] = V[LOC(0,0)] = tanh((m_theta - theta_c)*eps) * A + A;
        V[LOC(2,2)] = V[LOC(3,3)] = tanh(-(m_theta - theta_c)*eps) * A + A;
        
        std::complex<double> r = C * exp(-eps*eps*std::pow(m_theta - theta_c, 2)) * std::exp(std::complex<double>(0, W*m_r));
        double u = C2/2 * (std::tanh(eps2*(m_theta - theta_c2 + dtheta/2)) + std::tanh(-eps2*(m_theta - theta_c2 - dtheta/2)));
        V[LOC(2,0)] = V[LOC(1,3)] = std::conj(r);
        V[LOC(0,2)] = V[LOC(3,1)] = r;
        V[LOC(0,1)] = V[LOC(1,0)] = u;
        V[LOC(2,3)] = V[LOC(3,2)] = 0.0;
#undef LOC
    }

    void get_PE_edgeA(double x, std::vector<double>& V){

        if (x + R/2 < r || x + R/2 >= R){
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

    int dim()const{
	return 4;
    }
    static const char* get_name(){
        return "anglem";
    }

};

class AngularPotential2StateMagFull : public Potential {
public:

    double r = 4.0;
    double R = 8.0;
    double A = 0.02;
    double eps = 5.0;
    double C = 0.005;
    double W = 1.0;

    double B = 0.0;
    double thetaB = 0.0;
    double phiB = 0.0;

    const double theta_c = arma::datum::pi/4;
    const double forbidden = 1000.0;

    void init(const std::vector<double>& param){
        r = param[0];
        R = param[1];
        A = param[2];
        eps = param[3];
        C = param[4];
        W = param[5];
        B = param[6];
        thetaB = param[7];
        phiB = param[8];
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        double x1 = x + R/2, y1 = y + R/2;
        double m_r = sqrt(x1*x1 + y1*y1);
        double m_theta = std::atan2(y1, x1);
        
#define LOC(a,b) (a*4+b)
        if (m_r < r || m_r >= R){
            std::fill(V.begin(), V.end(), 0.0);
            V[LOC(0,0)] = V[LOC(1,1)] = V[LOC(2,2)] = V[LOC(3,3)] = forbidden;
            return;
        }
        double dA = tanh((m_theta - theta_c)*eps) * A; 
        V[LOC(0,0)] = A + dA + B;
        V[LOC(1,1)] = A + dA - B;
        V[LOC(2,2)] = A - dA + B;
        V[LOC(3,3)] = A - dA - B;
        
        double C0 =  C * exp(-eps*eps*std::pow(m_theta - theta_c, 2));
        std::complex<double> r = C0 * std::complex<double>(cos(W*m_r), cos(thetaB)*sin(W*m_r));
        std::complex<double> s = C0 * sin(W*m_r) * sin(thetaB) * std::exp(std::complex<double>(0, -phiB-arma::datum::pi/2)); 

        V[LOC(0,2)] = V[LOC(3,1)] = r;
        V[LOC(2,0)] = V[LOC(1,3)] = std::conj(r);
        V[LOC(0,3)] = s;
        V[LOC(1,2)] = -std::conj(s);
        V[LOC(2,1)] = -s;
        V[LOC(3,0)] = std::conj(s);
        V[LOC(0,1)] = V[LOC(1,0)] = V[LOC(2,3)] = V[LOC(3,2)] = 0.0;
#undef LOC
    }

    void get_PE_edgeA(double x, std::vector<double>& V){

        if (x + R/2 < r || x + R/2 >= R){
            V[0] = forbidden;
            V[1] = forbidden;
            V[2] = forbidden;
            V[3] = forbidden;
            return;
        }
        V[0] = B;
        V[1] = -B;
        V[2] = 2*A + B;
        V[3] = 2*A - B;
    }

    void get_PE_edgeB(double x, std::vector<double>& V){
        if (x + R/2 < r || x + R/2 >= R){
            V[0] = forbidden;
            V[1] = forbidden;
            V[2] = forbidden;
            V[3] = forbidden;
            return;
        }
        V[0] = 2*A + B;
        V[1] = 2*A - B;
        V[2] = B;
        V[3] = -B;

    }

    int symmetry()const{
        return 0;
    }

    int dim()const{
	return 4;
    }
    static const char* get_name(){
        return "angleb";
    }
};




class AngularPotential4State : public Potential {
public:

    double r = 0.5;
    double R = 2.0;
    double A = 1.0;
    double eps = 5.0;
    double C = 1.0;
    double W = 1.0;

    double C2 = 0.1;
    double W2 = 1.0;

    const double theta1 = 0;
    const double theta2 = arma::datum::pi/2;
    const double forbidden = 1000.0;

    void init(const std::vector<double>& param){
        r = param[0];
        R = param[1];
        A = param[2];
        eps = param[3];
        C = param[4];
        W = param[5];
        C2 = param[6];
        W2 = param[7];
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        double x1 = x + R/2, y1 = y + R/2;
        double m_r = sqrt(x1*x1 + y1*y1);
        double m_theta = std::atan2(y1, x1);
        
#define LOC(a,b) (a*4+b)
        if (m_r < r || m_r >= R){
            std::fill(V.begin(), V.end(), 0.0);
            V[LOC(0,0)] = V[LOC(1,1)] = V[LOC(2,2)] = V[LOC(3,3)] = forbidden;
            return;
        }
        V[LOC(1,1)] = V[LOC(0,0)] = tanh((m_theta - (theta1+theta2)/2)*eps) * A + A;
        V[LOC(2,2)] = V[LOC(3,3)] = tanh(-(m_theta - (theta1+theta2)/2)*eps) * A + A;
        
        std::complex<double> r = C * exp(-eps*eps*std::pow(m_theta - (theta1+theta2)/2, 2)) * std::exp(std::complex<double>(0, W*m_r));
        // largerst in the middle, with a phase depend on angle
        std::complex<double> s = C2 * exp(-eps*eps*std::pow(m_theta - (theta1+theta2)/2, 2)) * std::exp(std::complex<double>(0, W2*m_r));
        V[LOC(2,0)] = V[LOC(1,3)] = std::conj(r);
        V[LOC(0,2)] = V[LOC(3,1)] = r;
        V[LOC(3,0)] = std::conj(s);
        V[LOC(0,3)] = s;
        V[LOC(2,1)] = -s;
        V[LOC(1,2)] = -std::conj(s);
        V[LOC(0,1)] = V[LOC(1,0)] = V[LOC(2,3)] = V[LOC(3,2)] = 0.0;
#undef LOC
    }

    void get_PE_edgeA(double x, std::vector<double>& V){

        if (x + R/2 < r || x + R/2 >= R){
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

    int dim()const{
	return 4;
    }
    static const char* get_name(){
        return "angle4";
    }

};

class TestPotential : public Potential {
public:

    double A = 1.0;
    double eps = 5.0;
    double C = 1.0;
    double r = 0.5;
    double R = 2.0;
    const double theta1 = 0, theta2 = arma::datum::pi/2;
    const double forbidden = 100.0;

    void init(const std::vector<double>& param){
        A = param[0];
        eps = param[1];
        C = param[2];
        r = param[3];
        R = param[4];
    }

    void get_PE(double x, double y, std::vector<std::complex<double> >& V){
        double x1 = x + R/2, y1 = y + R/2;
        double m_r = sqrt(x1*x1 + y1*y1);
        double m_theta = std::atan2(y1, x1);
        
         if (m_r < r || m_r >= R){
            V[0] = forbidden;
            V[3] = forbidden;
            V[1] = V[2] = 0.0;
            return;
        }

        V[0] = 0;
        V[3] = 1;
        V[2] = V[1] = C * exp(-2*eps*std::pow(m_theta - (theta1+theta2)/2, 2));
    }
    void get_PE_edgeA(double x, std::vector<double>& V){

        if (x + R/2 < r || x + R/2 >= R){
            V[0] = forbidden;
            V[1] = forbidden;
            return;
        }
        V[0] = 0;
        V[1] = 1;
    }

    int symmetry()const{
        return 1;
    }
    static const char* get_name(){
        return "test";
    }
};

class ThreeStateAngularPotential : public Potential {
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
    const double theta0 = arma::datum::pi/4;

    void init(const std::vector<double>& params){
        A = params[0]; 
        eps = params[1];
        omega_a = params[2];
        omega_b = params[3];
        M = params[4];
        C = params[5];
        L = params[6];
        r0 = params[7];
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
        double x1 = x+r0, y1 = y+r0;
        double m_r = sqrt(x1*x1 + y1*y1);
        double m_theta = std::atan2(y1, x1);

        if (m_r < r0-L/2 || m_r > r0+L/2){
            std::fill(V.begin(), V.end(), 0.0);
            V[LOC(0,0)] = V[LOC(1,1)] = V[LOC(2,2)] = forbidden;
            return;
        }

        double Et_a = 0.5 * M * omega_a*omega_a * (m_r - r0)*(m_r - r0);
        double Et_b = 0.5 * M * omega_b*omega_b * (m_r - r0)*(m_r - r0);
        double v = C * exp(-std::pow((m_theta - theta0)*eps, 2));

        V[LOC(0,0)] = Et_a + A * tanh((m_theta - theta0)*eps) + A;
        V[LOC(1,1)] = Et_a - A * tanh((m_theta - theta0)*eps) + A;
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
            V[2] = 2*A + Et_b + lift;
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
            V[2] = 2*A + Et_b + lift;
        }
    }
    bool is_angular()const {
        return true;
    }
    int dim()const{
        return 3;
    }
    int symmetry()const{
        return 0;
    }
    static const char* get_name(){
        return "3statea";
    }

};

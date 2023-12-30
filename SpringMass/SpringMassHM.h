// HM for harmonic oscillator
#include <iostream>
#include <cmath>
#include <fstream>

class SpringMassHM {
    public:
        const double dt; // time step, should be .01 for usual use
        const double gamma; // damping constant
        const double w; //omega aka sqrt(k/m)
        const double k; //spring constant
        const double m; //mass of spring
        

        double get_t_final() {return t_final;}
        double get_num_timesteps() {return num_timesteps;}
        void eq_solver(double x, double v);
        SpringMassHM(double dt, double k, double m, double t_final, double gamma);
        ~SpringMassHM();

        //Need to Implement
        void change_inputs(double dt, double w, double k, double m, double t_final, double gamma); 
    
    private:
        double t_final; //end condition for simulation
        double a; // acceleration
        int num_timesteps;
        double *t, *x, *v, *dx_dt, *d2x_dt;
        double* RK4(double dt, int i);
        void euler(double dt, int i);

};
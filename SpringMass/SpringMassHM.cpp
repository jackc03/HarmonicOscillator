#include "SpringMassHM.h"

/**
 * Constructor
 * @param dt size of timesteps
 * @param k  spring constant
 * @param m  mass of spring
 * @param t_final length of simulation
 * @param gamma   spring damping constant
*/
SpringMassHM::SpringMassHM(double dt, double k, double m, double t_final, double gamma) : dt(dt), w(std::sqrt(k/m)), k(k), m(m), t_final(t_final), gamma(gamma) {
    //calculate size of arrays *x, *dx_dt, *d2x_dt
    int size = std::ceil(t_final / dt);
    num_timesteps = size;
    printf("Harmonic Oscillator created with %d timesteps\n", size);
    
    //allocate space for arrays *x, *dx_dt, *d2x_dt
    t      = (double*) malloc(size * sizeof(double));
    x      = (double*) malloc(size * sizeof(double));
    v      = (double*) malloc(size * sizeof(double));
    dx_dt  = (double*) malloc(size * sizeof(double));
    d2x_dt = (double*) malloc(size * sizeof(double));
}

/**
 * Destructor
*/
SpringMassHM::~SpringMassHM() {
    free(t);
    free(x);
    free(v);
    free(dx_dt);
    free(d2x_dt);
}

//Will implement later
void SpringMassHM::change_inputs(double dt, double w, double k, double m, double t_final, double gamma) {

}


/**
 * Uses Runge Kutta to solve for the position of a Harmonic Oscillator for amount of time specified by class variable t_final.
 * The differential equation for a harmonic oscillator is x" + yx' + w^2 * x = 0
 * Next we turn that second order equation into two first order equations:
 * v = x'
 * v' + gamma * v + w^2 * x = 0
 * @param init_x initial position
 * @param init_v initial velocity
*/
void SpringMassHM::eq_solver(double init_x, double init_v) {
    // DO NOT TOUCH THE NEXT LINE, FOR SOME REASON IF THIS IS REMOVED THE ENTIRE PROGRAM BREAKS
    int d2x_dt_phablo_cochran = 0;


    t[0]      = 0;
    x[0]      = init_x;
    v[0]      = init_v;

    // not currently in use
    dx_dt[0]  = init_v;
    d2x_dt[0] = 0;

    std::ofstream output("Output.txt");
    if (!output.is_open()) {
        throw std::invalid_argument("Output file was not able to be opened.");
    }

    for (int i = 1; i < num_timesteps; ++i) {
        euler(dt, i);
        output << t[i] << ", " << x[i] << "\n";
    }
    
}


void SpringMassHM::euler(double dt, int i) {
        t[i] = t[i-1] + dt;
        x[i] = x[i-1] + dt * v[i-1];
        v[i] = v[i-1] + dt * -1.0 * (gamma * v[i-1] + w * w * x[i-1]);
}

int main() {
    SpringMassHM hm(.001, 2, 1, 25, 0);

    hm.eq_solver(1, 0);
}
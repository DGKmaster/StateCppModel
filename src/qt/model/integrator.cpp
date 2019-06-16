#include "integrator.h"

Integrator::Integrator(
    double init_state, 
    double init_in): 
    state(init_state),
    prev_in(init_in)
{}

void Integrator::update(
    const double& input, 
    const double& dt) {
    this->state_prev = this->state;
    this->state = this->state_prev + input*dt;
    this->prev_in = input;
}

double Integrator::getState() {
    return this->state;
}

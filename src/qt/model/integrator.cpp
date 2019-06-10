#include "integrator.h"

Integrator::Integrator(double init_state, double init_in): 
    state(init_state),
    prev_in(init_in)
{}

void Integrator::update(double input, double dt) {
    this->state_prev = this->state;
    this->state = state_prev + input*dt;
    this->prev_in = input;
}

double Integrator::getState_prev() {
    return this->state_prev;
}

double Integrator::getState() {
    return this->state;
}

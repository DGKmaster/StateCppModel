class Integrator {

public:
	Integrator(double init_state, double init_in): 
        state(init_state),
        prev_in(init_in)
    {}

	void update(double input, double dt) {
        this->state_prev = this->state;
        this->state = state_prev + input*dt;
        this->prev_in = input;
	}

    double getState_prev() {
		return this->state_prev;
	}

	double getState() {
		return this->state;
	}

private:
	double state;
    double state_prev;
	double prev_in = 0;
};

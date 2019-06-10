class Integrator {

public:
	Integrator(double init_state, double init_in);

	void update(double input, double dt);

    double getState_prev();

	double getState();

private:
	double state;
    double state_prev;
	double prev_in = 0;
};

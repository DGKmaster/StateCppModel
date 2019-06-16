class Integrator {

public:
	Integrator(double init_state, double init_in);

	void update(
		const double& input, 
		const double& dt);

	double getState();

private:
	double state;
    double state_prev;
	double prev_in = 0;
};

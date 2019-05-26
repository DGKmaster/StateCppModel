
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iterations = 1000;
input = 1;

A = [0, 1, 0; 0, 0, 1; -1.5, -5, -2];
T = 0.01;
B = [0; 0; 1];
C = [0.5 0 0];

Ad = expm(A*T);
Bd = inv(A)*(Ad - eye(3))*B;
Cd = C;

x = zeros(3, iterations);
y = zeros(-1, iterations);

for i=2:1000
    x(:, i) = Ad*x(:, i-1) + Bd*input;
    y(:, i) = Cd*x(:, i);
end

plot(y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x_iminus1_up x_iminus1_down omega T x_init_up x_init_down;
x_iminus1 = [x_iminus1_up; x_iminus1_down];
x_init = [x_init_up; x_init_down];
A = [0 -1;omega^2 0];

MatExp = inv(eye(2)+A*T)*x_iminus1;
x_i = x_iminus1 - A* MatExp*T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

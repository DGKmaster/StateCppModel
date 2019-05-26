% A = [0, 1, 0; 0, 0, 1; -1.5, -5, -2];
A = [0, 1, 0; 0, 0, 1; -0.2, -1, -1];
B = [0; 0; 1];
C = [2 1 0];

iterations = 1000;
input = 1;

x = zeros(3, iterations);
y = zeros(1, iterations);
for i=2:1000
    x(:, i) = A*x(:, i-1) + B*input;
    y(:, i) = C*x(:, i);
end

plot(y);

clear
clc

h1 = [2; 6];
h2 = [-2; 4];
lambda1 = 10;
lambda2 = 0;
I = [1 0; 0 1];
variance = 1;

matrix = I + lambda1/variance*h1*transpose(h1) + lambda2/variance*h2*transpose(h2);
inverse_matrix = inv(matrix);

h1_new = inverse_matrix*h1/norm(inverse_matrix)
h2_new = inverse_matrix*h2/norm(inverse_matrix)
plot([0 h1(1)], [0 h1(2)]);
hold on
plot([0 h2(1)], [0 h2(2)]);
hold on
plot([0 h1_new(1)], [0 h1_new(2)]);
hold on
plot([0 h2_new(1)], [0 h2_new(2)]);
legend('h1', 'h2', 'h1new', 'h2new');
hold on
xlim([-20 20]);
ylim([-20 20]);
transpose(h1)*h2_new
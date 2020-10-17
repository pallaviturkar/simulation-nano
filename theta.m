V= 1;
R=1.5;
x = 39;
fun = @(x)((3*V*(sin(x)/R)^3)/((1-cos(x))^2.*(2+cos(x))));
x = fzero(fun,39);


c1 = (3*V)/(pi*R^3);
syms x
eqn = (((1-cos(x))^2)*(2+cos(x)))/(sin(x))^3 == c1;
S = solve(eqn, x, 'Real', true);
theta1 = double(S);


%x = fzero(fun,x0)
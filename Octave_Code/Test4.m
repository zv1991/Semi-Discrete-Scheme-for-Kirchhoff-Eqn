%%% Test %%%
% shifted Legendre polynomials of degrees 0, 1, 2, 3 %
%%% ----------- %%%
osc = 1;
%%% ----------- %%%

P0    = @(x) 1;
P1    = @(x) (2 / ell) * x - 1;
P2    = @(x) 0.5 * (3 * ((2 / ell) * x - 1) .* ((2 / ell) * x - 1) - 1);
P3    = @(x) 0.5 * (5 * ((2 / ell) * x - 1) .* ((2 / ell) * x - 1) .* ...
((2 / ell) * x - 1) - 3 .* ((2 / ell) * x - 1));
Phi1  = @(x) 0.5 * sqrt(ell / 3) * (P2(x) - P0(x));
% exact solution %
u     = @(x,t) (t / T) * Phi1(x);
% initial conditions %
psi0  = @(x) u(x,0);
psi1  = @(x) (1 / T) * Phi1(x);

% for test %
du    = @(x,t) (t / T) * sqrt(3 / ell) * P1(x);
intdu = @(t) (t * t) / (T * T);

% time-dependent coefficients %
alpha = @(t) 1;
beta  = @(t) 1;
% right-hand side of the equation %
f     = @(x,t) -(2 * t * sqrt(3)) / (ell * T * sqrt(ell)) * (alpha(t) + ...
intdu(t) * beta(t)) * P0(x);
%%%%%%%%%%%%

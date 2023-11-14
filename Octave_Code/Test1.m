%%% Test %%%
%  The oscillation number of test function %

##osc = 1;
##osc = 3;
##osc = 7;
##osc = 11;

% functions %
u1    = @(x) sin((pi * x) / ell);
u2    = @(t) sin((osc * pi * t) / T);
% exact solution %
u     = @(x,t) u1(x) * u2(t);
% initial conditions %
psi0  = @(x) u(x,0);
psi1  = @(x) ((osc * pi) / T) * u1(x);

% for test %
du    = @(x,t) (pi / ell) * cos((pi * x) / ell) * u2(t);
intdu = @(t) ((pi * pi) / (2 * ell)) * u2(t) * u2(t);

% time-dependent coefficients %
alpha = @(t) 2 + sin((10 * pi * t) / T);
beta  = @(t) 1 + t .* t;
% right-hand side of the equation %
f     = @(x,t) pi * pi * (-(osc * osc) / (T * T) + (1 / (ell * ell)) * ...
(alpha(t) + beta(t) * intdu(t))) * u(x,t);
%%%%%%%%%%%%

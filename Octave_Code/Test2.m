%%% Test %%%
%  The oscillation number of test function %
osc = 5;
% osc = 17;
% functions %
u1    = @(x) sin((osc * pi * x) / ell);
u2    = @(t) exp(pi * t);
% exact solution %
u     = @(x,t) u1(x) * u2(t);
% initial conditions %
psi0  = @(x) u1(x);
psi1  = @(x) pi * u1(x);

% for test %
du    = @(x,t) (osc * pi / ell) * cos((osc * pi * x) / ell) * u2(t);
intdu = @(t) ((osc * osc * pi * pi) / (2 * ell)) * u2(t) * u2(t);

% time-dependent coefficients %
alpha = @(t) 2 + sin((10 * pi * t) / T);
beta  = @(t) 1 + t .* t;
% right-hand side of the equation %
f     = @(x,t) pi * pi * (1 + (osc * osc) / (ell * ell) * (alpha(t) + ...
beta(t) * intdu(t))) * u(x,t);
%%%%%%%%%%%%

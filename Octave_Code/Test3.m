%%% Test %%%
%  The oscillation number of test function %
osc = 5;
##osc = 17;

amp = 1; %% Amplitude %%

% functions %
u1    = @(x) sin((osc * pi * x) / ell);
u2    = @(t) sqrt(1 + t);
% exact solution %
u     = @(x,t) amp * u1(x) * u2(t);
% initial conditions %
psi0  = @(x) amp * u1(x);
psi1  = @(x) 0.5 * amp * u1(x);

% for test %
du    = @(x,t) (amp * osc * pi) / ell * u2(t) * cos((osc * pi * x) / ell);
intdu = @(t) (amp * amp * osc * osc * pi * pi) / (2 * ell) * (1 + t);

% time-dependent coefficients %
coeff_alpha = (ell * ell * ell) / 2; %% Coefficient c1 in the handwriting %%
alpha0      = (ell * ell * ell - coeff_alpha) / (4 * ell * osc * osc * pi * pi);
beta0       = coeff_alpha / (2 * amp * amp * osc * osc * osc * osc * ...
pi * pi * pi * pi);
alpha = @(t) alpha0 / ((1 + t) .* (1 + t));
beta  = @(t) beta0 / ((1 + t) .* (1 + t) .* (1 + t));
% right-hand side of the equation %
f     = @(x,t) 0;
%%%%%%%%%%%%

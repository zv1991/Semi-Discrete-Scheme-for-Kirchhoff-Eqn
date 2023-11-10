%%% Test without knowing exact solution %%%

%%% Preliminary %%%
LegPoly = @(n, x) LegendrePoly (n, (2 * x) / ell - 1);
Phi     = @(x) LegPoly(33,x) - LegPoly(31,x);

%%%%%%%%%%%%%%%%%%%%

amp = 1; %% height of the curve's peak %%
lam = (4 * 4) / (ell * ell); %% (c * c) / (ell * ell) %%

% initial conditions %
% psi0  = @(x) amp * exp(-lam * (2 * x - ell) .* (2 * x - ell)) .* ...
% sin(8 * pi * x);
psi0  = @(x) amp * exp(-lam * (2 * x - ell) .* (2 * x - ell)) .* Phi(x);
psi1  = @(x) 0;

% time-dependent coefficients %
alpha = @(t) 1;
beta  = @(t) 1;

% right-hand side of the equation %
f     = @(x,t) 0;

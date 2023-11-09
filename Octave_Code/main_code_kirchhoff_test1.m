clc
clearAllMemoizedCaches
clearvars
clear all
close all

format longG

%%% Length of spatial and temporal intervals %%%
ell = 1;
T   = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Please input numbers, e.g. 1, 2, ...');
problem = 1;
run(sprintf('%s%d.m','Test',problem));

m = 1e3; % Division number of the spatial interval %
%%%  Discretization of the spatial interval %%%
x   = linspace(0,ell,m + 1);
h   = ell / m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eigcos = cos(pi / m); %%% Cosine appeared to eigenvalues %%%

%%% The second-order derivative of u_0(x) %%%
d2psi0 = zeros(1,length(x));
d2cfd  = [-1, 16, -30, 16, -1]; %%% fourth-order accurate scheme %%%
for i = 2:m
  temp = linspace(x(i - 1),x(i + 1),5);
  % temp = x(i - 1):h/2:x(i + 1);
  d2psi0(i) = dot(d2cfd,psi0(temp)) / (3 * h * h);
endfor
clear i

%%% Tridiagonal system coefficients %%%
b   = zeros(1,length(x) - 2); % diagonal entries %
phi = zeros(1,length(x)); % right-hand side %
w   = zeros(1,length(x)); % Solution of tridiagonal system %

coeffsys2 = h * h; % tridiagonal system coefficients %

error      = []; %%% error initialization %%%
n1         = [];
tau1       = [];
time       = [];
max_CondN2 = []; %%% Maximum Condition Number of Matrix %%%

n_max = m;
% n_max = 1024;

elaps_t = tic;
for j = 1e1:1e1:n_max
  tStart = cputime; % For CPU time determination %
  n      = j;
  n1     = [n1;n];
  %%% Discretization of the temporal interval %%%
  t      = linspace(0,T,n + 1);
  tau    = T / n;
  tau1   = [tau1;tau];
  CondN2 = []; %%% Condition Number of Matrix %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%% Initialization of approximation of the functions u(x,t) on the grid %%%
  val_u      = zeros(length(t),length(x)); % exact values of u(x,t) %
  app_u      = zeros(length(t),length(x)); % approximate values of u_k(x) %

  coeffsys1  = 1 / (tau * tau); % tridiagonal system coefficients %

  k = 1;
  while k <= n + 1
    fprintf('Layer, k = %d\n', k);
    switch k
      case 1
        app_u(k,:) = psi0(x);
        tempdiff1  = FinDiff (h, app_u(k,:));
        tempdiff1  = tempdiff1 .* tempdiff1;
        %%% q integral %%%
        q1         = alpha(t(k)) + beta(t(k)) * SimpsRule (h, tempdiff1);
        condN2     = NaN; %%% Condition Number %%%
        CondN2     = [CondN2;condN2]; %%% Condition Number %%%
        val_u(k,:) = u(x,t(k));
      case 2
        % app_u(k,:) = u(x,t(k));
        app_u(k,:) = app_u(k - 1,:) + tau * (psi1(x) + 0.5 * tau * ...
        (f(x,t(k - 1)) + q1 * d2psi0));
        tempdiff1  = FinDiff (h, app_u(k,:));
        tempdiff1  = tempdiff1 .* tempdiff1;
        %%% q integral %%%
        q2 = alpha(t(k)) + beta(t(k)) * SimpsRule (h, tempdiff1);
        val_u(k,:) = u(x,t(k));
      otherwise
        p      = coeffsys1 / q2;
        coeff1 = coeffsys2 * p;
        coeff2 = coeff1 / 6;
        v      = tau * tau * f(x,t(k - 1)) + 2 * app_u(k - 1,:);
        b(:)   = -2 * (1 + coeff1 * (1 + coeff2));
        for i = 2:m
          phi(i) = -coeff2 * (v(i - 1) + 2 * (5 + coeff1) * v(i) + v(i + 1));
        endfor
        clear i
        w(2:m)     = tridiagonal_algotirhm (b, phi(2:m));
        app_u(k,:) = w - app_u(k - 2,:);
        val_u(k,:) = u(x,t(k));
        tempdiff1  = FinDiff (h, app_u(k,:));
        tempdiff1  = tempdiff1 .* tempdiff1;
        q3         = alpha(t(k)) + beta(t(k)) * SimpsRule (h, tempdiff1);
        %%% Condition Number %%%
        condN2     = 1 + (2 * eigcos) / (1 + coeff1 * (1 + coeff2) - eigcos);
        %%%%%%%%%%%%%%%%%%%%%%%%
        CondN2     = [CondN2;condN2]; %%% Condition Number %%%
        q1 = q2;
        q2 = q3;
    endswitch
    k = k + 1;
  endwhile
  CondN2 = [CondN2;NaN]; %%% Condition Number %%%

  tEnd = cputime - tStart; % For CPU time determination %
  time = [time;tEnd];

  error1     = max(max(abs(app_u - val_u)));
  error      = [error;error1];
  max_CondN2 = [max_CondN2;max(CondN2(2:(end - 1)))];
endfor

total_t = toc(elaps_t);
if total_t == 1
  fprintf('Elapsed time is %.15f second\n',total_t);
elseif total_t < 60
  fprintf('Elapsed time is %.15f seconds\n',total_t);
else
  total_t1 = total_t / 60;
  switch total_t1
    case 1
      fprintf('Elapsed time is %.15f minute\n',total_t1);
    otherwise
      fprintf('Elapsed time is %.15f minutes\n',total_t1);
  endswitch
endif

n10     = log10(n1);
tau10   = log10(tau1);
error10 = log10(error);
n2      = log2(n1);
tau2    = log2(tau1);
error2  = log2(error);

%%% Simple Linear Regression %%%

S_xy = SLR_sampcov (tau10, error10);
S_xx = SLR_sampcov (tau10, tau10);
S_yy = SLR_sampcov (error10, error10);

[beta1, beta0, y2, res, min_res, max_res, r_squared] = SLR (tau10, error10);

% run("SimpleLinReg.m");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirname  = sprintf('Test%d',problem);
if not(isfolder(dirname))
  mkdir(dirname); % Folder creator %
endif

run("csvtable_slr.m");
run("csvtable_log10.m");
run("csvtable_time.m");
run("csvtable_maxCondN2.m");

%%% generate mat file %%%
filemat  = sprintf('log_log_graph_Test%d_osc=%d.mat',problem,osc);
matfile  = fullfile(dirname,filemat);
save(matfile,"h","n1","error","tau1","tau10","tau2","time","n10","error10",...
"max_CondN2","n2","error2","S_xy","S_xx","S_yy","beta1","beta0","y2","res",...
"min_res","max_res","r_squared","total_t");
%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except matfile
close all

load(matfile);

clearAllMemoizedCaches
clear matfile

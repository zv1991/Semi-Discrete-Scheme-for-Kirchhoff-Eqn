clc
clearAllMemoizedCaches
clearvars
clear all
close all

format longG

tStart = cputime; % For CPU time determination %

%%% Length of spatial and temporal intervals %%%
ell = 1;
T   = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Please input numbers, e.g. 1, 2, ...');
problem = 4;
run(sprintf('%s%d.m','Test',problem));
% Exact_Test
% Test1
% Test2

%%% Division numbers of the temporal and the spatial intervals %%%
m = 6; % Division number of the spatial interval %
n = 2;
% n = T * (m / ell); % Division number of the temporal interval %
% n = T * (m / ell) * (m / ell); % Division number of the temporal interval %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Discretization of the temporal interval %%%
t   = linspace(0,T,n + 1);
tau = T / n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Discretization of the spatial interval %%%
x   = linspace(0,ell,m + 1);
h   = ell / m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eigcos = cos(pi / m); %%% Cosine appeared to eigenvalues %%%
delta  = h / tau; %% Fraction of h and tau needed for condition numbers %%

%%% The second-order derivative of u_0(x) %%%
d2psi0 = zeros(1,length(x));
d2cfd  = [-1, 16, -30, 16, -1]; %%% fourth-order accurate scheme %%%
for i = 2:m
  temp = linspace(x(i - 1),x(i + 1),5);
  % temp = x(i - 1):h/2:x(i + 1);
  d2psi0(i) = dot(d2cfd,psi0(temp)) / (3 * h * h);
endfor
clear i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initialization of approximation of the functions u(x,t) on the grid %%%
val_u      = zeros(length(t),length(x)); % exact values of u(x,t) %
app_u      = zeros(length(t),length(x)); % approximate values of u_k(x) %

q      = []; %%% q integral %%%
CondN2 = []; %%% Condition Number of Matrix %%%
%%% maximum error %%%
err = zeros(length(t),1);

%%% Tridiagonal system coefficients %%%
b   = zeros(1,length(x) - 2); % diagonal entries %
phi = zeros(1,length(x)); % right-hand side %
w   = zeros(1,length(x)); % Solution of tridiagonal system %

k = 1;

while k <= n + 1
  fprintf('Layer, k = %d\n', k);
  switch k
    case 1
      app_u(k,:) = psi0(x);
      tempdiff1 = FinDiff (h, app_u(k,:));
      tempdiff1 = tempdiff1 .* tempdiff1;
      %%% q integral %%%
      q1     = alpha(t(k)) + beta(t(k)) * SimpsRule (h, tempdiff1);
      condN2 = NaN; %%% Condition Number %%%
      q      = [q;q1];
      CondN2 = [CondN2;condN2]; %%% Condition Number %%%
      val_u(k,:) = u(x,t(k));
      err(k) = max(abs(app_u(k,:) - val_u(k,:)));
    case 2
      % app_u(k,:) = u(x,t(k));
      app_u(k,:) = app_u(k - 1,:) + tau * (psi1(x) + 0.5 * tau * ...
      (f(x,t(k - 1)) + q1 * d2psi0));
      tempdiff1 = FinDiff (h, app_u(k,:));
      tempdiff1 = tempdiff1 .* tempdiff1;
      %%% q integral %%%
      q2     = alpha(t(k)) + beta(t(k)) * SimpsRule (h, tempdiff1);
      condN2 = NaN; %%% Condition Number %%%
      q      = [q;q2];
      CondN2 = [CondN2;condN2]; %%% Condition Number %%%
      val_u(k,:) = u(x,t(k));
      err(k) = max(abs(app_u(k,:) - val_u(k,:)));
    otherwise
      v      = tau * tau * f(x,t(k - 1)) + 2 * app_u(k - 1,:);
      coeff1 = (delta * delta) / q2; % tridiagonal system coefficients %
      coeff2 = coeff1 / 6; % tridiagonal system coefficients %
      b(:)   = -2 * (1 + coeff1 * (1 + coeff2));
      for i = 2:m
        phi(i) = -coeff2 * (v(i - 1) + 2 * (5 + coeff1) * v(i) + v(i + 1));
      endfor
      clear i
      w(2:m) = tridiagonal_algotirhm (b, phi(2:m));
      app_u(k,:) = w - app_u(k - 2,:);
      val_u(k,:) = u(x,t(k));
      err(k) = max(abs(app_u(k,:) - val_u(k,:)));
      tempdiff1 = FinDiff (h, app_u(k,:));
      tempdiff1 = tempdiff1 .* tempdiff1;
      %%% q integral %%%
      q3     = alpha(t(k)) + beta(t(k)) * SimpsRule (h, tempdiff1);
      %%% Condition Number %%%
      condN2 = 1 + (2 * eigcos) / (1 + coeff1 * (1 + coeff2) - eigcos);
      %%%%%%%%%%%%%%%%%%%%%%%%
      q      = [q;q3];
      CondN2 = [CondN2;condN2]; %%% Condition Number %%%
      q1 = q2;
      q2 = q3;
  endswitch
  k = k + 1;
endwhile
clear ('k');



tEnd = cputime - tStart; % For CPU time determination %
if tEnd == 1
  fprintf('CPU time is %.2f second\n', tEnd);
elseif tEnd < 60
  fprintf('CPU time is %.2f seconds\n', tEnd);
else
  tEnd = tEnd / 60;
  switch tEnd
    case 1
      fprintf('CPU time is %.2f minute\n', tEnd);
    otherwise
      fprintf('CPU time is %.2f minutes\n', tEnd);
  endswitch
endif
clear tStart tEnd

dirname  = sprintf('Test%d',problem);
if not(isfolder(dirname))
  mkdir(dirname); % Folder creator %
endif

%%% CSV Files %%%
my_table = [x',app_u(end,:)',val_u(end,:)'];
filecsv  = sprintf('data_fig_Test%d_osc%d_n%d_m%d_all.csv',problem,osc,n,m);
cHeader  = {'x','app_u','val_u'};
save_as_csv (my_table, filecsv, dirname, cHeader);

my_table = [x',app_u(end,:)'];
filecsv  = sprintf('data_fig_Test%d_osc%d_n%d_m%d_appr.csv',problem,osc,n,m);
cHeader  = {'x','app_u'};
save_as_csv (my_table, filecsv, dirname, cHeader);

my_table = [(0:n)',CondN2];
filecsv  = sprintf('CondN2_Test%d_osc%d_n%d_m%d.csv',problem,osc,n,m);
cHeader  = {'k','CondN2'};
save_as_csv (my_table, filecsv, dirname, cHeader);

my_table = [t',err];
filecsv  = sprintf('error_Test%d_osc%d_n%d_m%d_appr.csv',problem,osc,n,m);
cHeader  = {'t','error'};
save_as_csv (my_table, filecsv, dirname, cHeader);

##node_t = [n / 4 + 1, n / 2 + 1, (3 * n) / 4 + 1, n + 1];
##node_t = round(node_t);
##clear("i");
##max_err    = [];
##max_CondN2 = [];
##for i = 1:length(node_t)
##  max_err    = [max_err;max(err(1:node_t(i)))];
##  max_CondN2 = [max_CondN2;max(CondN2(3:node_t(i)))];
##endfor
##my_table = [(node_t - 1)',t(node_t)',max_err,max_CondN2];
##filecsv  = sprintf('errorTable_Test%d_osc%d_n%d_m%d.csv',problem,osc,n,m);
##cHeader  = {'k','t_k',sprintf('error_n%d_m%d',n,m),'CondN2'};
##save_as_csv (my_table, filecsv, dirname, cHeader);
%%%%%%%%%%%%%%%%%

%%% generate mat file %%%
filemat  = sprintf('Results_Test%d_osc=%d_n=%d_m=%d.mat',problem,osc,n,m);
matfile  = fullfile(dirname,filemat);
##save(matfile,"t","tau","x","h","delta","CondN2","val_u","app_u","q","err",...
##"ell","T","m","n","node_t","problem");
save(matfile,"t","tau","x","h","delta","CondN2","val_u","app_u","q","err",...
"ell","T","m","n","problem");
%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except matfile
close all

load(matfile);

clearAllMemoizedCaches
clear matfile

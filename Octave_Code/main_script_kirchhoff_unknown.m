## Author: Zurab Vashakidze <zurab.vashakidze@gmail.com>
## Created: 2023-06-16

clc
clearAllMemoizedCaches
clearvars
clear all
close all

format longG

tStart = cputime; % For CPU time determination %

r     = 1; %%% Power of the tolerance, e.g. 10^(-r) %%%
Tol   = sprintf('1e-%d',r);
Tol   = str2double(Tol);

%%% Length of spatial and temporal intervals %%%
ell   = 4;
T     = 1;
Test  = 'Test4_(UnkSoln)';
run(sprintf('%s.m',Test)); %% Run the test %%

################################################################################
m      = 1024;
%%%  Discretization of the spatial interval %%%
x      = linspace(0,ell,m + 1);
h      = ell / m;

eigcos = cos(pi / m); %%% Cosine appeared to eigenvalues %%%

%%% The second-order derivative of u_0(x) %%%
d2psi0 = zeros(1,length(x));
d2cfd  = [-1, 16, -30, 16, -1]; %%% fourth-order accurate scheme %%%
for i = 2:m
  temp = linspace(x(i - 1),x(i + 1),5);
  d2psi0(i) = dot(d2cfd,psi0(temp)) / (3 * h * h);
endfor
clear ('i')

%% Tridiagonal system coefficients %%
b         = zeros(1,length(x) - 2); %% diagonal entries %%
phi       = zeros(1,length(x)); %% right-hand side %%
w         = zeros(1,length(x)); %% Solution of tridiagonal system %%


n     = 32; %%% Initial division of the temporal interavl %%%
u_old = zeros(n + 1,m + 1);

n     = 2*n; %%% Division improvement %%%
u_new = ones(n + 1,m + 1);
################################################################################
N     = [];
Delta = [];

count = 1;

%%% For condition number %%%
CondN2   = cell();
max_cond = []; %%% Max condition number for each discrete time step %%%
error    = []; %%% For errors %%%

while max(max(abs(u_old - u_new(1:2:end,:)))) >= Tol
  u_old = u_new;
  n     = 2*n;
  cond  = []; %%% For condition number %%%
  N     = [N;n];

  %% Discretization of the temporal interval %%
  t         = linspace(0,T,n + 1);
  tau       = T / n;

  delta = h / tau;
  Delta = [Delta;delta];

  %% Initialization of approximation of the functions u(x,t) on the grid %%
  u_new     = zeros(length(t),length(x)); %% approximate values of u_k(x) %%

  k = 1;

  while k <= n + 1
    switch k
      case 1
        u_new(k,:) = psi0(x);
        tempdiff1  = FinDiff (h, u_new(k,:));
        tempdiff1  = tempdiff1 .* tempdiff1;
        %% q integral %%
        q1         = alpha(t(k)) + beta(t(k)) * SimpsRule (h, tempdiff1);
        condn      = NaN; %%% Condition Number %%%
        cond       = [cond;condn]; %%% Condition Number %%%
      case 2
        u_new(k,:) = u_new(k - 1,:) + tau * (psi1(x) + 0.5 * tau * ...
        (f(x,t(k - 1)) + q1 * d2psi0));
        tempdiff1  = FinDiff (h, u_new(k,:));
        tempdiff1  = tempdiff1 .* tempdiff1;
        %% q integral %%
        q2         = alpha(t(k)) + beta(t(k)) * SimpsRule (h, tempdiff1);
      otherwise
        v      = tau * tau * f(x,t(k - 1)) + 2 * u_new(k - 1,:);
        coeff1 = (delta * delta) / q2; % tridiagonal system coefficients %
        coeff2 = coeff1 / 6; % tridiagonal system coefficients %
        b(:)   = -2 * (1 + coeff1 * (1 + coeff2));
        for i = 2:m
          phi(i) = -coeff2 * (v(i - 1) + 2 * (5 + coeff1) * v(i) + v(i + 1));
        endfor
        clear ('i')
        w(2:m) = tridiagonal_algotirhm (b, phi(2:m));
        u_new(k,:) = w - u_new(k - 2,:);
        tempdiff1 = FinDiff (h, u_new(k,:));
        tempdiff1 = tempdiff1 .* tempdiff1;
        %% q integral %%
        q3 = alpha(t(k)) + beta(t(k)) * SimpsRule (h, tempdiff1);
        %%% Condition Number %%%
        condn  = 1 + (2 * eigcos) / (1 + coeff1 * (1 + coeff2) - eigcos);
        cond   = [cond;condn]; %%% Condition Number %%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        q1 = q2;
        q2 = q3;
    endswitch
    k = k + 1;
  endwhile
  error    = [error;max(max(abs(u_old - u_new(1:2:end,:))))];
  cond     = [cond;NaN]; %%% Condition Number %%%
  max_cond = [max_cond;max(cond(2:(end - 1)))];
  CondN2   = [CondN2,cond];
  count    = count + 1;
  fprintf('n = %d\n', n);
  disp('------------------------');
endwhile
error(1) = NaN;
count = count - 1;

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

TestFold = sprintf('%s',Test);
TolFold  = sprintf('Tol_%d',Tol);
dirname  = fullfile(TestFold,TolFold);
if not(isfolder(dirname))
  mkdir(dirname); % Folder creator %
endif

node_t = [1, n / 4 + 1, n / 2 + 1, (3 * n) / 4 + 1, n + 1];
node_t = round(node_t);

my_table_all = []; %% for all values of one table %%

for i = 1:length(node_t)
  my_table     = [x',u_new(node_t(i),:)'];
  %%% Generate CSV file for SLR parameters %%%
  filecsv      = sprintf('%s_ell_%.2f_T_%.2f_soln_at_node_%d_Tol%f.csv',...
  Test,ell,T,node_t(i) - 1,Tol);
  cHeader      = {'x',sprintf('u%.2f',t(node_t(i)))};
  save_as_csv (my_table, filecsv, dirname, cHeader);

  %% all values in one table %%%
  my_table_all = [my_table_all,u_new(node_t(i),:)'];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
endfor
my_table_all = [x',my_table_all];

%%% Generate CSV file %%%
filecsv      = sprintf('%s_ell_%.2f_T_%.2f_soln_at_all_nodes_Tol%f.csv',...
Test,ell,T,Tol);
cHeader      = {'x'};
for i = 1:length(node_t)
  cHeader = [cHeader,sprintf('u%.2f',t(node_t(i)))];
endfor
clear ('i');
save_as_csv (my_table_all, filecsv, dirname, cHeader);

%%% Errors with tolerance and Maximum of Condition Number in CSV file %%%
my_table = [(0:1:count - 1)',error,N,max_cond];
filecsv  = sprintf(...
'%s_ell_%.2f_T_%.2f_error_and_max_cond_number_for_N_%d_to_%d_Tol%f.csv',...
Test,ell,T,N(1),N(end),Tol);
cHeader  = {'count','Error','N','MaxCondN2'};
save_as_csv (my_table, filecsv, dirname, cHeader);
################################################################################

my_table = [Tol,N(end),count - 1,error(end)];
filecsv  = sprintf(...
'%s_ell_%.2f_T_%.2f_error_Tol%f_with_loop_numb_%d.csv',...
Test,ell,T,Tol,count - 1);
cHeader  = {'Tol','max_N','count','error'};
save_as_csv (my_table, filecsv, dirname, cHeader);
################################################################################

%%% Condition number CSV file %%%
my_table = [(0:n)',CondN2{count}];
filecsv  = sprintf('%s_ell_%.2f_T_%.2f_condition_number_n%d_Tol%f.csv',...
Test,ell,T,n,Tol);
cHeader  = {'k','CondN2'};
save_as_csv (my_table, filecsv, dirname, cHeader);

################################################################################
app_u = u_new(node_t,:)'; %%% If you need all data, please use u_new %%%

%%% generate mat file %%%
filemat  = sprintf('Results_%s_ell_%.2f_T_%.2f_Tol%f.mat',Test,ell,T,Tol);
matfile  = fullfile(dirname,filemat);
save(matfile,'Tol','h','m','n','N','node_t','t','tau','app_u','x','T','ell',...
'Delta','CondN2','count','error');
%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except matfile
close all

load(matfile);

clearAllMemoizedCaches
clear matfile

clear all; clc; close all;
%{
This is a toy example of ASSD. 
%}

%% Simulation parameters
p = 2000; %  dimension 
n = 300;  % the number of samples 

wvar = 0.5; % the variance of noise
s0 = 30; % the number of non-zero entries

%% Generate data 
supp = randperm(p,s0);
rand_num = unifrnd(0.5,1,s0,1);
beta0 = zeros(p,1);
beta0(supp,1) = rand_num;

% Generate measurement matrix

% Correlated Gaussian matrix
pi = 0; % \pi=0 and 0.7 means no and strong correlations respectively
CORR = zeros(p,p);
for corr = 1:size(CORR,1)
    for corrj = 1:size(CORR,2)
        CORR(corr,corrj) = pi^(abs(corr-corrj));
    end
end
mu = zeros(1,p);
X = mvnrnd(mu,CORR,n);  % each row of the matrix X are drawn independently from N(0,\bm{\Sigma})

% % Structured matrix
% r = n+5;  % r= n + 3000 and r = n + 5 corresponding to weak correlated and highly correlated (or structured) matrix, respectively. 
% X1 = randn(n,r);
% X2 = randn(r,p);
% X = X1*X2;

% generate observation
z = X*beta0;
w = sqrt(wvar)*randn(n,1);
y = z + w;
        
noiselev = norm(w);
% normalize 
[nX,DE] = normalize(X);
     
%%  ASSD 
t1 = cputime;
% decimation with early stop
[beta_ssd,support] = ssd(nX, y, noiselev);
nonzero_find = beta_ssd(support,1); 

% second-stage thresholding
[~,d] = sort(abs(nonzero_find));  % sort the nonzero elements in ascending order 
beta_sort = nonzero_find(d);  
h_half = beta_sort(1:round(length(beta_sort)/2));
sigma_hat = std(h_half);     % compute $\hat{\sigma}$
R = 20;
tau = 0:0.01:R;
beta_old = beta_ssd;
Bic_ssd = [];
for i = 1:size(tau,2)
    theta = tau(i)*sigma_hat*sqrt(2*log(p)); 
    beta_ssd(abs(beta_ssd)<theta) = 0;
    residual_ssd = norm(nX*beta_ssd-y); 
    Bickk = 0.5*(residual_ssd^2) + log(n)*length(find(beta_ssd));  % compute the BIC
    Bic_ssd = [Bic_ssd,Bickk];
    beta_ssd = beta_old;
end
% output the coefficient vector $\bm{\beta}$ if the corresponding BIC value is the minimum.
id = find(Bic_ssd == min(Bic_ssd), 1);
theta = tau(id)*sigma_hat*sqrt(2*log(p)); 
beta_ssd(abs(beta_ssd)<theta) = 0;    
beta_ssd = DE*beta_ssd;
t2 = cputime - t1;
fprintf('ASSD done using time (second) = %g\n',t2)

%% TPR,FPR
[TP_ASSD, FP_ASSD] = TP_FP(beta0, beta_ssd);

%% relative error
rel2error_ASSD = norm(beta_ssd-beta0)/norm(beta0);


%% plot result
stem(beta0,'-x');
hold on
stem(beta_ssd,'--ro');
hold off
h = legend('TRUE','ASSD','Location','Best');
set(h,'FontName','Times New Roman','FontSize',8,'Box','off')
xlabel('location') 
ylabel('\beta')

%% output result
fprintf('TP = %g, FP = %g, RE = %5.4f, CPU time = %g\n', TP_ASSD, FP_ASSD, rel2error_ASSD, t2)


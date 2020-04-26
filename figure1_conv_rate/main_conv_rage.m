clc;
clear;

%% construct adjacency matrix K
m = 3;
n = 4;


index = 1;
for i = 1:m*n
    if (mod(i,n)>0)
        K(index,i) = -1;
        K(index,i+1) = 1;
        index = index +1;
    end
end

for i = 1:m*n-n
    K(index,i) = -1;
    K(index,i+n) = 1;
    index = index+1;
end

weight = 1;

[numedges, numvex] = size(K);
K = weight.*K;

%% construct preconditioner
I = eye(numedges);

partition1 = 1:9;
partition2 = 10:17;

K1 = K(partition1,:);
K2 = K(partition2,:);
T1 = K1*K1';
T2 = K2*K2';
P1 = I(partition1,:);
P2 = I(partition2,:);
T = P1'*T1*P1 + P2'*T2*P2;

rng(4);
p_generate = rand(numedges,1)*4;
f = K'*-p_generate-1e-8;

%% unpreconditioned proximal gradient
A = K';
eigs = svd(A'*A);
tau = 1/eigs(1);
p = zeros(numedges,1);
p_tilde = p;
max_it = 500;

stat.p = zeros(numedges,max_it+1);
stat.active = zeros(numedges,max_it+1);
stat.rate = zeros(1, max_it);

stat.p(:,1) = p_tilde;
stat.acitve(:,1) = 0;
primal_energy = zeros(1,max_it);
dual_energy = zeros(1,max_it);

for l = 1:max_it
    p = p - tau*K*(K'*p+f);
    
    p = min(max(p,-1), 1);
    
    primal_eng = 1/2*(K'*p+f)'*(K'*p+f);
    
    stat.active(:,l+1) = abs(abs(p)-ones(numedges,1))<1e-10;
    Pt = eye(numedges);
    inds = find(stat.active(:,l+1)==0);
    
    Pt(inds,:) = 0;
    Pt_m =  eye(size(T))- Pt*pinv(Pt);
    At = A*Pt_m;
    
    
    eigs = svd(At'*At);
    delta_max = eigs(1);
    delta_min = min(eigs(eigs>1e-5));
    stat.rate(1,l) = max(abs(1-tau*delta_max), abs(1-tau*delta_min));
    
    
    
    u = K'*p+f;
    stat.p(:,l+1) = p;
    dual_eng = 1/2*norm(u)^2 - u'*f + sum(abs(K*u));
    if(abs(dual_eng+primal_eng)<1e-12)
        break;
    end
end
p_opt = stat.p(:,l+1);
p_emp = sqrt(sum((stat.p-p_opt).^2,1));
p_emp_1 = p_emp(2:end)./p_emp(1:end-1);
rate = stat.rate(1,l);

%% preconditioned proximal gradient
A = K'*T^(-1/2);
eigs = svd(A'*A);
tau = 1/eigs(1);
p = zeros(numedges,1);
p_tilde = p;
max_it = 500;

stat.p = zeros(numedges,max_it+1);
stat.active = zeros(numedges,max_it+1);
stat.rate = zeros(1, max_it);

stat.p(:,1) = p_tilde;
stat.acitve(:,1) = 0;
primal_energy = zeros(1,max_it);
dual_energy = zeros(1,max_it);

for l = 1:max_it
    p_tilde = p_tilde - tau*T^(-1/2)*K*(K'*T^(-1/2)*p_tilde+f);
    
    
    p_tilde(partition1) = backward_solver(K1, p_tilde, numvex, T1, partition1, weight);
    p_tilde(partition2) = backward_solver(K2, p_tilde, numvex, T2, partition2, weight);
    
    
    p = T^(-1/2)*p_tilde;
    primal_eng = 1/2*(K'*p+f)'*(K'*p+f);
    
    stat.active(:,l+1) = abs(abs(p)-ones(numedges,1))<1e-10;
    Pt = eye(numedges);
    inds = find(stat.active(:,l+1)==0);
    
    Pt(inds,:) = 0;
    Pt_m =  eye(size(T))- T^(-1/2)*Pt*pinv(T^(-1/2)*Pt);
    At = A*Pt_m;
    
    
    eigs = svd(At'*At);
    delta_max = eigs(1);
    delta_min = min(eigs(eigs>1e-5));
    stat.rate(1,l) = max(abs(1-tau*delta_max), abs(1-tau*delta_min));
    
    
    
    u = K'*p+f;
    stat.p(:,l+1) = p;
    dual_eng = 1/2*norm(u)^2 - u'*f + sum(abs(K*u));
    if(abs(dual_eng+primal_eng)<1e-12)
        break;
    end
end
p_opt = stat.p(:,l+1);
p_emp_precond = sqrt(sum((stat.p-p_opt).^2,1));
%p_emp_1 = p_emp(2:end)./p_emp(1:end-1);
rate_precond = stat.rate(1,l);

%% unpreconditioned Hoffman bound
C = [eye(numedges)];
hoff_c = hoffman_constant(C, K');
miu = 1/hoff_c^2 / norm(K)^2;
hoff_rate = sqrt((1-miu)/(1+miu));

%% preconditioned Hoffman bound
C = [eye(numedges)*T^(-1/2)];
hoff_c_precond = hoffman_constant(C,A);
miu = 1/hoff_c_precond^2 *tau;
hoff_rate_precond = sqrt((1-miu)/(1+miu));

%% visiualization
n = 100;
% precond practical rate
semilogy(p_emp_precond(1:n),'LineWidth',1.5);
hold on;
x = 1:100;

% precond theoretical local rate
local_rate_precond = rate_precond;
local_y_precond = local_rate_precond.^x*3;
semilogy(local_y_precond, '--', 'LineWidth',1.5);

% precond hoffman rate
hoff_y_precond = hoff_rate_precond.^x*5;
semilogy(hoff_y_precond, '--', 'LineWidth', 1.5);

% unprecond practical rate
semilogy(p_emp(1:n), 'LineWidth', 1.5);

% unprecond theoretical local rate
local_rate = rate;
local_y = local_rate.^x*2;
semilogy(local_y, '--', 'LineWidth',1.5);

% unprecond hoffman rate
hoff_y = hoff_rate.^x*5;
semilogy(hoff_y, '--', 'LineWidth', 1.5);

% plot
legend('Empirical (Precond.)', 'Theoretical (Precond.)', 'Hoffman (Precond.)', ...
       'Empirical (Unprecond.)', 'Theoretical (Unprecond.)', ...
       'Hoffman (Unprecond.)', 'Location', 'southwest')

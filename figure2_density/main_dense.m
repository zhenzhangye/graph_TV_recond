clc;
clear;
addpath('../')

%% initialization
rng(5);
tol = 1e-10;

numV = 512;
A = rand(numV, numV);
A = A + A';
org_A = A;

f = rand(numV,1);

% threshold = 1.8;
% lambda = 2e-3:1e-3:1e-2;
% lambda = [lambda 1.1e-2:2e-3:1e-1];
% lambda = [lambda 1.1e-1:1e-1:0.7];
% threshold = 1.9;
% lambda = 0.39:0.027:0.92;
% lambda = [lambda 0.1:0.007:0.379];
% lambda = [lambda 0.003:0.0022:0.09];
threshold = 1.85;
lambda = 0.11:0.006:0.433;
lambda = [lambda 0.08:0.001:0.108];
lambda = [lambda 0.005:0.004:0.07];
iters.unprecond = zeros(length(threshold),length(lambda));
iters.nested = zeros(length(threshold),length(lambda));
iters.inactive = zeros(length(threshold),length(lambda));
run_time.unprecond = zeros(1, length(lambda));
run_time.nested = zeros(1, length(lambda));
run_time.inactive = zeros(1, length(lambda));
active_percent = zeros(1,length(lambda));

for i = 1:length(lambda)
    clear mex;
    max_it = 20000;
   
    A = sparse(double(org_A >threshold));
    A = A - diag(diag(A));
    A_2D = A*lambda(i);

    G = precond_mex('CreateGraphFromAdjacency', A_2D);

    [nabla, weights] = precond_mex('GraphToNabla', G);

    num_edges = size(nabla, 1);
    num_vertices = size(nabla, 2);
    W = spdiags(weights, 0, num_edges, num_edges);
    K = W*nabla;
    
    % unpreconditioned case
    u = zeros(num_vertices, 1);
    p = zeros(num_edges, 1);
    
    primal_energy = 0.5*sum((u-f).^2) + sum(abs(K*u));
    dual_energy = 0.5*sum((K'*p).^2)-f'*(K'*p);
    gap_zero = primal_energy+dual_energy;

    tic
    normK = normest(K,1e-6);
    tau = 1/normK^2;
    run_time.unprecond(i) = run_time.unprecond(i)+toc;
    
    for it = 1:max_it
        tic
        p = p - tau*K*(K'*p-f);
        p = max(min(p,1),-1);
        run_time.unprecond(i) = run_time.unprecond(i)+toc;
        
        u = f-K'*p;
        primal_energy = 0.5*sum((u-f).^2) + sum(abs(K*u));
        dual_energy = 0.5*sum((K'*p).^2)-f'*(K'*p);
        
        if(mod(it,10)==0)
            fprintf('%d: duality gap: %e\n', it, (primal_energy+dual_energy)/gap_zero);
        end
        if((primal_energy+dual_energy)/gap_zero<tol)
            break;
        end
    end
    
    iters.unprecond(i) = it;
    
    % nested_forest
    u = zeros(num_vertices, 1);
    a = ones(num_vertices,1);
    p = zeros(num_edges, 1);
    tic;
    forest_precond_nest = precond_mex('PartitionInitialGraph', G, true, false, zeros(num_edges, 1)); 
    L = precond_mex('GetNumberOfForests', forest_precond_nest);
    tau = L;
    run_time.nested(i) = run_time.nested(i)+toc;
    
     for it=1:max_it
        tic;
        p = precond_mex('ForestBackwardPDHG', forest_precond_nest, p, -K'*p + f, weights, tau, a);
        run_time.nested(i) = run_time.nested(i)+toc;
        
        u = f-K'*p;
        primal_energy = 0.5*sum((u-f).^2) + sum(abs(K*u));
        dual_energy = 0.5*sum((K'*p).^2)-f'*(K'*p);
        
        if(mod(it,10)==0)
            fprintf('%d: duality gap: %e\n', it, (primal_energy+dual_energy)/gap_zero);
        end
        if((primal_energy+dual_energy)/gap_zero<tol)
            break;
        end
     end
    
     iters.nested(i) = it;
     
     %inactive_forest
      u = zeros(num_vertices, 1);
    p = zeros(num_edges, 1);
    tic;
    forest_precond = precond_mex('PartitionInitialGraph', G, true, false, zeros(num_edges, 1)); 
    L = precond_mex('GetNumberOfForests', forest_precond);
    tau = L;
    run_time.inactive(i) = run_time.inactive(i)+toc;
    
    for it=1:max_it
        tic;
        p = precond_mex('ForestBackwardPDHG', forest_precond, p, -K'*p + f, weights, tau, a);
      
        forest_precond = precond_mex('PartitionGraphActiveSets', G, forest_precond, p, true);
        L = precond_mex('GetNumberOfForests', forest_precond);
        tau = L;
        run_time.inactive(i) = run_time.inactive(i)+toc;
        
        u = f-K'*p;
        primal_energy = 0.5*sum((u-f).^2) + sum(abs(K*u));
        dual_energy = 0.5*sum((K'*p).^2)-f'*(K'*p);
        
        if(mod(it,10)==0)
            fprintf('%d: duality gap: %e\n', it, (primal_energy+dual_energy)/gap_zero);
        end
        if((primal_energy+dual_energy)/gap_zero<tol)
            break;
        end
    end
    active_percent(i) = sum(abs(1-abs(p))<1e-13)/num_edges;
    iters.inactive(i) = it;
end
figure(1);
hold on;
[~,idx] = sort(active_percent);
plot(active_percent(idx), log10(iters.unprecond), 'k', 'LineWidth', 1.5);  
plot(active_percent(idx), log10(iters.nested), 'r', 'LineWidth', 1.5);
plot(active_percent(idx), log10(iters.inactive), 'b', 'LineWidth', 1.5);
axis([0 1 1 4.5])
legend('none', 'nested', 'inactive');
xlabel('$\frac{|\mathcal{A}|}{|\mathcal{E}|}$','Interpreter', 'latex', 'FontSize',20);
ylabel('number of iteration');
dense = num_edges/num_vertices;
clc;
clear;
addpath('../')

%% Initialization
max_it = 500;

% construct graph
rng(1);

nx = 100;
ny = 100;
f = rand(nx*ny,1);
tol = 1e-13;

A1 = spdiags(ones(nx,3), [-1 0 1], nx, nx);
A2 = spdiags(ones(ny,3), [-1 0 1], ny, ny);
A_2D = kron(A1, speye(ny)) + kron(speye(nx), A2);

A_2D = A_2D + A_2D';

G = precond_mex('CreateGraphFromAdjacency', A_2D);

nabla = precond_mex('GraphToNabla', G);

num_edges = size(nabla, 1);
num_vertices = size(nabla, 2);
weights = rand(num_edges,1)*0.55;
W = spdiags(weights, 0, num_edges, num_edges);
nabla = W*nabla;
a = ones(num_vertices,1);


%% unpreconditioned proximal gradient

tau = 1/normest(nabla)^2;
p = zeros(num_edges,1);
u = f;
primal0= 0.5 * sum( (u-f).^2 ) + sum(abs(nabla * u));
dual0 = 0.5 * sum((nabla' * p).^2) - f' * (nabla' * p);
gap0 = primal0 + dual0;

primal_eng = zeros(1,max_it);
dual_eng = zeros(1,max_it);
dual_gap = zeros(1,max_it);

for l = 1:max_it
    p = p - tau*nabla*(nabla'*p-f);
    
    p = min(max(p,-1), 1);
    
    u = f-nabla'*p;
    primal_eng(l) = 0.5 * sum( (u-f).^2 ) + sum(abs(nabla * u));
    dual_eng(l) = 0.5 * sum((nabla' * p).^2) - f' * (nabla' * p);
    dual_gap(l) = abs(primal_eng(l) + dual_eng(l))/gap0;
       
    if(dual_gap(l)<tol)
        break;
    end
end
non_precond_gap = [1 dual_gap];
non_precond_it = l;

%% Chains precondition
chain_precond = precond_mex('PartitionInitialGraph', G, true, true, zeros(num_edges, 1));
L = precond_mex('GetNumberOfForests', chain_precond);
tau = L;

p = zeros(num_edges,1);
primal_eng = zeros(1,max_it);
dual_eng = zeros(1,max_it);
dual_gap = zeros(1,max_it);


for l = 1:max_it
    p = precond_mex('ForestBackwardPDHG', chain_precond, p, -nabla'*p + f, weights, tau, a);
    
    u = f-nabla'*p;
    primal_eng(l) = 0.5 * sum( (u-f).^2 ) + sum(abs(nabla * u));
    dual_eng(l) = 0.5 * sum((nabla' * p).^2) - f' * (nabla' * p);
    dual_gap(l) = abs(primal_eng(l) + dual_eng(l))/gap0;
       
    if(dual_gap(l)<tol)
        break;
    end
    
    if (mod(l,50)==0)
        fprintf('iteration: %d, gap: %.4e, active edges: %.2f\n', l, dual_gap(l), sum(abs(1-abs(p))<1e-10)/length(p))
    end
       
end
chain_precond_gap = [1 dual_gap];
chain_precond_it = l;

%% Nested forest precondition
nested_precond = precond_mex('PartitionInitialGraph', G, true, false, zeros(num_edges, 1)); 
L = precond_mex('GetNumberOfForests', nested_precond);
tau = L;

p = zeros(num_edges,1);
primal_eng = zeros(1,max_it);
dual_eng = zeros(1,max_it);
dual_gap = zeros(1,max_it);

for l = 1:max_it
    p = precond_mex('ForestBackwardPDHG', nested_precond, p, -nabla'*p + f, weights, tau, a);
    
    u = f-nabla'*p;
    primal_eng(l) = 0.5 * sum( (u-f).^2 ) + sum(abs(nabla * u));
    dual_eng(l) = 0.5 * sum((nabla' * p).^2) - f' * (nabla' * p);
    dual_gap(l) = abs(primal_eng(l) + dual_eng(l))/gap0;
       
    if(dual_gap(l)<tol)
        break;
    end
    
    if (mod(l,50)==0)
        fprintf('iteration: %d, gap: %.4e, active edges: %.2f\n', l, abs(dual_gap(l)), sum(abs(1-abs(p))<1e-10)/length(p))
    end

end
nested_precond_gap = [1 dual_gap];
nested_precond_it = l;

%% inactive forest precond n = 20
inactive_precond_n20 = precond_mex('PartitionInitialGraph', G, true, false, zeros(num_edges, 1)); 
L = precond_mex('GetNumberOfForests', inactive_precond_n20);
tau = L;

p = zeros(num_edges,1);
primal_eng = zeros(1,max_it);
dual_eng = zeros(1,max_it);
dual_gap = zeros(1,max_it);

for l = 1:max_it
    p = precond_mex('ForestBackwardPDHG', inactive_precond_n20, p, -nabla'*p + f, weights, tau, a);
    
    u = f-nabla'*p;
    primal_eng(l) = 0.5 * sum( (u-f).^2 ) + sum(abs(nabla * u));
    dual_eng(l) = 0.5 * sum((nabla' * p).^2) - f' * (nabla' * p);
    dual_gap(l) = abs(primal_eng(l) + dual_eng(l))/gap0;
       
    if(dual_gap(l)<tol)
        break;
    end
    
    if (mod(l,50)==0)
        fprintf('iteration: %d, gap: %.4e, active edges: %.2f\n', l, abs(dual_gap(l)), sum(abs(1-abs(p))<1e-10)/length(p))
    end
    
    if(mod(l,20)==0)
        inactive_precond_n20 = precond_mex('PartitionGraphActiveSets', G, inactive_precond_n20, p, true);
        L = precond_mex('GetNumberOfForests', inactive_precond_n20);
        tau = L;
    end 
end
inactive_precond_n20_gap = [1 dual_gap];
inactive_precond_n20_it = l;

%% inactive forest precond n = 10
inactive_precond_n10 = precond_mex('PartitionInitialGraph', G, true, false, zeros(num_edges, 1)); 
L = precond_mex('GetNumberOfForests', inactive_precond_n10);
tau = L;

p = zeros(num_edges,1);
primal_eng = zeros(1,max_it);
dual_eng = zeros(1,max_it);
dual_gap = zeros(1,max_it);

for l = 1:max_it
    p = precond_mex('ForestBackwardPDHG', inactive_precond_n10, p, -nabla'*p + f, weights, tau, a);
    
    u = f-nabla'*p;
    primal_eng(l) = 0.5 * sum( (u-f).^2 ) + sum(abs(nabla * u));
    dual_eng(l) = 0.5 * sum((nabla' * p).^2) - f' * (nabla' * p);
    dual_gap(l) = abs(primal_eng(l) + dual_eng(l))/gap0;
       
    if(dual_gap(l)<tol)
        break;
    end
    
    if (mod(l,50)==0)
        fprintf('iteration: %d, gap: %.4e, active edges: %.2f\n', l, dual_gap(l), sum(abs(1-abs(p))<1e-10)/length(p))
    end
    
    if(mod(l,10)==0)
        inactive_precond_n10 = precond_mex('PartitionGraphActiveSets', G, inactive_precond_n10, p, true);
        L = precond_mex('GetNumberOfForests', inactive_precond_n10);
        tau = L;
    end 
end
inactive_precond_n10_gap = [1 dual_gap];
inactive_precond_n10_it = l;

%% inactive forest precond n = 5
inactive_precond_n5 = precond_mex('PartitionInitialGraph', G, true, false, zeros(num_edges, 1)); 
L = precond_mex('GetNumberOfForests', inactive_precond_n5);
tau = L;

p = zeros(num_edges,1);
primal_eng = zeros(1,max_it);
dual_eng = zeros(1,max_it);
dual_gap = zeros(1,max_it);

for l = 1:max_it
    p = precond_mex('ForestBackwardPDHG', inactive_precond_n5, p, -nabla'*p + f, weights, tau, a);
    
    u = f-nabla'*p;
    primal_eng(l) = 0.5 * sum( (u-f).^2 ) + sum(abs(nabla * u));
    dual_eng(l) = 0.5 * sum((nabla' * p).^2) - f' * (nabla' * p);
    dual_gap(l) = abs(primal_eng(l) + dual_eng(l))/gap0;
       
    if(dual_gap(l)<tol)
        break;
    end
    
    if (mod(l,50)==0)
        fprintf('iteration: %d, gap: %.4e, active edges: %.2f\n', l, dual_gap(l), sum(abs(1-abs(p))<1e-10)/length(p))
    end
    
    if(mod(l,5)==0)
        inactive_precond_n5 = precond_mex('PartitionGraphActiveSets', G, inactive_precond_n5, p, true);
        L = precond_mex('GetNumberOfForests', inactive_precond_n5);
        tau = L;
    end 
end
inactive_precond_n5_gap = [1 dual_gap];
inactive_precond_n5_it = l;

%% inactive forest precond n = 1
inactive_precond_n1 = precond_mex('PartitionInitialGraph', G, true, false, zeros(num_edges, 1)); 
L = precond_mex('GetNumberOfForests', inactive_precond_n1);
tau = L;

p = zeros(num_edges,1);
primal_eng = zeros(1,max_it);
dual_eng = zeros(1,max_it);
dual_gap = zeros(1,max_it);

for l = 1:max_it
    p = precond_mex('ForestBackwardPDHG', inactive_precond_n1, p, -nabla'*p + f, weights, tau, a);
    
    u = f-nabla'*p;
    primal_eng(l) = 0.5 * sum( (u-f).^2 ) + sum(abs(nabla * u));
    dual_eng(l) = 0.5 * sum((nabla' * p).^2) - f' * (nabla' * p);
    dual_gap(l) = abs(primal_eng(l) + dual_eng(l))/gap0;
       
    if(dual_gap(l)<tol)
        break;
    end
    
    if (mod(l,50)==0)
        fprintf('iteration: %d, gap: %.4e, active edges: %.2f\n', l, dual_gap(l), sum(abs(1-abs(p))<1e-10)/length(p))
    end
    
    if(mod(l,1)==0)
        inactive_precond_n1 = precond_mex('PartitionGraphActiveSets', G, inactive_precond_n1, p, true);
        L = precond_mex('GetNumberOfForests', inactive_precond_n1);
        tau = L;
    end 
end
inactive_precond_n1_gap = [1 dual_gap];
inactive_precond_n1_it = l;

%% plot figure
figure(1)
total = 250;
start = 1;
plot(log10(non_precond_gap(start:non_precond_it)))
hold on;
plot(log10(chain_precond_gap(start:chain_precond_it)))
%hold on;

plot(log10((nested_precond_gap(start:nested_precond_it))))
h20 = plot(log10(inactive_precond_n20_gap(start:inactive_precond_n20_it)));
points_x = start:20:inactive_precond_n20_it;
points_y = h20.YData(points_x);
plot(points_x, points_y, '.', 'Color', h20.Color, 'MarkerSize',20)

h10 = plot(log10(inactive_precond_n10_gap(start:inactive_precond_n10_it)));
points_x = start:10:inactive_precond_n10_it;
points_y = h10.YData(points_x);
plot(points_x, points_y, '.', 'Color', h10.Color, 'MarkerSize',20)

h5 = plot(log10(inactive_precond_n5_gap(start:inactive_precond_n5_it)));
points_x = start:5:inactive_precond_n5_it;
points_y = h5.YData(points_x);
plot(points_x, points_y, '.', 'Color', h5.Color, 'MarkerSize',20)

h1 = plot(log10(inactive_precond_n1_gap(start:inactive_precond_n1_it)));
points_x = start:1:inactive_precond_n1_it;
points_y = h1.YData(points_x);
%plot(points_x, points_y, '.', 'Color', h1.Color, 'MarkerSize',20)
xlim([0. 250])
legend
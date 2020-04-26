%% 
clc;
clear all;
close all;
addpath('../');

%% Generate blurred image
img = im2double(rgb2gray(imresize(imread('img.jpg'), 0.1)));

[m,n,c] = size(img);

scale = 0.1;
kernel = fspecial('motion', 7, 45);
%kernel = fspecial('gaussian', 10);

B = convmtx2(kernel, m,n);
B = kron(speye(c), B);

kx = size(kernel, 1);
ky = size(kernel, 2);

m2 = m+kx-1;
n2 = n+ky-1;

f = B*img(:) + 0.01*randn(m2*n2*c, 1);

%% Generate linear operator

A1 = spdiags(ones(m,3), [-1 0 1], m, m);
A2 = spdiags(ones(n,3), [-1 0 1], n, n);
A_2D = kron(A2, speye(m)) + kron(speye(n), A1);

A_2D = A_2D + A_2D';

A_2D = A_2D/2;

A_2D = kron(speye(c), A_2D);

G  = precond_mex('CreateGraphFromAdjacency', A_2D * 0.1);

[nabla, weights] = precond_mex('GraphToNabla', G);

num_edges = size(nabla,1);

num_vertices = size(nabla,2);

W = spdiags(weights, 0, size(weights, 1), size(weights, 1));
nabla_w = W * nabla;

%% parameters
max_iters = 2500;
primal_energy_unprecond = zeros(max_iters,1);
primal_energy_chain = zeros(max_iters,1);
primal_energy_nested = zeros(max_iters,1);
primal_energy_inactive = zeros(max_iters,1);
primal_energy_diag = zeros(max_iters,1);

time_unprecond = zeros(max_iters,1);
time_chain = zeros(max_iters,1);
time_nested = zeros(max_iters,1);
time_inactive = zeros(max_iters,1);
time_diag = zeros(max_iters,1);



%% unpreconditioned
if true
fprintf('unprecond...\n');
normK = normest_precise(nabla_w,1e-6);
s = normK * scale;
t = normK / scale;
A = B'*B + s*speye(num_vertices);
u = zeros(num_vertices, 1);
p = zeros(num_edges, 1);
primal_energy_unprecond(1) = 1/2*norm(B*u-f)^2 + sum(sum(abs(nabla_w*u)));

for it = 1:max_iters*3
        % primal update
        tic;
        u_prev = u;
        
        [u,~,~,~,~] = pcg(A,  (B'*f - nabla_w' * p + s*u), 1e-12, 1000, [], [], u_prev);
    
        % over-relaxation
        u_bar = 2*u - u_prev;
    
        % dual update
        p = min(max(p + 1/t*nabla_w*u_bar,-1) , 1);        
        time_unprecond(it + 1) = time_unprecond(it) + toc;
        primal_energy_unprecond(it + 1) = 1/2*norm(B*u-f)^2 + sum(sum(abs(nabla_w*u)));

end
end

%% inactively nested forest
fprintf('inactive...\n');
u = zeros(num_vertices, 1);
a = ones(num_vertices,1);
p = zeros(num_edges, 1);
primal_energy_inactive(1) = 1/2*norm(B*u-f)^2 + sum(sum(abs(nabla_w*u)));
%forest_precond = precond_mex('create_chains', G);
forest_precond = precond_mex('PartitionInitialGraph', G, true, false, zeros(num_edges, 1)); 
%forest_precond = precond_mex('create_nested_forests', G,  'dfs', false);
%precond_mex('initialize_forest', forest_precond, true);
L = precond_mex('GetNumberOfForests', forest_precond);
%precond_mex('initialize_forest', forest_precond, true);
%[L, sigmaM, sigmam] = precond_mex('get_number_of_forests', forest_precond);
s = scale;
t = L / scale;
A = B'*B;
for it = 1:max_iters
        % primal update
        tic;
        u_prev = u;
        [u,~,~,~,~] = pcg(A+s*speye(num_vertices),(B'*f - nabla_w' * p + s*u), 1e-12, 1000, [], [], u_prev);
    
        % over-relaxation
        u_bar = 2*u - u_prev;
        
        % dual update
        p = precond_mex('ForestBackwardPDHG', forest_precond, p, u_bar, weights, t, a);

        %p = precond_mex('forest_backward_pdhg', forest_precond, p, u_bar, weights, t, false);
        
        % reconditioning
        if(mod(it, 5)==0)
            %forest_precond = precond_mex('partition_graph_active_sets', G, p, forest_precond);
            %precond_mex('initialize_forest', forest_precond, true);
            %forest_precond = precond_mex('partition_initial_graph_active_sets', G, p, true, forest_precond);
            %[L, sigmaM, sigmam] = precond_mex('get_number_of_forests', forest_precond);
            
            forest_precond = precond_mex('PartitionGraphActiveSets', G, forest_precond, p, true);
            L = precond_mex('GetNumberOfForests', forest_precond);
            
            s = scale;
            t = L / scale;
        end
        time_inactive(it + 1) = time_inactive(it) + toc;
        
        primal_energy_inactive(it+1) = 1/2*norm(B*u-f)^2 + sum(sum(abs(nabla_w*u)));
        
end


%% nested forest
if true
fprintf('nested...\n');
u = zeros(num_vertices, 1);
p = zeros(num_edges, 1);
primal_energy_nested(1) = 1/2*norm(B*u-f)^2 + sum(sum(abs(nabla_w*u)));
forest_precond = precond_mex('PartitionInitialGraph', G, true, false, zeros(num_edges, 1)); 
%forest_precond = precond_mex('create_nested_forests', G, 'dfs', false);
%precond_mex('initialize_forest', forest_precond, true);
L = precond_mex('GetNumberOfForests', forest_precond);
s = scale;
t = L / scale;
A = B'*B + s*speye(num_vertices);
for it = 1:max_iters*2
        % primal update
        tic;
        u_prev = u;
        [u,~,~,~,~] = pcg(A,  (B'*f - nabla_w' * p + s*u), 1e-12, 1000, [], [], u_prev);
    
        % over-relaxation
        u_bar = 2*u - u_prev;
        
        %dual update
        p = precond_mex('ForestBackwardPDHG', forest_precond, p, u_bar, weights, t, a);

        %p = precond_mex('forest_backward_pdhg', forest_precond, p, u_bar, weights, t, false);    
                
        time_nested(it + 1) = time_nested(it) + toc;

        primal_energy_nested(it+1) = 1/2*norm(B*u-f)^2 + sum(sum(abs(nabla_w*u)));
        
end
end

%% diagonal preconditioner
if true
fprintf('diag...\n');
u = zeros(num_vertices, 1);
p = zeros(num_edges, 1);
primal_energy_diag(1) = 1/2*norm(B*u-f)^2 + sum(sum(abs(nabla_w*u)));
s = scale * max(sum(abs(nabla_w),1));
t = max(sum(abs(nabla_w),2)) / scale;
A = B'*B + s*speye(num_vertices);

for it = 1:max_iters*3
        tic;
        % primal update
        u_prev = u;
        
        [u,~,~,~,~] = pcg(A,  (B'*f - nabla_w' * p + s*u), 1e-12, 1000, [], [], u_prev);
    
        % over-relaxation
        u_bar = 2*u - u_prev;
    
        % dual update
        p = min(max(p + (1/t)*nabla_w*u_bar,-1) , 1);        
            
        time_diag(it + 1) = time_diag(it) + toc;
        primal_energy_diag(it+1) = 1/2*norm(B*u-f)^2 + sum(sum(abs(nabla_w*u)));

end
end

%% chains
fprintf('chains...\n');
u = zeros(num_vertices, 1);
p = zeros(num_edges, 1);
primal_energy_chain(1) = 1/2*norm(B*u-f)^2 + sum(sum(abs(nabla_w*u)));
forest_precond = precond_mex('PartitionInitialGraph', G, true, true, zeros(num_edges, 1));
%forest_precond = precond_mex('create_chains', G);
%precond_mex('initialize_forest', forest_precond, true);
L = precond_mex('GetNumberOfForests', forest_precond);
s = scale;
t = L / scale;
A = B'*B + s*speye(num_vertices);
for it = 1:max_iters*2
    tic;
        % primal update
        u_prev = u;
        [u,~,~,~,~] = pcg(A,  (B'*f - nabla_w' * p + s*u), 1e-12, 1000, [], [], u_prev);
    
        % over-relaxation
        u_bar = 2*u - u_prev;
        
        %dual update
        p = precond_mex('ForestBackwardPDHG', forest_precond, p, u_bar, weights, t, a);

        %p = precond_mex('forest_backward_pdhg', forest_precond, p, u_bar, weights, t, false);    
        
        time_chain(it + 1) = time_chain(it) + toc;
        primal_energy_chain(it+1) = 1/2*norm(B*u-f)^2 + sum(sum(abs(nabla_w*u)));
        
end


save results
return
%% iter vs energy
opti = min(min(primal_energy_inactive),min(primal_energy_chain))-1e-12;
%opti = min(primal_energy_inactive);

for i=1:(max_iters+1)
    primal_energy_inactive2(i) = min(primal_energy_inactive(1:i));
end

figure(1);
hold on;
plot(log10(primal_energy_unprecond-opti),'LineWidth', 3);
plot(log10(primal_energy_diag-opti),'LineWidth', 3);
plot(log10(abs(primal_energy_inactive2-opti)),'LineWidth', 3);
plot(log10(abs(primal_energy_nested-opti)),'LineWidth', 3);
plot(log10(abs(primal_energy_chain-opti)),'LineWidth', 3);
legend('Unpreconditioned','Diagonally preconditioned','Inactively Nested Forests','Nested Forests','Chains');
%legend('Inactive', 'Chains');
xlim([0, 1100]);


%% time vs energy
figure;
hold on;
plot(time_unprecond, log10(primal_energy_unprecond-opti),'LineWidth', 3);
plot(time_diag, log10(primal_energy_diag-opti),'LineWidth', 3);
plot(time_inactive, log10(abs(primal_energy_inactive2-opti)),'LineWidth', 3);
plot(time_nested, log10(abs(primal_energy_nested-opti)),'LineWidth', 3);
plot(time_chain, log10(abs(primal_energy_chain-opti)),'LineWidth', 3);
legend('Unpreconditioned','Diagonally preconditioned','Inactively Nested Forests','Nested Forests','Chains');
%legend('Inactive', 'Chains');
xlim([0, 40]);

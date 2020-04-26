function [u_star, time, iters] = pdhg_deblur(G, f, B, nabla, weights, alg, forest_option, optimal_eng, tol)
iters = 0;
time = 0;

max_iters = 2000;

primal_energy = zeros(1, max_iters);

num_edges = size(nabla, 1);
num_vertices = size(nabla, 2);
W = spdiags(weights, 0, size(weights, 1), size(weights, 1));
nabla_w = W * nabla;

u = zeros(num_vertices,1);
p = zeros(num_edges,1);

%% unpreconditioned PDHG
if(strcmp(alg, 'pdhg'))
    tic;
    normK = normest_precise(nabla_w,1e-6);
    s = 0.11;
    t = 10*normK^2;
    time = time+toc;
    
    A = B'*B + t*speye(num_vertices);
    
    for it = 1:max_iters
        
        tic;
        
        % primal update
        u_prev = u;
        
        [u,~,~,~,~] = pcg(A,  (B'*f - nabla_w' * p + s*u), 1e-12, 1000, [], [], u_prev);
    
        % over-relaxation
        u_bar = 2*u - u_prev;
    
        % dual update
        p = min(max(p + 1/t*nabla_w*u_bar,-1) , 1);        
        
        time = time+toc;
    
        primal_energy(it) = 1/2*norm(B*u-f)^2 + sum(sum(abs(nabla_w*u)));
    
        fprintf('%d-th iteration, optimal gap: %e\n', it, primal_energy(it)-optimal_eng);
        
        if(primal_energy(it)-optimal_eng<tol)
            break;
        end
    end
end


%% inactive forests
if(strcmp(alg, 'inactive_forests'))
    tic;
    forest_precond = precond_mex('create_chains', G);
    %forest_precond = precond_mex('create_nested_forests', G, 'dfs', true);
    precond_mex('initialize_forest', forest_precond, true);
    [L, sigmaM, sigmam] = precond_mex('get_number_of_forests', forest_precond);
    s = 0.11;
    t = 10*L;
    time = time+toc;
    A = B'*B + t*speye(num_vertices);
    for it = 1:max_iters
        tic;
        % primal update
        u_prev = u;
        [u,~,~,~,~] = pcg(A,  (B'*f - nabla_w' * p + s*u), 1e-12, 1000, [], [], u_prev);
    
        % over-relaxation
        u_bar = 2*u - u_prev;
        
        % dual update
        p = precond_mex('forest_backward_pdhg', forest_precond, p, u_bar, weights, t, false);
        
        % reconditioning
        if(mod(it,40)==0)
            forest_precond = precond_mex('partition_graph_active_sets', G, p, forest_precond);
            precond_mex('initialize_forest', forest_precond, true);
            %forest_precond = precond_mex('partition_initial_graph_active_sets', G, p, true, forest_precond);
            [L, sigmaM, sigmam] = precond_mex('get_number_of_forests', forest_precond);
            s = 0.11;
            t = 10*L;
        end
        
        time = time+toc;

        primal_energy(it) = 1/2*norm(B*u-f)^2 + sum(sum(abs(nabla_w*u)));
        
        if(mod(it,1) == 0)
            fprintf('%d-th iteration, optimal gap: %e\n', it, primal_energy(it)-optimal_eng);
        end
        
        if(primal_energy(it)-optimal_eng<tol)
            break;
        end
    end
end


%% forests
if(strcmp(alg, 'forests'))
    tic;
    if(forest_option)
        forest_precond = precond_mex('create_chains',G);
    else
        forest_precond = precond_mex('create_nested_forests', G, 'dfs', forest_option);
    end
    precond_mex('initialize_forest', forest_precond, true);
    [L, sigmaM, sigmam] = precond_mex('get_number_of_forests', forest_precond);
    s = 0.5*sqrt(L);
    t = 2.1*sqrt(L);
    time = time+toc;
    A = B'*B + t*speye(num_vertices);
     for it = 1:max_iters
        tic;
        % primal update
        u_prev = u;
        [u,~,~,~,~] = pcg(A,  (B'*f - nabla_w' * p + s*u), 1e-12, 1000, [], [], u_prev);
    
        % over-relaxation
        u_bar = 2*u - u_prev;
        
        %dual update
        p = precond_mex('forest_backward_pdhg', forest_precond, p, u_bar, weights, t, false);    
        
        time = time+toc;

        primal_energy(it) = 1/2*norm(B*u-f)^2 + sum(sum(abs(nabla_w*u)));
        
        if(mod(it,1) == 0)
            fprintf('%d-th iteration, optimal gap: %e\n', it, primal_energy(it)-optimal_eng);
        end
        
        if(primal_energy(it)-optimal_eng<tol)
            break;
        end
    end
end

%% return
u_star = u;
iters = it;
end
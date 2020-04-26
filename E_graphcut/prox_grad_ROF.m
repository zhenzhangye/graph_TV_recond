function [u_star, inactive, iters, total_time, num_forests, made_chains, sigmam, gap] = prox_grad_ROF(G, unaries, precond, tol, force_chain)
%% Initial settings
u_star = unaries;
inactive = 0;
iters = 0;
sigmam = 0;
gap = 0;

q = 2;

max_iters = 500000;
num_forests = -1;
made_chains = -1;

[nabla, weights] = precond_mex('GraphToNabla', G);
thresh = 1e-5;
f = -thresh - unaries;

%scale = max(max(abs(f)), max(weights));
scale = 1;
weights = weights / scale; % scale weights to [0, 1]
f = f / scale; % scale f to [-1, 1]

num_edges = size(nabla, 1);
num_vertices = size(nabla, 2);
W = spdiags(weights, 0, size(weights, 1), size(weights, 1));
nabla_w = W * nabla;

% initialize
total_time = 0;
u = zeros(num_vertices, 1);
p = zeros(num_edges, 1);
p_prev = p;
en_primal = 0.5 * sum( (u-f).^2 ) + sum(abs(nabla_w * u));
en_dual = 0.5 * sum((nabla_w' * p).^2) - f' * (nabla_w' * p);
gap_zero = en_primal + en_dual;

%% none precondition
if strcmp(precond, 'none')
    tic;
    norm_K = normest_precise(nabla_w, 1e-6);
    tau = 1/norm_K^2;
    t = 1;
    total_time = toc;
    
    for it=1:max_iters
        % check gap
        en_primal = 0.5 * sum( (u-f).^2 ) + sum(abs(nabla_w * u));
        en_dual = 0.5 * sum((nabla_w' * p_prev).^2) - f' * (nabla_w' * p_prev);
        gap = (en_primal + en_dual) / gap_zero;
        
        if(mod(it,50)==0)
            fprintf('none: %d-th iteration, duality gap:%e\n', it, gap);
        end
        
        if gap < tol
            break
        end
        
        tic;
        
        %gradient step
        p = p - tau*(nabla_w*(nabla_w'*p-f));
        
        %proximal step
        p = min(max(p,-1),1);
        
        %t_prev = t;
        %t = 0.5*(1+ sqrt(1+gamma*t^2));
        t = (it-1)/(it+q);
        xk = p;
        %p = p+(t_prev-1)*(p-p_prev)/t;
        p = p + t*(p-p_prev);
        p_prev = xk;
                
        total_time = total_time + toc;
        
        %recover primal variable
        u = f-nabla_w'*p_prev;
    end
end

%% forests
if strcmp(precond, 'forest')
     try
            if force_chain == true
                forest_precond = precond_mex('PartitionInitialGraph', G, true, true, zeros(num_edges, 1));
                made_chains = 1;
                tic;
            else
                tic;
                forest_precond = precond_mex('PartitionInitialGraph', G, true, false, zeros(num_edges, 1)); 
                made_chains = 0; 
            end
        catch
            tic;
            forest_precond = precond_mex('PartitionInitialGraphForceChain', G);
            made_chains = 0;
        end

        L = precond_mex('GetNumberOfForests', forest_precond);
        num_forests = L;
        tau = L;
        p_prev = p;
        total_time = toc;
        a = ones(num_vertices,1); % the weight in front of unary. Assume 1 here.
        
     for it=1:max_iters
        % check gap
        en_primal = 0.5 * sum( (u-f).^2 ) + sum(abs(nabla_w * u));
        en_dual = 0.5 * sum((nabla_w' * p_prev).^2) - f' * (nabla_w' * p_prev);
        gap = (en_primal + en_dual) / gap_zero;
        
        if(mod(it,50)==0)
            fprintf('forest: %d-th iteration, duality gap:%e, time: %.3f\n', it, gap, total_time);
        end
            
        if gap < tol
            break
        end
        
        tic;
        
        %gradient step
        
        p = precond_mex('ForestBackwardPDHG', forest_precond, p, -nabla_w'*p + f, weights, tau, a);

        
        t = (it-1)/(it+q);
        xk = p;
        p = p + t * (p-p_prev);
        p_prev = xk;
        
        total_time = total_time + toc;
        
        %recover primal variable
        u = f-nabla_w'*p_prev;
    end

end

%% inactive forests
if strcmp(precond, 'inactive_forest')
    tic;
    forest_precond = precond_mex('PartitionInitialGraph', G,  true, false, zeros(num_edges,1));

    L = precond_mex('GetNumberOfForests', forest_precond);
    tau = L;
    p_prev = p;
    total_time = toc;
    a = ones(num_vertices,1); % the weight in front of unary. Assume 1 here.

    
    for it=1:max_iters
        % check gap
        en_primal = 0.5 * sum( (u-f).^2 ) + sum(abs(nabla_w * u));
        en_dual = 0.5 * sum((nabla_w' * p_prev).^2) - f' * (nabla_w' * p_prev);
        gap = (en_primal + en_dual) / gap_zero;
        
  
        if(mod(it,50)==0)
            fprintf('inactive: %d-th iteration, duality gap:%e, time: %.3f\n', it, gap, total_time);
        end
        
        if gap < tol
            break
        end
        
        tic;
       
        
        p = precond_mex('ForestBackwardPDHG', forest_precond, p, -nabla_w'*p + f, weights, tau, a);
        
        %
        if(mod(it,30)==0)
            forest_precond = precond_mex('PartitionGraphActiveSets', G, forest_precond, p, true);
            L = precond_mex('GetNumberOfForests', forest_precond);
            tau = L;
            p_prev = p;
        else
           t = (it-1)/(it+q);
           xk = p;
           p = p + t * (p-p_prev);
           p_prev = xk;
       end
        
        total_time = total_time + toc;
        
        %recover primal variable
        u = f-nabla_w'*p_prev;
    end
end

%% block diag
if strcmp(precond, 'block_diag')
    tic;
    diag_T = diag(nabla_w*nabla_w');
    diag_T = 1./sqrt(diag_T);
    diag_T(diag_T==0) = 1;
    p_diag_T = 1./diag_T;
    T = spdiags(diag_T, 0, num_edges, num_edges);
    norm_K = normest_precise(nabla_w'*T, 1e-6);
    tau = 1/norm_K^2;
    total_time = toc;
    
    p_tilde = p;
    
    for it=1:max_iters
        % check gap
        en_primal = 0.5 * sum( (u-f).^2 ) + sum(abs(nabla_w * u));
        en_dual = 0.5 * sum((nabla_w' * p).^2) - f' * (nabla_w' * p);
        gap = (en_primal + en_dual) / gap_zero;
        
        if(mod(it,50)==0)
            fprintf('block_diag: %d-th iteration, duality gap:%e\n', it, gap);
        end
        
        if gap < tol
            break
        end
        
        tic;
        
        %gradient step
        p_tilde = p_tilde - tau*(T*nabla_w*(nabla_w'*T*p_tilde-f));
        
        %proximal step
        p_tilde(p_tilde > p_diag_T) = p_diag_T(p_tilde > p_diag_T);
        p_tilde(p_tilde < -p_diag_T) = -p_diag_T(p_tilde< -p_diag_T);
        
        %t_prev = t;
        %t = 0.5*(1+ sqrt(1+gamma*t^2));
        t = (it-1)/(it+q);
        xk = p_tilde;
        p_tilde = p_tilde+t*(p_tilde-p_prev);
        p_prev = xk;
        
        
        total_time = total_time + toc;
        
        %recover primal variable
        p = p_prev.*diag_T;
        u = f-nabla_w'*p;
    end
end


u_star = double(scale * u > thresh);
iters = it;
inactive = 1 - sum(abs(1-abs(p))<1e-10)/num_edges;
end

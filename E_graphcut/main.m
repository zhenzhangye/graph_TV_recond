%% 
close all;
clear all;

%% read inputs
tol = 1e-10;

addpath('..');

files = {    
    '../datasets/rmf-long.n2.3142.shuf.bk', ...
    '../datasets/rmf-wide.n2.3141.shuf.bk', ...
    '../datasets/horse-48112-0037-153.bk', ...
    '../datasets/alue7065-33338-0001-042.bk', ...
    '../datasets/lux-106214-0003-019.bk', ...
    '../datasets/punch-us18-09-p.bk', ...
    '../datasets/BVZ-venus1.bk', ...
    '../datasets/lazybrush-mangagirl.max.bk', ...
    '../datasets/KZ2-sawtooth1.bk', ...
    '../datasets/ferro.bk', ...
};

c = clock;
filename = sprintf('results_%d_%02d_%02d_%02d_%02d.txt', c(1), c(2), c(3), c(4), c(5));
fileID = fopen(filename, 'w');
fprintf(fileID, '%e %s %d it=%d q=%d\n', tol, 30, 2);

%% Optimization part
for i = 1:size(files, 2)
    clear mex;
    fprintf('%s\n', files{i});
		if ~isfile(files{i})
			printf([files{i},' not found\n']);
			continue;
		end
    
    [G, unaries, en_bias] = precond_mex('CreateGraphFrombkFile', files{i});
    
    [nabla,weights] = precond_mex('GraphToNabla', G);
    numV = size(nabla, 2);
    numE = size(nabla, 1);
    
    avg_weights = weights/max(abs(unaries));
    avg_weights = sum(avg_weights)/length(weights);
    fprintf(fileID, '%s (n/1024=%.1f, n/m=%f, avg_weights=%.8f):\n', files{i}, numV/1024, numE/numV, avg_weights);
 
    [u_star_1, inactive_1, iters_1, time_1, ~, ~, ~, gap_1] = prox_grad_ROF(G, unaries, 'none', tol);
    cut_1 = en_bias + unaries' * u_star_1 + sum(weights .* abs(nabla * u_star_1));
    fprintf(fileID, ' - none                                              : iter=%6d, time=%.2fs, cut=%.1f, ratio=%.5f, gap=%e\n', ...
        iters_1, time_1, cut_1, inactive_1, gap_1);
 
    [u_star_2, inactive_2, iters_2, time_2, nf, ch, ~, gap_2] = prox_grad_ROF(G, unaries, 'forest', tol, false);
    cut_2 = en_bias + unaries' * u_star_2 + sum(weights .* abs(nabla * u_star_2));
    fprintf(fileID, ' - forest (nf=%3d, ch=%d)                            : iter=%6d, time=%.2fs, cut=%.1f, ratio=%.5f, gap=%e\n', ...
        nf, ch, iters_2, time_2, cut_2, inactive_2, gap_2);
 
    [u_star_3, inactive_3, iters_3, time_3, nf, ch, ~, gap_3] = prox_grad_ROF(G, unaries, 'inactive_forest', tol);
    cut_3 = en_bias + unaries' * u_star_3 + sum(weights .* abs(nabla * u_star_3));
    fprintf(fileID, ' - inactive_forest_FISTA (nf=%3d, ch=%d)             : iter=%6d, time=%.2fs, cut=%.1f, ratio=%.5f, gap=%e\n', ...
        nf, ch, iters_3, time_3, cut_3, inactive_3, gap_3);

    [u_star_4, inactive_4, iters_4, time_4, nf, ch, ~, gap_4] = prox_grad_ROF(G, unaries, 'block_diag', tol);
    cut_4 = en_bias + unaries' * u_star_4 + sum(weights .* abs(nabla * u_star_4));
    fprintf(fileID, ' - diagonal (nf=%3d, ch=%d)                          : iter=%6d, time=%.2fs, cut=%.1f, ratio=%.5f, gap=%e\n',...
        nf, ch, iters_4, time_4, cut_4, inactive_4, gap_4);
  
    [u_star_5, inactive_5, iters_5, time_5, nf, ch, ~, gap_5] = prox_grad_ROF(G, unaries, 'forest', tol,  true);
    cut_5 = en_bias + unaries' * u_star_5 + sum(weights .* abs(nabla * u_star_5));
    fprintf(fileID, ' - linear_forest (nf=%3d, ch=%d)                     : iter=%6d, time=%.2fs, cut=%.1f, ratio=%.5f, gap=%e\n\n', ...
        nf, ch, iters_5, time_5, cut_5, inactive_5, gap_5);
  
end

fclose(fileID);

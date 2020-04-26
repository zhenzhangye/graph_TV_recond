%% compute the hoffman constant with exhaustive search
%  Ax = b, Cx <= d
function result = hoffman_constant(C,A)
    [p,n] = size(A);
    [m,n] = size(C);
    r = rank([A',C']);
    rank_A = rank(A);
    
    result = -1;
    rows_A = nchoosek(1:p, rank_A);
    for i = 1:length(rows_A(:,1))
        A_I = A(rows_A(i,:), :);
        rows_C = nchoosek(1:m, r-rank(A_I));
        for j = 1:length(rows_C(:,1))
            C_I = C(rows_C(j,:), :);
            if(rank([A_I', C_I'])<r)
                continue;
            end
            eigs = svd([A_I', C_I']');
            result = max(result, 1/eigs(r));
        end
    end  
end
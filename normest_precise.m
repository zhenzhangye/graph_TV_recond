function [e,cnt] = normest_precise(S,tol)
%NORMEST_PRECISE Estimate the matrix 2-norm.
%   NORMEST_PRECISE(S) is an estimate of the 2-norm of the matrix S.
%   NORMEST_PRECISE(S,tol) uses relative error tol instead of 1.e-6.
%   [nrm,cnt] = NORMEST_PRECISE(..) also gives the number of iterations used.
%
%   This function is intended primarily for sparse matrices,
%   although it works correctly and may be useful for large, full
%   matrices as well.  Use NORMEST when your problem is large
%   enough that NORM takes too long to compute and an approximate
%   norm is acceptable.
%
%   See also NORM, COND, CONDEST.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.9 $  $Date: 1997/11/21 23:38:38 $

if nargin < 2, tol = 1.e-6; end
x = sum(abs(S))';
cnt = 0;
e = norm(x);
if e == 0, return, end
x = x/e;
e0 = 0;
while abs(e-e0) > tol*e
   e0 = e;
   Sx = S*x;
   if nnz(Sx) == 0
      Sx = rand(size(x));
   end
   e = norm(Sx);
   x = S'*Sx;
   x = x/norm(x);
   cnt = cnt+1;
end

end

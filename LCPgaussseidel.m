function [z,err] = LCPgaussseidel(M,q,tol,z0)
% syntax: [z,err] = lemke(M,q,tol,z0)
% LEMKE    Solves linear complementarity problems (LCPs).
% An LCP solves
%   Mz+q >=0, z>=0, z'(Mz+q)=0.
% The input z0 defines a starting basis; it can be either
% an initial guess of the solution or a vector of zeros and ones 
% with ones representing those z(i) thought to be non-zero in the
% solution.  For example, passing z=[1.5;0;2.2] has the same 
% effect as passing z=[1;0;1]. 
% If z0 is omitted the origin is used as a starting basis.
% ERR returns an error condition:
%   0: Solution found
%   1: Maximum iterations exceeded

% If NARGOUT==1, a warning message is displayed instead.
%
% ALGORITHM
%   Uses a modified Gauss-Seidel method for solution of LCPs

n = length(q);

maxiter = min([1000 25*n]);
err=0;


z = z0;

% Trivial solution exists
if all(q >= 0.)
    z = zeros(n,1);
    return;
end



for i = 1:maxiter
    for j = 1:n
        z(j) = max(0,z(j)-(q(j)+sum(M(j,:)*z))/M(j,j));
    end
    
    residue = max(abs(min(z,q+M*z)));

    if(residue < tol)
        return;
    end
end

err = 1;
% Display warning messages if no error code is returned
if nargout<2 
  Display('Warning: solution not found - Iterations exceeded limit');
  
end

end
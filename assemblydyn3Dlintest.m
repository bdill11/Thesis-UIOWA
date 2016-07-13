function [A,M,b] = assemblydyn3Dlintest(A,M,b,pde,elt,p,t,fht,intmethod)
% function [A,M,b] = assemblydyn2dlin(A,M,b,pde,elt,p,t,fht,intmethod)
% Adds the assembled matrices representing the
% given PDE (pde) to the A matrix & M matrix & b vector.
% This uses a given element (elt) with the triangulation given by (p,t).
% The feature hash table (fht) is used to obtain variable indexes
% for given features.  This is obtained by create_fht().

% This is for construction in dynamic 2d linear elasticity problems.
%
% A must be nv x nv and M must be nv x nv and b must be nv x 1 where nv is the total
% number of variables (as returned by fht_num_vars()).

% Reference triangle has vertices (0,0), (1,0), (0,1).
[p_int,w_int] = intmethod(); % points and weights for reference triangle

% np is the total number of points in the triangulation
np = size(p,1);

% compute nv = total number of variables
nv = fht_num_vars(fht);

% nv_elt is the number of variables in one element
nv_elt = sum(elt.nvars);

% order is the order of derivatives used in the assembly;
% we need 0 <= order <= 2
order = pde.order;

intval1 = zeros(nv_elt,nv_elt);
intval2 = zeros(nv_elt,nv_elt);
intval3 = zeros(nv_elt,size(b,2));

% Save get_Aphihat() values for all the integration points
% on the reference element
Aphihatvals = cell(length(w_int),1);
Mphihatvals = cell(length(w_int),1);
for k = 1:length(w_int)
    Aphihatvals{k} = elt.Aphihat(p_int(k,:),order);
end

for i = 1:size(t,1) % for all triangles ...
    % obtain variable list and signs for this triangle
    [vlist,slist] = get_var_triangle(t(i,:),fht,elt,np)
    % set up affine transformation xhat :-> x = T.xhat + b
    i1 = t(i,1);  i2 = t(i,2);  i3 = t(i,3);
    T = [p(i2,:)'-p(i1,:)', p(i3,:)'-p(i1,:)'];
    b0 = p(i1,:)';
    % form weighted sum of integrand at integration points
    intval1 = 0;
    intval2 = 0;
    intval3 = 0;
    for k = 1:length(w_int)
        Aphival      = elt.trans_Aphihat(T,Aphihatvals{k},order)
        %Mphival      = elt.trans_Aphihat(T,Mphihatvals{k},order-1);
        Dmat         = pde.coeffs(T*p_int(k,:)'+b0)
        Mmat         = pde.rhs1(T*p_int(k,:)'+b0);
        rhsvec       = pde.rhs2(T*p_int(k,:)'+b0);
        integrand_val1 = Aphival*Dmat*Aphival'
        integrand_val2 = Aphival*Mmat*Aphival';
        integrand_val3 = Aphival*rhsvec;
        intval1 = intval1 + w_int(k)*integrand_val1;
        intval2 = intval2 + w_int(k)*integrand_val2;
        intval3 = intval3 + w_int(k)*integrand_val3;
        pause
    end
    detT = abs(det(T));
    intval1 = intval1*detT;                % scale by Jacobian
    intval2 = intval2*detT;
    intval3 = intval3*detT;
    intval1 = diag(slist)*intval1*diag(slist);  % change signs if needed
    intval2 = diag(slist)*intval2*diag(slist);
    intval3 = bsxfun(@times,intval3,slist');
    A(vlist,vlist) = A(vlist,vlist) + intval1    % add to matrix & vec
    M(vlist,vlist) = M(vlist,vlist) + intval2;
    b(vlist,:) = b(vlist,:) + intval3;
    pause
end
function [A,b,bvlist] = assembly2dbdrycontact(pde,elt,p,t,bedges,tidx,fht,norms,intmethod)
% function [A,b,bvlist] = assembly2dbdry(pde,elt,p,t,bedges,tidx,fht,intmethod)
% Adds the assembled matrix and vector representing the
% given PDE (pde) to the A matrix & b vector.
% This uses a given element (elt) with the triangulation given by (p,t,bedges,tidx)
% for the boundary. Note that bedges(i,:) is in triangle t(tidx(i),:).
% The feature hash table (fht) is used to obtain variable indexes
% for given features. This is obtained by create_fht().
%
% A must be nv x nv and b must be nv x 1 where nv is the total
% number of variables (as returned by fht_num_vars()).
% Reference edge has vertices 0 and 1.
[p_int,w_int] = intmethod(); % points and weights for reference triangle
% np is the total number of points in the triangulation
np = size(p,1);
% compute nv = total number of variables
nv = fht_num_vars(fht);
% nv_edge is the number of variables in one edge (and associated points)
%nv_elt = sum(elt.nvars);
nv_edge = 0;
for i = 1:size(elt.flist,1)
    if sum(elt.flist(i,:) ~= 0) <= 2
        nv_edge = nv_edge + elt.nvars(i);
    end
end % for
% order is the order of differentiation used in the "PDE"
order = pde.order;
A = sparse(nv,nv);
b = zeros(nv,1);
bvlist = [];
for i = 1:size(bedges,1) % for all boundary edges ...
    % obtain variable list and signs for this triangle & boundary edge
    bedge = bedges(i,:);
    triangle = t(tidx(i),:);
    [tvlist,slist] = get_var_triangle(t(tidx(i),:),fht,elt,np);
    bvlist1 = get_var_edge(bedges(i,:),fht,np);
    bvlist  = [bvlist,bvlist1];
    match = match_edge_triangle(bedges(i,:),t(tidx(i),:));
    % set up affine transformation xhat :-> x = T.xhat + b0
    i1 = t(tidx(i),1); i2 = t(tidx(i),2); i3 = t(tidx(i),3);
    T = [p(i2,:)'-p(i1,:)', p(i3,:)'-p(i1,:)'];
    b0 = p(i1,:)';
    % Turn p_int on the interval [0,1] to points on the appropriate
    % edge of the reference triangle
    p_ref = [0 0; 1 0; 0 1];
    p_ref0 = p_ref(match(1),:);
    p_ref1 = p_ref(match(2),:);
    % form weighted sum of integrand at integration points
    intval1 = zeros(length(tvlist),length(tvlist));
    intval2 = zeros(length(tvlist),1);
    for k = 1:length(w_int)
        p_int_ref = (1-p_int(k))*p_ref0+p_int(k)*p_ref1;
        % p_int_val = T*p_int_ref'+b0;
        Aphivalhat = elt.Aphihat(p_int_ref,order);
        Aphival = elt.trans_Aphihat(T,Aphivalhat,order);
        Dmat = pde.coeffs(T*p_int_ref'+b0,norms(i,:));
        rhsvec = pde.rhs(T*p_int_ref'+b0);
        integrand_val1 = Aphival*Dmat*Aphival';
        integrand_val2 = Aphival*rhsvec;
        intval1 = intval1 + w_int(k)*integrand_val1;
        intval2 = intval2 + w_int(k)*integrand_val2;
    end
    detT = norm(p(t(tidx(i),match(1)),:)-p(t(tidx(i),match(2)),:),2);
    intval1 = intval1*detT;
    intval2 = intval2*detT;
    intval1 = diag(slist)*intval1*diag(slist); % change signs if needed
    intval2 = slist'.*intval2;
    
    A(tvlist,tvlist) = A(tvlist,tvlist) + intval1; % add to matrix & vec
    %A(triangle,triangle) = A(triangle,triangle) + intval1;
    b(tvlist)        = b(tvlist) + intval2;
end
bvlist = unique(sort(bvlist));
v_array = ones(nv,1);
v_array(bvlist) = 0;
cbvlist = find(v_array ~= 0);
Ab(cbvlist,:) = 0;
Ab(:,cbvlist) = 0;
b(cbvlist) = 0;

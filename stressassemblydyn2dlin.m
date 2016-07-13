function [sig11,sig12,sig21,sig22] = stressassemblydyn2dlin(ux,uy,elt,p,t,fht,lamb,mu)


% Reference triangle has vertices (0,0), (1,0), (0,1).

% np is the total number of points in the triangulation
np = size(p,1);

% compute nv = total number of variables
%nv = fht_num_vars(fht);

% nv_elt is the number of variables in one element
%nv_elt = sum(elt.nvars);

axmat = zeros(np);
bxmat = zeros(np);
aymat = zeros(np);
bymat = zeros(np);

index = zeros(np);

for i = 1:size(t,1) % for all triangles ...
    % obtain variable list and signs for this triangle
    %[vlist,slist] = get_var_triangle(t(i,:),fht,elt,np);
    % set up affine transformation xhat :-> x = T.xhat + b
    i1 = t(i,1);  i2 = t(i,2);  i3 = t(i,3);
    T = [p(i2,:)'-p(i1,:)', p(i3,:)'-p(i1,:)'];
    b0 = p(i1,:)';
  
    index(t(i,1),i) = i;
    index(t(i,2),i) = i;
    index(t(i,3),i) = i;
    
    ax = 1/det(T)*(T(2,1)*(ux(i1)-ux(i3))+T(2,2)*(ux(i2)-ux(i1)));
    bx = 1/det(T)*(T(1,2)*(ux(i1)-ux(i2))+T(1,1)*(ux(i3)-ux(i1)));
    
    ay = 1/det(T)*(T(2,1)*(uy(i1)-uy(i3))+T(2,2)*(uy(i2)-uy(i1)));
    by = 1/det(T)*(T(1,2)*(uy(i1)-uy(i2))+T(1,1)*(uy(i3)-uy(i1)));
    
    axmat(t(i,:),i) = ax;
    bxmat(t(i,:),i) = bx;
    
    aymat(t(i,:),i) = ay;
    bymat(t(i,:),i) = by;
  
end

du1dx1 = zeros(np,1);
du1dx2 = zeros(np,1);
du2dx1 = zeros(np,1);
du2dx2 = zeros(np,1);

sig = zeros(2,2,np);

for i=1:np
    indices = find(index(i,:)); %Find triangles where each node is actually involved
    
    %Now average all the derivatives for all the triangles that a given
    %node is part of
    du1dx1(i,1) = mean(axmat(i,indices));
    du1dx2(i,1) = mean(bxmat(i,indices));
    
    du2dx1(i,1) = mean(aymat(i,indices));
    du2dx2(i,1) = mean(bymat(i,indices));
    
    %And build the stress tensor for each node
    sig(:,:,i) = [(lamb+2*mu)*du1dx1(i) + lamb*du2dx1(i),mu*(du1dx2(i)+du2dx1(i));mu*(du1dx2(i)+du2dx1(i)),(lamb+2*mu)*du2dx2(i) + lamb*du1dx1(i)];
end

eigs = eig(sig(:,:,1));
abs(eigs(1)-eigs(2))






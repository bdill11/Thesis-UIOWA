function [integral] = bodyintegrator3d(z,p,t,intmethod)
%Given function values z at each node of the triangulation described by p
%and t, this returns the integral of the z function (treated as linear 
%between nodes)over the entire region described by the triangulation.


% Reference triangle has vertices (0,0), (1,0), (0,1).

% np is the total number of points in the triangulation
np = size(p,1);

[p_int,w_int] = intmethod(); % points and weights for reference triangle

zvals = zeros(length(w_int),1);

integral = 0;
for i = 1:size(t,1) % for all triangles ...
    % set up affine transformation xhat :-> x = T.xhat + b
    i1 = t(i,1);  i2 = t(i,2);  i3 = t(i,3);
    T = [p(i2,:)'-p(i1,:)', p(i3,:)'-p(i1,:)'];
    b0 = p(i1,:)';
    
    %compute z values at each integration point 
    for k = 1:length(w_int)
        zvals(k) = ((b0(1)+T(1,1))*z(i2)-b0(1)*z(i1))*p_int(k,1)^2+((b0(1)+T(1,2))*z(i3)-b0(1)*z(i1))*p_int(k,2)+b0(1)*z(i1);
    end
  
    detT = abs(det(T));
    integral = integral + 2*pi*detT*w_int'*zvals;
  
end
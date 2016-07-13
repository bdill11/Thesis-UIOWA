function [n] = edgenormal(p,t,bedges,tidx)

%This function computes a list of outward normals for the edges input.  We
%compute a perpendicular to the edge and then correct the sign by checking
%the dot product of the perpendicular with a vector from one of the
%associated nodes to the interior node of the associated triangle.  If this
%dot product is positive, we need to reverse the sign

n = zeros(size(bedges,1),2);

for i = 1:size(bedges,1)
   
    v = [p(bedges(i,2),1)-p(bedges(i,1),1);p(bedges(i,2),2)-p(bedges(i,1),2)];
    ni = [-v(2);v(1)];
    ni = ni/norm(ni);
    
    
    for j = 1:3
       if(t(tidx(i),j) ~= bedges(i,1) && t(tidx(i),j) ~= bedges(i,2))
           intpointidx = t(tidx(i),j);
       end
    end
    
    w = [p(intpointidx,1)-p(bedges(i,1),1);p(intpointidx,2)-p(bedges(i,1),2)];
    
    if(dot(ni,w) > 0)
        ni = -ni;
    end
    
    n(i,:)=ni;
end
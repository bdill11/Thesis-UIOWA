function [gaps,nodereturn] = gap3(p,bnodes,bedges,normals,height,fht)

%This function computes the gap function for a given set of boundary nodes,
%boundary edges and associated outward normal vectors.
%We first compute the outward normal for each Bnode by averaging the
%normals for the Bedges touching it (at most 2).  Then we use this normal
%to compute how far the Bnode is from the foundation in the normal
%direction, returning this value for each bnode.  The foundation is taken
%to be at y = height




    




np = size(p,1);
pvlist = get_pvlist(fht,np);
nodenorms = zeros(length(bnodes),2);
gapf = zeros(length(pvlist),1);
tempnodenorms = zeros(length(pvlist),2);
gaps = [];
orderednodenorms = [];

for i=1:length(bnodes)
    
    bnodenorms = [];
    for j=1:length(bedges(:,1))
       if(bedges(j,1) == bnodes(i) || bedges(j,2) == bnodes(i))
           bnodenorms = [bnodenorms;normals(j,:)];
       end
    end
    nodenorms(i,:) = mean(bnodenorms,1)./norm(mean(bnodenorms,1));
    
    if(nodenorms(i,2) < -.01)
        gapf(pvlist(bnodes(i))) = (height - p(bnodes(i),2))/nodenorms(i,2);
        
        tempnodenorms(pvlist(bnodes(i)),:) = nodenorms(i,:);
        
    else
        gapf(i) = Inf;
    end
    
end

for i=1:length(pvlist)
     if(gapf(i) ~= 0) 
         gaps = [gaps; gapf(i)];
         orderednodenorms = [orderednodenorms; tempnodenorms(i,:)];
     end
end

nodereturn = sparse(tempnodenorms);
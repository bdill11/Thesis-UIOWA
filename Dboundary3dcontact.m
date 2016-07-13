function [cbedges,cbnodes,ct_index,fbedges,fbnodes,ft_index] = Dboundary3dcontact(t,p,tol)
%Pass in a condition which will identify a portion of boundary
%TEST FIRST PART, IT GIVES SOME EDGES THAT ARE NOT BEDGES

% function [bedges,bnodes,t_index] = boundary2d(t)
%
% Construct boundary edge list from triangle list t
% t is ntriangles x 3, bedges = nedges x 2
% Edge k joins points bd(k,1) and bd(k,2).
%
% Simply check when edges only appear once in the triangle list.
%
% Also returns the triangle index for each boundary edge
t = sort(t,2); % sort each row of t
bd1 = sortrows([t(:,1),t(:,2),(1:size(t,1))';
                t(:,2),t(:,3),(1:size(t,1))';
                t(:,1),t(:,3),(1:size(t,1))']);
[bd2,idx1] = unique(bd1(:,1:2),'rows','first');
[bd2,idx2] = unique(bd1(:,1:2),'rows','last');
eqlist = find(idx1 == idx2);
bedges = bd1(idx1(eqlist),1:2);
t_index = bd1(idx1(eqlist),3);
bnodes = unique(sort(bedges(:)));

%%%%%%%%%%%
%c stands for conditional, these are the bedges that meet a specified
%condition.  f stands for failed, these are bedges that fail to meet the
%condition
cbedges = [];
ct_index = [];
cbnodes = [];

fbedges = [];
ft_index = [];
fbnodes = [];

% 
% for i=1:length(bnodes)-1
%     if(p(bnodes(i),2) > 0.7 )
%         cbedges = [cbedges;bedges(i+1,:)];
%         ct_index = [ct_index;t_index(i+1)]; 
%         cbnodes = [cbnodes;bnodes(i)];
%     else
%         fbedges = [fbedges;bedges(i,:)];
%         ft_index = [ft_index;t_index(i)]; 
%         fbnodes = [fbnodes;bnodes(i)];
%     end 
% end

%Find bnodes subject to a condition

for i=1:length(bnodes)
    if(p(bnodes(i),1) < tol/10)
        if(p(bnodes(i),2) == -.01)
            fbnodes = [fbnodes;bnodes(i)];
        end
        cbnodes = [cbnodes;bnodes(i)];
    elseif(p(bnodes(i),2) < -0.0095)
        fbnodes = [fbnodes;bnodes(i)];
    end 
end

%Find bedges and t_index subject to a condition (same condition as above
%generally)

for i=1:length(bnodes)
    if(any(bedges(i,1) == cbnodes) && any(bedges(i,2) == cbnodes))
        cbedges = [cbedges;bedges(i,:)];
        ct_index = [ct_index;t_index(i)]; 
    elseif(any(bedges(i,1) == fbnodes) && any(bedges(i,2) == fbnodes))
        fbedges = [fbedges;bedges(i,:)];
        ft_index = [ft_index;t_index(i)]; 
    end 
end

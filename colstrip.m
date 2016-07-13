function [B] = colstrip(A)

%This function removes any columns which consist entiredly of 0 from a
%matrix A
B=sparse(0,0);
[row,col]=find(A);
col = unique(col);
sumA = sum(A);
for i=1:size(col,1)
    if(abs(sumA(col(i)))>1e-9)
        B = [B, A(:,col(i))];
    end
end
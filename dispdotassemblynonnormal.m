function [dispdot] = dispdotassemblynonnormal(Nbnodenorms, matsize)

j=1;
rows = [];
cols = [];
vals = [];
for i=1:size(Nbnodenorms,1)
    if((Nbnodenorms(i,1)*Nbnodenorms(i,2))~=0)
        rows = [rows, j,j];
        cols = [cols, 2*i-1,2*i];
        vals = [vals, Nbnodenorms(i,1), Nbnodenorms(i,2)];
        j = j+1;
    end
end

dispdot = sparse(rows, cols, vals, matsize(1),matsize(2));

end
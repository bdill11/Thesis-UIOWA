function [] = surfy(N, gap, numt, numx,dispfrac)


[t,x] = size(N);
t = t/dispfrac;
stept = floor(t/numt);

[X,Y] = meshgrid(1:numx,1:stept:t);


figure
surf(X,Y,N(1:stept:t,1:numx));
savefig(strcat('/space/bdill/Desktop/Code 6-14-16/Animations/Normals',int2str(1),'.fig'));

figure
surf(X,Y,gap(1:stept:t,1:numx));
savefig(strcat('/space/bdill/Desktop/Code 6-14-16/Animations/Gaps',int2str(1),'.fig'));





end
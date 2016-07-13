nu = .3;
mu = 10000;
R1 = 1;
R2 = 1;
A = (2-2*nu)/mu;
k = 1/R1 + 1/R2;
a = .01;
p0 = 4*k*a/(pi*A);
z=0:a/10:3*a;
ra = r/a;

tresca = zeros(length(z),1);
sigtest = zeros(length(z),1);
vonMise = zeros(length(z),1);

for i = 1:length(z)
    sigthth = -p0*((1+nu)*(1-z(i)/a*atan(a/z(i)))-1/(2*(1+z(i)^2/a^2)));
    sigzz = -p0/(1+z(i)^2/a^2);
    sigrr = sigthth;
   
    
    sig = [(sigrr),0,0;0,(sigthth),0;0,0,(sigzz)];
    
    eigens = eig(sig);
    
    %sigtest(i) = sigthth;
    sigtest(i) = abs(eigens(3)/2-eigens(1)/2);
    tresca(i) = max([abs(eigens(1)-eigens(2)), abs(eigens(1)-eigens(3)),abs(eigens(3)-eigens(2))]);
    vonMise(i) = sqrt(1/2*(((eigens(1)-eigens(2))^2 + (eigens(1)-eigens(3))^2 + (eigens(3)-eigens(2))^2)));
end

plot(z/a,sigtest/p0)
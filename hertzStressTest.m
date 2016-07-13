nu = .3;
mu = 10000;
R1 = 1;
R2 = 1;
A = (2-2*nu)/mu;
k = 1/R1 + 1/R2;
a = .01;
p0 = -4*k*a/(pi*A);
[r,z] = meshgrid(-1.0005*a:.001:a,.001:.0005:2*a);
[len,wid] = size(r);
ra = r/a;
za = z/a;

% P = 8*k*a^3/(3*A)

maxshear = zeros(len*wid,1);
sigtest = zeros(len*wid,1);

for i = 1:len*wid
    u = (.5*(ra(i)^2+za(i)^2-1+sqrt((ra(i)^2+za(i)^2-1)^2+4*za(i)^2)));
    sigthth = -p0*((1-2*nu)/3*1/ra(i)^2*(1-(za(i)/sqrt(u))^3)+za(i)/sqrt(u)*(2*nu+(1-nu)*u/(1+u)-(1+nu)*sqrt(u)*atan(1/sqrt(u))));
    sigzz = -p0*((za(i)/sqrt(u))^3*u/(u^2+za(i)^2));
    sigrr = p0*((1-2*nu)/3*1/ra(i)^2*(1-(za(i)/sqrt(u))^3)+(za(i)/sqrt(u))^3*u/(u^2+za(i)^2)+za(i)/sqrt(u)*((1-nu)*u/(1+u)+(1+nu)*sqrt(u)*atan(1/sqrt(u))-2));
    taurz = -p0*(ra(i)*za(i)^2/(u^2+za(i)^2)*sqrt(u)/(1+u));
    
    sig = [(sigrr),0,(taurz);0,(sigthth),0;(taurz),0,(sigzz)];
    
    eigens = eig(sig);
    
   sigtest(i) = sigthth;
    maxshear(i) = max([abs(eigens(1)-eigens(2)), abs(eigens(1)-eigens(3)),abs(eigens(3)-eigens(2))]);
end

[C,h]= contour(r/a,z/a,reshape((maxshear/(-p0)),len,wid));
clabel(C,h)

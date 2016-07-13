addpath('C:\Users\Ben\Thesis Iowa\pde2\pde2\noweb');
addpath('C:\Users\Ben\Thesis Iowa\distmesh-master');
addpath('C:\Users\Ben\Thesis Iowa\Code 6-15-16\3D');
addpath('C:\Users\Ben\Thesis Iowa\Code 6-15-16');

fd = @(p) dintersect(drectangle(p,0,.01,-.01,.01),dcircle(p,0,0,.01));
meshsize = .00032; %in meters
%fh = @huniform;
fh = @(p) meshsize + 100*meshsize*max(dcircle(p,0,-.0075,.0025),0);
h = .25*meshsize/6400; %timestep in seconds (CFL condition for elastic wave propogation at speed of sound in steel)
[p,t]=distmesh2d(fd,fh,meshsize,[-.01,-.01;.01,.01],[0,.01;0,-.01]);

np = size(p,1)

%set Lame coefficients
lamb = 121*10^9; %Pascals
mu = 80.8*10^9; %Pascals

lamb2 = lamb*10^(-7);
mu2 = mu*10^(-7);

yieldstr = 2000*10^6; %Pascals

lin2d = lin2d_elt();
lin2dx2 = eltx2_elt(lin2d);
fhtpts = create_fht(p,t,lin2d);
fht = create_fht(p,t,lin2dx2);
nv = fht_num_vars(fht);

grav = 0; %m/s^2
%Body forces and their integrals
f1 = @(x)(0)
f2 = @(x)(-grav)
pde = struct('coeffs',@(x)2*pi*[(lamb+2*mu)/x(1),0,lamb,0,0,lamb;0,0,0,0,0,0; lamb,0,x(1)*(lamb+2*mu),0,0,lamb*x(1);0,0,0,mu*x(1),mu*x(1),0;0,0,0,mu*x(1),mu*x(1),0;lamb,0,lamb*x(1),0,0,x(1)*(lamb+2*mu)],'rhs1',@(x)[2*pi*x(1),0,0,0,0,0;0,2*pi*x(1),0,0,0,0;zeros(4,6)],'rhs2',@(x)[2*pi*x(1)*f1(x);2*pi*x(1)*f2(x);0;0;0;0],'order',1);
% Initialize A and b 
A = sparse(nv,nv); 
M = sparse(nv,nv);
b = zeros(nv,1); 
% Assemble matrices (A and M) and vector b 
[A,M,b] = assemblydyn3Dlin(A,M,b,pde,lin2dx2,p,t,fht,@int2d_radon7);


%Ensure A and M spd by making small adjustments
[row,col] = find(abs(A-A') < abs(A*10e-6) & abs(A-A')>0);
[row2,col2] = find(abs(M-M') < M*10e-9 & abs(M-M')>0);

for i = 1:length(row)
    A(col(i),row(i)) = A(row(i),col(i));
end

for i = 1:length(row2)
    M(col2(i),row2(i)) = M(row2(i),col2(i));
end


pdeD = struct('coeffs',@(x)2*pi*[(lamb2+2*mu2)/x(1),0,lamb2,0,0,lamb2;0,0,0,0,0,0; lamb2,0,x(1)*(lamb2+2*mu2),0,0,lamb2*x(1);0,0,0,mu2*x(1),mu2*x(1),0;0,0,0,mu2*x(1),mu2*x(1),0;lamb2,0,lamb2*x(1),0,0,x(1)*(lamb2+2*mu2)],'rhs1',@(x)[2*pi*x(1),0,0,0,0,0;0,2*pi*x(1),0,0,0,0;zeros(4,6)],'rhs2',@(x)[2*pi*x(1)*f1(x);2*pi*x(1)*f2(x);0;0;0;0],'order',1);
% Initialize A and b 
D = sparse(nv,nv); 
W = sparse(nv,nv);
% Assemble matrices (A and M) and vector b 
[D,W,c] = assemblydyn3Dlin(D,W,b,pdeD,lin2dx2,p,t,fht,@int2d_radon7);



rho = 7850; %density in kg/m^3
M = rho*M;

[Dbedges,Dbnodes,Dt_index,Nbedges,Nbnodes,Nt_index] = Dboundary3dcontact(t,p,2*meshsize); %find boundary specific to a condition
[fullbedges,fullbnodes,fullt_index] = boundary2d(t); %all parts of boundary

%construct list of bv where we want clamping (displacement = 0)
pde2 = struct('coeffs',@(x)[1],'rhs',@(x)0,'order',0);
[AbD,bbD,bvlistD] = assembly2dbdry(pde2,lin2d,p,t,Dbedges,Dt_index,fhtpts,@int1d_gauss5);
[Abfull,bbfull,fullbvlist] = assembly2dbdry(pde2,lin2d,p,t,fullbedges,fullt_index,fhtpts,@int1d_gauss5);



%compute unit outward normals for each Nbedge
normals = edgenormal(p,t,Nbedges,Nt_index);


%Compute the distance in the normal direction from each potential contact node to the
%foundation
[gapf,Nbnodenorms] = gap3(p,Nbnodes,Nbedges, normals,-.01001,fhtpts);


%Compute matrix which will be equivalent to dotting the displacement vector
%with the appropriate normal vector.  This will be premultiplied by
%displacement and velocity vectors when solving

dispdot = dispdotassembly(Nbnodenorms, [length(Nbnodes), nv]);


%%%%TESTING ABSOLUTE DISPLACEMENT IN Y DIRECTION RATHER THAN NORMAL
%%%%DISPLACEMENT
[row,col] = find(dispdot);
for i = 1:length(col)
    if mod(col(i),2)
        dispdot(row(i),col(i))=0;
    else
        dispdot(row(i),col(i))=-1;
    end
end



%define boundary traction conditions and find the integral values on
%Neumann part of boundary
% g1 = @(x)(0);
% g2 = @(x)(0);


%assembly3D assumes that g1 and g2 are 0
pde6 = struct('coeffs1',@(x,n)[2*pi*x(1)*n(1)],'coeffs2',@(x,n)[2*pi*x(1)*n(2)],'order',0);
B2 = assembly3Dbdrycontact(pde6,lin2d,p,t,Nbedges,Nt_index,fhtpts,normals,@int1d_gauss5);


% find non Dirichlet boundary variables (cbvlistx2) and Dboundary variables 
v_array = ones(nv,1); v_array(2*bvlistD-1) = 0; %v_array(2*bvlistD) = 0; 
cbvlistx2 = find(v_array ~= 0);
bvlistx2 = (find(v_array == 0))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u_int = A(cbvlistx2,cbvlistx2) \ (b(cbvlistx2) + bbN(cbvlistx2));
% 
% u = zeros(nv,1); 
% u(cbvlistx2) = u_int; 
% u( bvlistx2) = 0; %these are the places with 0 displacement condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tend = .00006; %end time

initdownv = 2 %in m/s
icond = [zeros(nv,1);(reshape([zeros(1,nv/2);-1*initdownv*ones(1,nv/2)],1,[]))'];


pvlist = get_pvlist(fhtpts,np); %gives list of positions for each node
daspect([1,1,1])


%dispfrac = floor(tend/h/80); %stresses will be computed and plot made for 1 out of each dispfrac timesteps
dispfrac = 5;

%solve for the displacements and max shear stress and make plots  
%solveElast_CGMwithcompl3Dv3(A,M,b,B2,gapf,dispdot,bvlistx2,cbvlistx2,icond,tend,h,pvlist,lamb,mu,p,t,fullbedges,yieldstr,dispfrac,Nbnodes);
[E,delE,maxshear,Ntest,gaptest,EPE2] = solveElast_CGM_KV(A,D,M,b,B2,gapf,dispdot,bvlistx2,cbvlistx2,icond,tend,h,pvlist,lamb,mu,p,t,fullbedges,yieldstr,dispfrac,Nbnodes,rho,0,1);


%surfy(Ntest,gaptest,50,8,dispfrac);
%hertzComp(meshsize, Ntest, dispfrac,h);
%hertzStresses(meshsize,Ntest, dispfrac,h,p,t,fullbedges,yieldstr)

%save N.mat Ntest


%Make Energy Graph
% x=linspace(1,tend/h,tend/h);
% plot(x,E,x,KE,x,EPE,x,GPE)

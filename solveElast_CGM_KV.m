function [E,delE,maxshear, Ntest, gaptest,EPE2] = solveElast_CGM_KV(A,D,M,f,B,gapf,bnodenorms,bvlist,cbvlist,icond,tend,h,pvlist,lamb,mu,p,t,fullbedges,yieldstr, dispfrac,Nbnodes,rho,vidtf,sheartf)
%This function numerically solves the system Mu'' = -Au -Du'+ f - BN using the conjugate
%gradient method build into matlab.  f represents the body forces, BN gives surface forces
%due to contact with a foundation.  These forces are computed by solving a complementarity
%problem.  The resulting displacements are then used to compute stress at
%each node and these stresses are used to make plots of the results

%This solver imposes dirichlet
%conditions passed in as part of icond (the values at time 0 on the
%boundary are taken to be the permanent, fixed boundary values) on the
%nodes in bvlist.

%icond is a vector containing the initial displacements and velocities for
%each node

%tend is the final time (there will be tend/h timesteps)

%stresses will be calculated and a plot will be made for 1 out of each dispfrac timesteps 

nv = length(icond)/2;

shear = zeros(nv/2,1);

u0 = icond(1:nv); %IC for u and u' 
v0 = icond(nv+1:2*nv);

u1 = zeros(nv,1); 
v1 = zeros(nv,1);

u1(bvlist) = u0(bvlist); %dirichlet boundary


 LM=h^2/2*bnodenorms*((M+h^2/4*A)\B); %Trapezoid!!!
 imagesc(LM)
 colorbar
%  savefig(strcat('/space/bdill/Desktop/Code 6-14-16/Animations/LM',int2str(1),'.fig'));
savefig(strcat('C:\Users\Ben\Thesis Iowa\Animations\LM',int2str(1),'.fig'));
 %LM=h^2/2*B'*((M+h^2/4*A)\B);
 
 %LM = h^2*bnodenorms*(M\B); %Implicit Euler!!
 
N = zeros(length(Nbnodes),1);
%dottest = zeros(floor(dispfrac*tend/h)-1,length(N));
Ntest = zeros(floor(dispfrac*tend/h)-1,length(N));
gaptest = zeros(floor(dispfrac*tend/h)-1,length(N));
%Ftest = zeros(floor(dispfrac*tend/h)-1,1);
k = 0;

[sorted,I] = sort(pvlist(Nbnodes));
[sorted,I2] = sort(I);


% filename = strcat('/space/bdill/Desktop/Code 6-14-16/Animations/ball',int2str(1),'.avi');
% filename2 = strcat('/space/bdill/Desktop/Code 6-14-16/Animations/energy',int2str(1),'.fig');
% filename3 = strcat('/space/bdill/Desktop/Code 6-14-16/Animations/maxshear',int2str(1),'.fig');
filename = strcat('C:\Users\Ben\Thesis Iowa\Animations\ball',int2str(1),'.avi');
filename2 = strcat('C:\Users\Ben\Thesis Iowa\Animations\energy',int2str(1),'.fig');
filename3 = strcat('C:\Users\Ben\Thesis Iowa\Animations\maxshear',int2str(1),'.fig');
%filename = strcat('Animations\ball',int2str(1),'.avi');
%filename2 = strcat('Animations\energy',int2str(1),'.fig');


if(vidtf == 1)
vid = VideoWriter(filename);
vid.FrameRate = 60;
open(vid);
end

elasenergy = zeros(nv/2,1);
EPE2 = zeros(floor(dispfrac*tend/h)-1,1);
maxshear = zeros(floor(dispfrac*tend/h)-1,1);
EPE = zeros(floor(dispfrac*tend/h)-1,1);
KE = zeros(floor(dispfrac*tend/h)-1,1);
x=linspace(1,floor(dispfrac*tend/h)-1,floor(dispfrac*tend/h)-1);

delE = zeros(floor(dispfrac*tend/h)-1,1);

%Lqtest=zeros(length(gapf),floor(dispfrac*tend/h)-1);
close all
figure('units','normalized','outerposition',[0 0 1 1])
for i=1:floor(tend/h)-1

%   
     %Lq = gapf - bnodenorms*u(i,:)' - h/2*bnodenorms*v(i,:)' - h/2*bnodenorms*((M+h^2/4*A)\((M-h^2/4*A)*v(i,:)' - h*A*u(i,:)' + h*f));
     
     %Lq = gapf - bnodenorms*(pcg((M+h^2/4*A),(M-h^2/4*A)*u0(:) + h*M*v0(:) + h^2/2*f,[],1000));
     
     
     
     
     
     [dist,flag2] = pcg((M+h^2/4*A),(M-h^2/4*A)*u0(:) + h*M*v0(:) + h^2/2*f,[1e-8],1000,diag(diag(M+h^2/4*A)));  %Trapezoid!!
     %[dist,flag2] = pcg((eye(nv)+h^2*(M\A)),u0+h*v0-h*2*(M\f),[1e-8],1000,diag(diag(eye(nv)+h^2*(M\A)))); %Implicit Euler
     if(flag2)
         flag2
         pause
     end
     Lq = gapf - bnodenorms*(dist);
     
     
    

     
     %Lqtest(:,i) = Lq(I2);
     
     %Lq = gapf - B'*u(i,:)' - h/2*B'*v(i,:)' - h/2*B'*((M+h^2/4*A)\((M-h^2/4*A)*v(i,:)' - h*A*u(i,:)' + h*f));

    N0 = N;
    [N,err] = LCPgaussseidel(LM,Lq,1e-18,N);
    %[N,err] = lemke(LM,Lq,N);
    N1 = N;
    

    %pause
    if(err)
        err
        pause
    end
   
    
    
    
   %[v1,flag] = pcg(M+h^2/4*A, (M-h^2/4*A)*v - h*A*u(i,:)'  + h*f - h*B*N);
   [v1(:),flag] = pcg(M+h/2*D+h^2/4*A, (M-h/2*D-h^2/4*A)*v0(:) - h*A*u0(:)  + h*f - h*B*N,[1e-10],[1000],diag(diag(M+h^2/4*A)));
   if(flag)
       flag
       pause
   end

   v1(bvlist) = 0;
   
   delE(i) = -h/4*(v1+v0)'*B*(N0+N1);
   %Ftest(i)=(sum(M*v1)-sum(M*v0))/h;
   
   u1(cbvlist) = u0(cbvlist)+h/2*(v0(cbvlist)+v1(cbvlist));
   %u1 = u0+h/2*(v0+v1);
   
   v0(:) = v1(:);
   u0(:) = u1(:);
   
   %comptest = [gapf - bnodenorms*u1, N];
   %dottest(i,:) = gapf-bnodenorms*u1.*N;
   Ntest(i,:) = N(I2);
   gaptest0 = gapf-bnodenorms*u1;
   gaptest(i,:) = gaptest0(I2);
   
   %compute stresses and make plots
   if(~mod(i,dispfrac))
       k = k+1;
       
       if(sheartf == 1)
           [sig,eps] = nodestressdyn2dlinv3(u1(2*pvlist-1),u1(2*pvlist),p,t,lamb,mu);
           [sig2,eps2] = nodestressdyn2dlinv3(v1(2*pvlist-1),v1(2*pvlist),p,t,lamb,mu);
            for j = 1:nv/2  %for each node
                [V,Dsig] = eig(sig(:,:,j));
                %shear(i) = abs(eigs(1)-eigs(2))/2;
                shear(j) = max([abs(Dsig(1,1)-Dsig(2,2)), abs(Dsig(1,1)-Dsig(3,3)),abs(Dsig(2,2)-Dsig(3,3))]);
                %shear(j) = sqrt(1/2*((Dsig(1,1)-Dsig(2,2))^2+(Dsig(2,2)-Dsig(3,3))^2+(Dsig(3,3)-Dsig(1,1))^2));
                %shear(j) = sqrt(.5*((eigs(1)-eigs(2))^2+(eigs(2)^2)+(eigs(1)^2)));
                %could also output where value is large at each timestep
                %rather than saving 
                
                xx = Dsig(1,1);
                yy = Dsig(2,2);
                zz = Dsig(3,3);
                R = (xx+yy+zz)/3;
                
                test = 1/2*((xx-yy)^2+(yy-zz)^2+(zz-xx)^2) - yieldstr^2;

                if(test < 0)
                    Psig = [xx,yy,zz];
    
                else
                    Psig = [R,R,R]+[xx-R,yy-R,zz-R]*sqrt(2/3)*yieldstr/norm([xx-R,yy-R,zz-R]);
                end

                
%                 elasenergy(j) = .5*sum(sum( sig(:,:,j).*eps2(:,:,j)));
                elasenergy(j) = .5*sum(sum( (V*(diag(Psig)-Dsig)/V).*eps(:,:,j)));
                %elasenergy(j) = .5*sum(sum( (V*(Dsig)/V).*eps(:,:,j)));

            end
            maxshear(k) = max(shear);
%             maxshear(k) = max(shear(find(shear ~= max(shear))));
%             if(maxshear(k) > .1*yieldstr)
%                 xx = find(shear==maxshear(k));
%                 posx = p(xx,1) + u1(2*xx-1);
%                 posy = p(xx,2) + u1(2*xx);
%             end
       end
        
        if(vidtf == 1)
            %vnodes = sqrt(v1(2*pvlist-1).^2+v1(2*pvlist).^2);
            %subplot(1,2,1)
            scatter(p(:,1)+u1(2*pvlist-1),p(:,2)+u1(2*pvlist),[.1],shear,'filled'); %shear stresses plot!
            %%velocity plot
            %scatter(p(:,1)+u1(2*pvlist-1),p(:,2)+u1(2*pvlist),[10],vnodes,'filled');
        
        
            %scatter(p(:,1)+S(k,2*pvlist-1)',p(:,2)+S(k,2*pvlist)',[10],'filled');

            hold on
%             if(maxshear(k) > .1*yieldstr)
%                 scatter(posx,posy,5,'black');
%             end
%             hold on

        %     quiver(p(fullbedges(:,1),1),p(fullbedges(:,1),2),p(fullbedges(:,2),1)-p(fullbedges(:,1),1),p(fullbedges(:,2),2)-p(fullbedges(:,1),2),0,'MaxHeadSize',0);
        %     hold on


            quiver(p(fullbedges(:,1),1)+u1(2*pvlist(fullbedges(:,1))-1) ,p(fullbedges(:,1),2)+u1(2*pvlist(fullbedges(:,1))) ...
                ,p(fullbedges(:,2),1)+u1(2*pvlist(fullbedges(:,2))-1)-(p(fullbedges(:,1),1)+u1(2*pvlist(fullbedges(:,1))-1)),...
                p(fullbedges(:,2),2)+u1(2*pvlist(fullbedges(:,2)))-(p(fullbedges(:,1),2)+u1(2*pvlist(fullbedges(:,1)))),0,'MaxHeadSize',0);
            hold on
            %tricontour(p,t,shear,[yieldstr/6,  yieldstr/3, yieldstr/2, 2*yieldstr/3, 5*yieldstr/6, yieldstr]);
            tricontour(p,t,shear,[yieldstr/10000, yieldstr]);
            %tricontour(p,t,sqrt(v1(2*pvlist-1).^2 + v1(2*pvlist).^2),0:.5:2);
            %tricontf(p(:,1),p(:,2),t,shear,.4e9:.4e9:4e9);
            axis([-.0005, .0105, -.0105, .01]); %Full half-disc
            axis equal tight
            %axis([-.0001, .0005, -.01001 - .0000025, -0.01001 + .0000025]); %For close up of contact surface
            %daspect([1,1,1]);
            caxis([0,yieldstr]); %Colors of nodes based on yield strength in Tresca sense
            %caxis([0,1]); %Colors of nodes based on mag. of velocity of nodes
            
            grid off;
            title(strcat('t = ',sprintf('%0.7f',i*h)));
            xlabel('r (meters)')
            ylabel('z (meters)')
            drawnow

    %%%%%%% For making animations of plot output
%         frame = getframe(1);
%         im = frame2im(frame);
%         [AA,map] = rgb2ind(im,256); 
%         
%     	%imwrite(AA,map,filename);
% 
%         
%         if k == 1;
%           imwrite(AA,map,filename,'gif', 'Loopcount',inf,'DelayTime',.001);
%         else
%           imwrite(AA,map,filename,'gif','WriteMode','append','DelayTime',.001);
%         end
        
            hold off

%             subplot(1,2,2)
%             scatter(p(Nbnodes),N(I2));
%             axis([0 .5e-3 0 1e10])
            %pause

            drawnow

%         frame = getframe(1);
%         im = frame2im(frame);
%         [AA,map] = rgb2ind(im,256); 
        
            if k == 1;
              set(gca,'nextplot','replacechildren'); 
            else
              frame = getframe(1);
              writeVideo(vid,frame);
            end
            
        end
        
         KE(k) = .5*v1'*M*v1;
         EPE(k) = .5*u1'*A*u1;

%         KE(k) = .5*rho*bodyintegrator3d(v1(2*pvlist-1).^2 + v1(2*pvlist).^2,p,t,@int2d_radon7);
         EPE2(k) = bodyintegrator3d(elasenergy,p,t,@int2d_radon7);
         %EPE2 is the attempt at quantifying energy loss to plastic def.
         
        
                
   end
   
end

if(vidtf == 1)
    close(vid);
end

figure(2)
plot(x(1:k),KE(1:k)+EPE(1:k), x(1:k), KE(1:k), x(1:k), EPE(1:k)) 

savefig(filename2)
E = KE + EPE;

figure(3)
plot(x(1:k),maxshear(1:k),x(1:k),yieldstr*ones(k,1))
savefig(filename3)

figure(4)
plot(x(1:k),EPE2(1:k))

% figure(3)
% plot(x(1:i),Ftest(1:i))
% frame = getframe(3);
% im = frame2im(frame);
% [AA,map] = rgb2ind(im,256); 
% imwrite(AA,map,filename3,'png');

% imagesc(Lqtest(1:2,:));
% colorbar
end
function[] = hertzStresses(meshsize,N, dispfrac,h,p,t,fullbedges,yieldstr)
    close all
    figure('units','normalized','outerposition',[0 0 1 1])
    
    nu = .3;
    k = 100;
    A = .7/(80.8*10^9);
    
    [len,wid] = size(N);
    first = zeros(wid,1);
    last = zeros(wid,1);
    maxshear = zeros(length(p),1);
    
    for i = 1:wid
        if(any(N(:,i)))
            maxnodes = i;
            first(i) = find(N(:,i),1,'first');
            last(i) = find(N(:,i),1,'last');
        end
        
        
    end
    
    P = zeros(floor((first(1)+last(1))/25),1);
    
    mid = (first(1) + last(1))/2;
    

    maxa = maxnodes*meshsize;
    maxp = max([max(N),14693*10^9*maxa]);
    
    %filename = strcat('/space/bdill/Desktop/Code 6-14-16/Animations/hertzstress',int2str(1),'.avi');
    filename = strcat('C:\Users\Ben\Thesis Iowa\Animations\hertzStress',int2str(1),'.avi');
    
    vid = VideoWriter(filename);
    vid.FrameRate = 60;
    open(vid);
    
    
    m = 0;
    for i = first(1):last(1)
        if(~mod(i,dispfrac))
            m = m+1;
            numContactNodes = find(N(i,:),1,'last');

            if(numContactNodes > 0)
                if(numContactNodes < maxnodes && i < mid)
                    a = meshsize/(first(numContactNodes+1) - first(numContactNodes))*(i - first(numContactNodes)) + meshsize*(numContactNodes-1)+.00001; %meters (radius of contact area)
                elseif(numContactNodes == maxnodes)
                    a = (numContactNodes-1)*meshsize + .00001;
                else
                    a = meshsize/(last(numContactNodes+1) - last(numContactNodes))*(i - last(numContactNodes)) + meshsize*(numContactNodes-1)+.00001;
                end
           
                P(m) = 8*k*a^3/(3*A);
                p0 = 14693*10^9*a;
                
                for j = 1: length(p(:,1))
                
                    ra = (p(j,1)+.0001)/a;
                    za = (p(j,2)+.0101)/a;
                    
                    %start fixing hereeeee
                    u = (.5*(ra^2+za^2-1+sqrt((ra^2+za^2-1)^2+4*za^2)));
                    sigthth = -p0*((1-2*nu)/3*1/ra^2*(1-(za/sqrt(u))^3)+za/sqrt(u)*(2*nu+(1-nu)*u/(1+u)-(1+nu)*sqrt(u)*atan(1/sqrt(u))));
                    sigzz = -p0*((za/sqrt(u))^3*u/(u^2+za^2));
                    sigrr = p0*((1-2*nu)/3*1/ra^2*(1-(za/sqrt(u))^3)+(za/sqrt(u))^3*u/(u^2+za^2)+za/sqrt(u)*((1-nu)*u/(1+u)+(1+nu)*sqrt(u)*atan(1/sqrt(u))-2));
                    taurz = -p0*(ra*za^2/(u^2+za^2)*sqrt(u)/(1+u));

                    
                    
                    
%                     sigthth = -p0*((1-2*nu)/3*1./ra.^2.*(1-(za./sqrt(u)).^3)+za./sqrt(u).*(2*nu+(1-nu).*u./(1+u)-(1+nu)*sqrt(u).*atan(1./sqrt(u))));
%                     sigzz = -p0*((za./sqrt(u)).^3.*u./(u.^2+za.^2));
%                     sigrr = p0*((1-2*nu)/3*1./ra.^2.*(1-(za./sqrt(u)).^3)+(za./sqrt(u)).^3.*u./(u.^2+za.^2)+za./sqrt(u).*((1-nu)*u./(1+u)+(1+nu)*sqrt(u).*atan(1./sqrt(u))-2));
%                     taurz = -p0*(ra.*za.^2./(u.^2+za.^2).*sqrt(u)./(1+u));
    
                    
                    
                        sig = [(sigrr),0,(taurz);0,(sigthth),0;(taurz),0,(sigzz)];
    
                        eigens = eig(sig);
  
                        maxshear(j) = max([abs(eigens(1)-eigens(2)), abs(eigens(1)-eigens(3)),abs(eigens(3)-eigens(2))]);
                    
                

                
                end
%                 -.01+.48*a
%                 p(find(maxshear == max(maxshear)),:)
%                 pause
                
                
%               scatter3(p(:,1),p(:,2),maxshear)
%               pause
                %scatter([p(:,1);-p(:,1)],[p(:,2);p(:,2)],[1],[maxshear;maxshear],'filled');
                scatter(p(:,1),p(:,2),[1],maxshear,'filled');
                %scatter(r,z,[10],maxshear,'filled');
                hold on
%                 quiver(p(fullbedges(:,1),1) ,p(fullbedges(:,1),2) ...
%                 ,p(fullbedges(:,2),1)-(p(fullbedges(:,1),1)),...
%                 p(fullbedges(:,2),2)-(p(fullbedges(:,1),2)),0,'MaxHeadSize',0);
                hold on
                %tricontour([p;-p(:,1),p(:,2)],[t;t+length(p)],[maxshear;maxshear],[yieldstr/6,  yieldstr/3, yieldstr/2, 2*yieldstr/3, 5*yieldstr/6, yieldstr]);
                tricontour(p,t,maxshear,[yieldstr/6,  yieldstr/3, yieldstr/2, 2*yieldstr/3, 5*yieldstr/6, yieldstr]);
                axis([-.0005, .0105, -.0105, .01]); %Full half-disc
                axis equal tight
                caxis([0,yieldstr]);
                grid off;
                title(strcat('Hertzian Max Shear Stress Contours at t = ',sprintf('%0.7f',i*h)));
                xlabel('r (meters)')
                ylabel('z (meters)')
%                 ylabel('Pressure (Pascals)')
%                 legend('Simulation','Hertzian')
                drawnow 

                if i == first(1);
                  set(gca,'nextplot','replacechildren'); 
                else
                  frame = getframe(1);
                  writeVideo(vid,frame);
                end

                hold off
            else
                continue
            end
        end
    end
    close(vid)
    
    length(P)
    length(h*dispfrac*(0:(length(P)-1)))
    
    plot(h*dispfrac*(0:(length(P)-1)),P)
end



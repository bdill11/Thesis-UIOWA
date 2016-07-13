function[] = hertzComp(meshsize,N, dispfrac,h)
    close all
    
    [len,wid] = size(N);
    first = zeros(wid,1);
    last = zeros(wid,1);
    
    for i = 1:wid
        if(any(N(:,i)))
            maxnodes = i;
            first(i) = find(N(:,i),1,'first');
            last(i) = find(N(:,i),1,'last');
        end
        
        
    end
    
    mid = (first(1) + last(1))/2;
    

    maxa = maxnodes*meshsize;
    maxp = max([max(N),14693*10^9*maxa]);
    
    vid = VideoWriter(strcat('C:\Users\Ben\Thesis Iowa\Animations\hertzcomp',int2str(1),'.avi'));
    vid.FrameRate = 60;
    open(vid);
    
    

    for i = first(1):last(1)
        if(~mod(i,dispfrac))
            numContactNodes = find(N(i,:),1,'last');

            if(numContactNodes > 0)
                if(numContactNodes < maxnodes && i < mid)
                    a = meshsize/(first(numContactNodes+1) - first(numContactNodes))*(i - first(numContactNodes)) + meshsize*(numContactNodes-1)+.00001; %meters (radius of contact area)
                elseif(numContactNodes == maxnodes)
                    a = (numContactNodes-1)*meshsize + .00001;
                else
                    a = meshsize/(last(numContactNodes+1) - last(numContactNodes))*(i - last(numContactNodes)) + meshsize*(numContactNodes-1)+.00001;
                end
                r = 0:a/30:a;
                p = 14693*10^9*a*(1-r.^2/a^2).^.5;
                x = 0:meshsize:meshsize*numContactNodes;
                plot(x,N(i,1:length(x)),r,p)
                axis([0 1.2*maxa 0 1.2*maxp])
                title(strcat('Pressure Distributions at t = ',sprintf('%0.7f',i*h)));
                xlabel('r (meters)')
                ylabel('Pressure (Pascals)')
                legend('Simulation','Hertzian')
                drawnow 

                if i == first(1);
                  set(gca,'nextplot','replacechildren'); 
                else
                  frame = getframe(1);
                  writeVideo(vid,frame);
                end


            else
                continue
            end
        end
    end
    close(vid)
end



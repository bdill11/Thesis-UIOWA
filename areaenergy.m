totarea = 0;
for i=1:959
    area = (EPE2(i)+EPE2(i+1))/2*(h*dispfrac);
    
    totarea = totarea + area;
    
end
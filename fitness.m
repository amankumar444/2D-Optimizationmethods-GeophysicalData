function [fit]= fitness(pop,popsize,check)
f=[383.7 1828 8600 41290 192650];
h=63.72;
r=[6.87000000 6.73000000 6.59000000 6.68000000 6.64000000];
low = 0.3;
high = 0.02;
for i = 1:popsize
    d = pop(i,1:2);
    rho= pop(i,3:5);
    out(i,:)= forward_HEM(rho,d,h,r,f);

fit(i,:) = 1/(((mean((out(i,:)-check)./check).^2)+0.001))^(1/4);

end

end
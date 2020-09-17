function [a]=delta(n)
%n=length(m);
a=eye(n);
for i=1:n-1
    a(i+1,i)=-1;
end

a(1,1)=0;
end
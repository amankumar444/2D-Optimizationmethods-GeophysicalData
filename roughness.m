function [a]=roughness(m,q)
al = length(m)-1 ;
a=zeros(q*al,al+1);
ini = 1;
fin = q+1;
for i=1:q*al
    a(i,ini)=1;
    a(i,fin) = -1;
    ini = ini+q;
    fin = fin+q;
    if(fin>=length(m))
        ini = 1;
        fin = 1+q;
    end
end
size(a)
end
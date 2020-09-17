function [pop] = cross(pop,s,crossrate)
newpop = pop;
for i = 1:2:length(s)
    
        a= pop(s(i),:);
        b= pop(s(i+1),:);
        x=a;
        y=b;
        if rand(1)<crossrate
        
        r=rand(1,5)<0.5;
        for j= 1:5
               if r(j)==1
                   x(j)=b(j);
                    y(j)=a(j);
                    
               end
        end
        
        end
        newpop(s(i),:)=x;
        newpop(s(i+1),:)=y;
        
        
        
        
        
end
pop = newpop;
end
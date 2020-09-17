function [s] = select(fit)

    s= 1:length(fit);
    
    df=cumsum(fit);
    
    r = rand(length(fit));
    r = r*sum(fit);
    for i = 1:length(fit)
    
        s(i)=sum(r(i)>df)+1;
        
    end
end
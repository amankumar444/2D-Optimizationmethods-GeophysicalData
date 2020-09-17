%function to calculate Jacobian Matrix 
function [Resistivity,J,Re] = Jack(rho,d,h,r,f,m)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Re=forward_HEM(rho,d,h,r,f);
J=zeros(length(Re),length(m));


for i =1:length(m)
    if i<=length(rho)
        res=m(i);
        rho(i)=1.05*res;
    else
        dep=m(i);
        d(i)=1.05*dep;
    end
    
    Re1=forward_HEM(rho,d,h,r,f);
    J(:,i)= (Re1-Re)/0.05;
end
    Resistivity = rho;
end
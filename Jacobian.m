
function J = Jacobian(rho,d,h,r,f,m)

m_ini = m;
Re=forward_HEM(rho,d,h,r,f);

J=zeros(length(Re),length(m));


for i =1:length(m)
    m = m_ini;
    m(i) = m(i)*1.05;
    
    Re1=forward_HEM(m(1:length(rho)),m(length(rho+1):length(m)),h,r,f);
    
    J(:,i)= (Re1-Re)/(0.05*m(i));
end
end
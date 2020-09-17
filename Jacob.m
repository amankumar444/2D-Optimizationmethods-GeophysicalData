
function [J, K] = Jacob(h,r,f,m,n,q,p_1,p_2,p_3)
K = zeros(10*n);
h0 = p_1;
r0 = p_2;
f0 = p_3;
m_ini = m;
dp =[];
rho = ones(1,(q+1)/2);
for i = 1:n
    h = h0(i);
    r = r0(:,i)';
    f = f0(:,i)';
    m_f = m((i-1)*q+1:i*q);
    g = forward_HEM(m_f(1:length(rho)),m_f(length(rho)+1:length(m_f)),h,r,f);
    dp = [dp
          g];
end
Re=dp;

K = Re;
Re1 = [];
J=zeros(length(Re),length(m));
for i =1:n
    for j = 1:q
        m = m_ini;
        Re1 = [];
        m_f_1 = m((i-1)*q+1:i*q);
        m_f_1(j) = m_f_1(j)*1.05;
        Re1=forward_HEM(m_f_1(1:length(rho)),m_f_1(length(rho+1):length(m_f_1)),h0(i),r0(:,i)',f0(:,i)');
    
        J((i-1)*10+1:i*10,(i-1)*q+j)= (-1*(Re((i-1)*10+1:i*10)-Re1))/(0.05*m_f_1(j));
     end
end
end

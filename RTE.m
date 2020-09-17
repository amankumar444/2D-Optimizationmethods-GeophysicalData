%Fuction to calculate B(N) values for given input 
%n= number of layers
%l=lambda values
%d = layer thickness
%w2 = frequency
%a = alpha-0
function [B] = RTE(n,lamda,mew,ep,rho,d,w2)


i=sqrt(-1);

%defining reflection coefficient 
% B=zeros(n,length(l));
% A=zeros(n,length(l));

B(n)=sqrt(lamda.^2-(w2).^2*mew(n)*ep(n)+i*(mew(n)/rho(n)*w2));%Appendix-A (Formula A9) B(n)for last layer
%A(n)=B(n);%A(n) for last layer 

for k1=n-1:-1:1
		A(k1)=sqrt(lamda.^2-(w2).^2*mew(k1)*ep(k1)+i*(mew(k1)/rho(k1)*w2));% Appendix-A (Formula A9)for A(n-1) to A(1)
		B(k1)=A(k1).*(B(k1+1)+A(k1).*tanh(A(k1).*d(k1)))./(A(k1)+B(k1+1).*tanh(A(k1).*d(k1)));%for B(n-1) to B(1)
%     R(k)=(mew(1)*B(k)-mew(2)*a)./(mew(1)*B(k)+mew(2)*a);
end
end





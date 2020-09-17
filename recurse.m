
function [B] = recurse(n,lamda,mew,ep,rho,d,w2)


i=sqrt(-1);


B(n)=sqrt(lamda.^2-(w2).^2*mew(n)*ep(n)+i*(mew(n)/rho(n)*w2));
A(n)=B(n);

for k1=n-1:-1:1
		A(k1)=sqrt(lamda.^2-(w2).^2*mew(k1)*ep(k1)+i*(mew(k1)/rho(k1)*w2));
		B(k1)=A(k1).*(B(k1+1)+A(k1).*tanh(A(k1).*d(k1)))./(A(k1)+B(k1+1).*tanh(A(k1).*d(k1)));

end


end
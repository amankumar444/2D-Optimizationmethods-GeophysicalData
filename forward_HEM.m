function P = forward_HEM(rho,d,h,r,f)

rho0=10^9;%air resistivity value

J= [1.15E-06
1.45E-06
1.83E-06
2.30E-06
2.90E-06
3.65E-06
4.59E-06
5.78E-06
7.28E-06
9.17E-06
1.15E-05
1.45E-05
1.83E-05
2.30E-05
2.90E-05
3.65E-05
4.59E-05
5.78E-05
7.28E-05
9.17E-05
1.15E-04
1.45E-04
1.83E-04
2.30E-04
2.90E-04
3.65E-04
4.59E-04
5.78E-04
7.28E-04
9.17E-04
1.15E-03
1.45E-03
1.83E-03
2.30E-03
2.90E-03
3.65E-03
4.59E-03
5.78E-03
7.28E-03
9.16E-03
1.15E-02
1.45E-02
1.83E-02
2.30E-02
2.89E-02
3.63E-02
4.55E-02
5.69E-02
7.10E-02
8.81E-02
1.08E-01
1.31E-01
1.55E-01
1.76E-01
1.86E-01
1.70E-01
1.03E-01
-3.03E-02
-2.28E-01
-3.62E-01
-2.06E-01
3.37E-01
3.18E-01
-5.14E-01
3.09E-01
-1.27E-01
4.62E-02
-1.81E-02
8.35E-03
-4.47E-03
2.62E-03
-1.60E-03
9.98E-04
-6.26E-04
3.94E-04
-2.49E-04
1.57E-04
-9.89E-05
6.24E-05
-3.94E-05
2.48E-05
-1.57E-05
9.89E-06
-6.24E-06
3.94E-06
-2.48E-06
1.57E-06
-9.89E-07
6.24E-07
-3.94E-07
2.48E-07
-1.57E-07
];%filter coefficients

w=2*pi*f;
n=length(rho);%number of layers >2

%c0=2.99792458*10^8;%speed of EM wave in vacuum
mew0=4*pi*10^-7; % Vacuum permeability
mew=mew0*ones(1,n);%can be changed later
ep0=8.854187817*10^-12; %Vacuum permittivity

ep=ep0*ones(1,n);%can be changed later

%defining Hankel transform

j0=54;
m1=10;%number of nodes per decade (Formula 10 Fast Henkel Transform)

lamda=zeros(1,length(J));

% l = lamda;
R=zeros(1,length(f));


F=zeros(n,length(f));
Z=zeros(1,length(f));
Z1=zeros(1,length(f));
Ip=zeros(1,length(f));
Qp=zeros(1,length(f));

i=sqrt(-1);
for m=1:length(f) %frequency loop 
    w2=w(m);%Individual frequency for each iteration
    temp1=0;% temporary to perform summation
    for k=1:length(J)%filter loop
        lamda(k)=(1/8)*(10.^((k-j0)/m1));

         [B] = RTE(n,lamda(k),mew,ep,rho,d,w2);%calculation of B(1)for formula A8 
         alpha0=sqrt(lamda(k).^2-(w2).^2.*mew0.*ep0+i.*(w2).*mew0./rho0);
         
         R=(mew0*B(1)-mew(1)*alpha0)/(mew0*B(1)+mew(1)*alpha0);%to determine Complex Reflection coefficient Formula A-8
         Z(k) = R.*lamda(k).^3.*J(k)/(alpha0).*exp(-2*alpha0*h);%using Formula 11            
		  %Z(k)=F.*J(k);%using formula 10 for each j values
         
          temp1=temp1+Z(k);%temporary variable for summation using Formula 10

    end 
  %end of filter loop
          Z1(m)=r(m).^2.*temp1.*10^6;% using formula 10 & 11
 temp1 = 0; %initializing for new iteration
		
          

Ip(m)=real(Z1(m));%Inphase
Qp(m)=imag(Z1(m));%Quadrature
end

 P=[Ip Qp]';
end
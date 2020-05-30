
% LPF stub impedance calculation

fc = 2*10^9;            %cutoff freq
wc = (2*pi)*(fc);
R0 = 50;                %normalized to 50 Ohms
N=5;                    %5th order
a = 1:N;
b = 1:N;
g = 1:N;
ripple = .2;
beta = log(coth(ripple/17.37));
y = sinh(beta/(2*N));
%calculate a, b, and g values
for i = 1:N
   a(i) = sin(((2*i-1)*pi)/(2*N)); 
   b(i) = (y^2) + (sin((i*pi)/N))^2;
   if(i == 1)
    g(1) = (2*a(1))/y;
   else
    g(i) = ((4*a(i-1)*a(i))/(b(i-1)*g(i-1)));
   end
end

%lumped elements
C1= g(1)/(R0*wc);
L1 = (g(2)*R0)/wc;
C2= g(3)/(R0*wc);
L2 = (g(4)*R0)/wc;
C3= g(5)/(R0*wc);

%kuroda identities
n1 = 1 + (g(1));                    %1 + cap/tline
n2 = 1 + n1;                        %1 + tline/ind
n3 = 1 + g(1)/(n1*g(2));            %1 + tline/ind
Tc1 = 1/n2;                         %first stub
Ttl1 = n2/n1;                       %first tline
Tc2 = n1/(n3*g(1));                 %second stub
Ttl2 = n3*g(2);                     %second tline

%final impedances
Z1 = 50*Tc1;
Z2 = 50*Ttl1;
Z3 = 50*Tc2;
Z4 = 50*Ttl2;
Z5 = 50/g(3);
%because of symmetry...
Z6 = Z4;
Z7 = Z3;
Z8 = Z2;
Z9 = Z1;






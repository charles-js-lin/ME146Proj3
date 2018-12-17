%% Task 1 Problem Statement
% Find molar air-to-fuel ratio alpha for i) pure propane, ii) pure methane,
%   iii) 50% mixture

%% Initializing Constants/Assumptions

Tr = 550; %[K]
Tp = 1900; %[K]
Tavg = (Tr + Tp)/2;
To = 25+273; %[K] from 25 oC

% Fuel Propane ratios for i) pure propane, ii) pure methane, 
%   iii) 50 percent mixture
g1 = 1;
g2 = 0;
g3 = 0.5;

%% Finding h_i(T) Molar Enthalpy for Gas Species i
% cp = a + bT + cT^2 + dT^3, for T = Tavg
% 6 Gas species: 
%   1) C3H8
%   2) CH4
%   3) O2 
%   4) H2O 
%   5) CO2 
%   6) N2

a = [-4.04 19.89 25.48 32.24 22.26 28.90];
b = [30.48 5.024 1.520 0.1923 5.981 -0.1571].*10^-2;
c = [-15.72 1.269 -0.7155 1.055 -3.501 0.8081].*10^-5;
d = [31.74 -11.01 1.312 -3.595 7.469 -2.873].*10^-9;

cp=zeros(1,6);
for i=1:length(cp)
    cp(i) = a(i) + b(i)*Tavg + c(i)*Tavg^2 + d(i)*Tavg^3; %[kJ/(kmol*K)]
end

ho = [-103850 -74850 0 -241820 -393520 0]; %[kJ/kmol]
h = @(n,T) ho(n) + cp(n)*(T-To); %Molar Enthalpy

%% Finding alpha_i for Cases (i), (ii), (iii)

a1 = 4.76*(-g1*h(1,Tr) - (1-g1)*h(2,Tr) - (3*g1+2)*h(3,Tp) + ...
    (2+2*g1)*h(4,Tp) + (1+2*g1)*h(5,Tp)) / ...
    (h(3,Tr) + 3.76*h(6,Tr) - 3.76*h(6,Tp) - h(3,Tp));

a2 = 4.76*(-g2*h(1,Tr) - (1-g2)*h(2,Tr) - (3*g2+2)*h(3,Tp) + ...
    (2+2*g2)*h(4,Tp) + (1+2*g2)*h(5,Tp)) / ...
    (h(3,Tr) + 3.76*h(6,Tr) - 3.76*h(6,Tp) - h(3,Tp));

a3 = 4.76*(-g3*h(1,Tr) - (1-g3)*h(2,Tr) - (3*g3+2)*h(3,Tp) + ...
    (2+2*g3)*h(4,Tp) + (1+2*g3)*h(5,Tp)) / ...
    (h(3,Tr) + 3.76*h(6,Tr) - 3.76*h(6,Tp) - h(3,Tp));

results = [g1 g2 g3; a1 a2 a3]'; %Left is fuel propane fraction, right is alpha
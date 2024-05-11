% UNIQUAC PROBLEM
% DATE : 19-01-21
% NAME : SANTOSH M. ZOL

clear all
clc
format short

P = 100000; % in pascals
T = 345; % in kelvin
R = 8.314; % in joule.mol^(-1).K^(-1)
z = 10; % co-ordination number

x = [0.2 0.3 0.5]; % mole fraction
r = [2.5755 3.1878 2.7694]; % volume parameter
q = [2.588 2.634 2.400]; % area paramter
u = [1.0000 1.2160 0.2030;0.6170 1.0000 0.0480;0.8380 0.6120 1.0000]; % interaction parameter

% calculating value of tao
[row,col] = size(u);
tao = zeros(row,col);

for i=1:row
    for j=1:col
        tao(j,i) = exp(-(u(j,i)-u(i,i))/(R*T));
    end
end

theta = zeros(1,row);
phi = zeros(1,row);
l = zeros(1,row);

% calculating value of theta, phi and l
for i = 1:row
    theta(1,i) = (x(1,i)*q(1,i)/sum(x.*q));
    phi(1,i) = ((x(1,i)*r(1,i))/sum(x.*r));
    l(1,i) = (z/2)*(r(1,i) - q(1,i)) - (r(1,i)-1); 
end


% calculating residual part
for i=1:row
    sumj1 = 0; sumj2 = 0;
    for j=1:row
        sumk = 0;
        for k = 1:row
            sumk = sumk + theta(1,k)*tao(k,j);
        end
        sumj1 = sumj1 + (theta(1,j)*tao(i,j))/sumk;
        sumj2 = sumj2 + theta(1,j)*tao(j,i);
    end
    lngama_R(1,i)= q(1,i)*(1-log(sumj2)-sumj1);
end

% calculating combinatorial part
xl = sum(x.*l);
for i=1:row
    lngama_C(1,i) = log(phi(1,i)/x(1,i)) + (0.5*z)*q(1,i)*log(theta(1,i)/phi(1,i))+l(1,i)-(phi(1,i)/x(1,i))*xl;
end

for i=1:row
    lngama(1,i) = lngama_C(1,i) + lngama_R(1,i);
end

fprintf('Activity coefficient of Benzene is  %f\n',exp(lngama(1,1)));
fprintf('Activity coefficient of Toluene is  %f\n',exp(lngama(1,2)));
fprintf('Activity coefficient of Water is  %f\n',exp(lngama(1,3)));

clear; close all; clc
ks = linspace(0,10,5);
x = linspace(-10,10,100);

figure;
for i = 1:length(ks)
    k = ks(i);
    rho = WelshRho(x,k);
    plot(x,rho, 'DisplayName', ['k = ', num2str(k)]);
    hold on
    legend
end


function rho = WelshRho(x, k)
%Welsh function (simplified GGW function)
rho = 1-exp( (-(x/k).^2)/2 );
    
end


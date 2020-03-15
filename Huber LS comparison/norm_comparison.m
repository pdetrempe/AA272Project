% Script for plotting norm comparisons
% Paul DeTrempe
clear;clc; close all;

x = linspace(-3,3,1000);
k_Huber = 1.345;
k_Welsh = 1;

rho_Huber = HuberRho(x, k_Huber);
rho_Welsh = WelshRho(x, k_Welsh);
rho_l2 = TwoNormRho(x);
rho_l1 = OneNormRho(x);

%% Plots
figure;
plot(x, rho_l2,'Linewidth',2)
hold on
plot(x, rho_l1,'Linewidth',2)
plot(x, rho_Huber,'Linewidth',2)
plot(x, rho_Welsh,'Linewidth',2)
xlabel('x')
ylabel('\rho')
legend('L_2', 'L_1', 'Huber', 'Welsh')


%% Functions
function rho = HuberRho(x, k)
    rho = zeros(size(x));
    for i = 1:length(x)
       if abs(x(i)) < k
           rho(i) = x(i)^2;
       else
           rho(i) = 2*k*(abs(x(i))-k/2);
       end
    end
end

function rho = WelshRho(x, k)
%Welsh function (simplified GGW function)
rho = (1-exp( (-(x/k).^2) ));
    
end

function rho = OneNormRho(x)
    rho = abs(x);
end

function rho = TwoNormRho(x)
    rho = x.^2;
end
% Demonstration of Huber norm vs. naive LS
% Paul DeTrempe
clear; close all; clc;

%% Load data
clear;
load('satpos.mat');

%% Part a) Write the code for iterative least-squares using all measurements from a timestep. The output should be the final estimate of position or position+clockbias (ECEF system).
x0 = zeros(3,1); % meters
b0 = 0;

% repeat for pranges200.mat
measurement200 = load('pranges200.mat');
time = 200;

% Form rho0 and rhoC vectors
for i = 1:satpos.numSats
    rho0(i,:) = norm([satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]-x0);
    x_sat_matrix(i,:) = [satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]';
end

rhoC = measurement200.prange';
[x_LS,b_LS,dx_history,~] = getPositionLS(x_sat_matrix, rhoC, rho0, x0, b0);
x_LS
b_LS
format short


% % Plot residual history
% figure;
% semilogy(dx_history, 'DisplayName', '\delta x')
% hold on
% legend
% title('Q4a) Residuals vs. Iteration, Least Squares Positioning')
% xlabel('Iteration')
% ylabel('Norm Error, $\mid \mid \delta* \mid \mid_2$', 'Interpreter', 'LaTeX')

%% Part a.ii) Vary initial conditions

% repeat for pranges200.mat
measurement200 = load('pranges200.mat');
time = 200;

% Form rho0 and rhoC vectors
for i = 1:satpos.numSats
    rho0(i,:) = norm([satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]-x0);
    x_sat_matrix(i,:) = [satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]';
end

rhoC = measurement200.prange';
[x_Huber,b_Huber,dx_history,~] = getPositionHuber(x_sat_matrix, rhoC, rho0, x0, b0);
x_Huber
b_Huber
format short

%% Functions
function [x_LS,b_LS,dx_history,db_history] = getPositionLS(x_sat_matrix, rhoC, rho0, x0, b0)
y0 = [x0;b0];
d_rho = rhoC-rho0;

% iterate until convergence to obtain position and bias information
iter = 1;
max_iter = 8;
epsilon = 1e-4; % meters

dx_history = [];
db_history =[];
while iter<max_iter && norm(d_rho) > epsilon
    % form matrix
    G = getGeometryMatrix(x_sat_matrix, x0);
    
    % update guess
    dy = G\d_rho;
    
    y0 = y0 + dy;
    x0 = y0(1:end-1);
    b0 = y0(end);
    
    for i = 1:size(x_sat_matrix, 1)
        rho0(i) = norm(x_sat_matrix(i,:)'-x0) + b0;
    end
    d_rho = rhoC-rho0;
    
    % Store history for plotting
    dx = norm(dy(1:end-1));
    db = abs(dy(end));
    dx_history = [dx_history;dx];
    db_history = [db_history;db];
    
    % break criteria
    iter = iter + 1;
end
x_LS = x0;
b_LS = b0;

end

function [x_Huber,b_Huber,dx_history,db_history] = getPositionHuber(x_sat_matrix, rhoC, rho0, x0, b0)
y0 = [x0;b0];
d_rho = rhoC-rho0;

% iterate until convergence to obtain position and bias information
iter = 1;
max_iter = 8;
epsilon = 1e-4; % meters

dx_history = [];
db_history =[];
while iter<max_iter && norm(d_rho) > epsilon
    % form matrix
    G = getGeometryMatrix(x_sat_matrix, x0);
    
    % update guess
    cvx_begin quiet
        variable dy(4)
        M = 1;
        f = huber_circ(d_rho-G*dy,M);
        minimize(f);
    cvx_end
    
    y0 = y0 + dy;
    x0 = y0(1:end-1);
    b0 = y0(end);
    
    for i = 1:size(x_sat_matrix, 1)
        rho0(i) = norm(x_sat_matrix(i,:)'-x0) + b0;
    end
    d_rho = rhoC-rho0;
    
    % Store history for plotting
    dx = norm(dy(1:end-1));
    db = abs(dy(end));
    dx_history = [dx_history;dx];
    db_history = [db_history;db];
    
    % break criteria
    iter = iter + 1;
end
x_Huber = x0;
b_Huber = b0;

end



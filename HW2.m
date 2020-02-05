% AA272C HW #2
% Paul DeTrempe
clear; close all; clc;

%% Question 3
load('satpos.mat');

%% Part a) Obtain positions for all satellites across all timesteps in ENU frame w.r.t. given reference position
x_ref_ECEF = [-2699.908; -4292.867; 3855.203]; % [km]
x_ref_LLH = [37.427611; -122.166971; 25.592];  % Lat [deg], East-Long [deg], altitude [meters]
x_ref_LLH(1:2) = deg2rad(x_ref_LLH(1:2));       % radians
x_ref_LLH(3) = x_ref_LLH(3)/1000;               % km

R_ecef2enu = getECEF2ENU(x_ref_LLH);

% Initialize new fields for storing data
r_ENU = zeros(3, satpos.numSats, satpos.numEpochs);
elevations = zeros(satpos.numEpochs, satpos.numSats);
inView = zeros(satpos.numEpochs, satpos.numSats);

for i = 1:satpos.numEpochs
    for j = 1:satpos.numSats
        x_ECEF = [satpos.x(i,j); satpos.y(i,j); satpos.z(i,j)];
        
        % rotate to ENU coordinates
        r_ENU(:,j,i) = R_ecef2enu*x_ECEF;
        
        % calculate elevation
        [~,elev,~] = enu2aer(r_ENU(1,j,i), r_ENU(2,j,i), r_ENU(3,j,i));
        elevations(i,j) = elev;
        
    end
end

%% Part b) Calculate elevation angles for all satellites across all timesteps w.r.t. given reference position. Plot elevation angles for satellite 7, 14, 21 and 28.
figure;
plot(satpos.time, elevations(:,7), 'DisplayName', 'SV 7')
hold on;
plot(satpos.time, elevations(:,14), 'DisplayName', 'SV 14')
plot(satpos.time, elevations(:,21), 'DisplayName', 'SV 21')
plot(satpos.time, elevations(:,28), 'DisplayName', 'SV 28')
legend
title('Question 3b, Elevation vs. Time')
xlabel('Time, [hours]')
ylabel('Elevation, [degrees]')

%% Part c) Calculate and plot number of visible satellites for masks of 10, 15, 20 and 25 degrees
masking_angles = [10; 15; 20; 25]; % degrees
PDOP = zeros(satpos.numEpochs, length(masking_angles));
HDOP = zeros(satpos.numEpochs, length(masking_angles));
VDOP = zeros(satpos.numEpochs, length(masking_angles));
TDOP = zeros(satpos.numEpochs, length(masking_angles));
figure;
for l = 1:length(masking_angles)
    mask = masking_angles(l);
    % find SVs in view
    inView = elevations>mask;
    numInView = sum(inView,2);
    
    % plot number in view
    subplot(4,1,l)
    plot(satpos.time, numInView)
    title(['SVs in view vs. time, masking angle = ', num2str(mask), ' degrees'])
    xlabel('Time, [hours]')
    ylabel('Number of visible SVs')
    
    % Compute PDOP, HDOP, and VDOP using visible satellites
    for i = 1:satpos.numEpochs
        % Go through each sat, if it's in view, add an entry to the G
        % matrix
        x_sat_matrix = [];
        for j = 1:satpos.numSats
            if inView(i,j)
                x_sat_matrix = [x_sat_matrix; r_ENU(:,j,i)'];
            end
        end
        x0 = zeros(3,1);
        G = getGeometryMatrix(x_sat_matrix, x0);
        
        % Compute PDOP HDOP and VDOP from G matrix
        H = inv(G'*G);
        PDOP(i,l) = sqrt(norm([H(1,1); H(2,2); H(3,3)]));
        HDOP(i,l) = sqrt(norm([H(1,1); H(2,2)]));
        VDOP(i,l) = sqrt(abs(H(3,3)));
        TDOP(i,l) = sqrt(abs(H(4,4)));
    end
end

%% Part d) Compute and plot PDOP, HDOP and VDOP for each mask
% Plots
figure;
for l = 1:length(masking_angles)
    subplot(4,1,l)
    plot(satpos.time, PDOP(:,l), 'DisplayName','PDOP')
    hold on
    plot(satpos.time, HDOP(:,l),  'DisplayName','HDOP')
    plot(satpos.time, VDOP(:,l),  'DisplayName','VDOP')
%     plot(satpos.time, TDOP(:,l),  'DisplayName','TDOP')
    legend
    title(['Mask angle = ', num2str(masking_angles(l)), 'degrees'])
    xlabel('Time, hours')
    ylabel('DOP')
end

%% Part e: Correlation between DOPs and # of SVs
for l = [1,4]
    figure;
    yyaxis left;
    plot(satpos.time, PDOP(:,l), 'DisplayName','PDOP')
    hold on
    plot(satpos.time, HDOP(:,l),  'DisplayName','HDOP')
    plot(satpos.time, VDOP(:,l),  'DisplayName','VDOP')
%     plot(satpos.time, TDOP(:,l),  'DisplayName','TDOP')
    legend
    title(['Mask angle = ', num2str(masking_angles(l)), 'degrees'])
    xlabel('Time, hours')
    ylabel('DOP')
    
    yyaxis right;
    mask = masking_angles(l);
    % find SVs in view
    inView = elevations>mask;
    numInView = sum(inView,2);
    
    % plot number in view
    plot(satpos.time, numInView, 'DisplayName', '# of SVs')
    ylabel('Visible SVs')
end

%% Question 4
clear;
measurement10 = load('ranges10.mat');
load('satpos.mat');

%% Part a) Write the code for iterative least-squares using all measurements from a timestep. The output should be the final estimate of position or position+clockbias (ECEF system).
x0 = zeros(3,1) % meters
b0 = 0

% rho0 = zeros(satpos.numSats,1);
% x_sat_matrix = zeros(satpos.numSats,3);
time = 10;

% Form rho0 and rhoC vectors
for i = 1:satpos.numSats
    rho0(i,:) = norm([satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]-x0);
    x_sat_matrix(i,:) = [satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]';
end

rhoC = measurement10.prange';
format long
[x_LS,dx_history] = getPositionLS3D(x_sat_matrix, rhoC, rho0, x0);
x_LS

% repeat for ranges100.mat
measurement100 = load('ranges100.mat');
time = 100;

% Form rho0 and rhoC vectors
for i = 1:satpos.numSats
    rho0(i,:) = norm([satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]-x0);
    x_sat_matrix(i,:) = [satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]';
end

rhoC = measurement100.prange';
[x_LS,dx_history] = getPositionLS3D(x_sat_matrix, rhoC, rho0, x0);
x_LS

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

x0 = 1e7*ones(3,1) % meters
b0 = 1e4

rhoC = measurement10.prange';
format long
time = 10;

% Form rho0 and rhoC vectors
for i = 1:satpos.numSats
    rho0(i,:) = norm([satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]-x0);
    x_sat_matrix(i,:) = [satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]';
end
[x_LS,dx_history] = getPositionLS3D(x_sat_matrix, rhoC, rho0, x0);
x_LS

% repeat for ranges100.mat
measurement100 = load('ranges100.mat');
time = 100;

% Form rho0 and rhoC vectors
for i = 1:satpos.numSats
    rho0(i,:) = norm([satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]-x0);
    x_sat_matrix(i,:) = [satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]';
end

rhoC = measurement100.prange';
[x_LS,dx_history] = getPositionLS3D(x_sat_matrix, rhoC, rho0, x0);
x_LS

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

%% Part b) Repeat above with 5, 10, 20 and 31 measurements at a time.
b0 = 0;

numSats = [5;10;20;31];
SVIDs = 1:satpos.numSats;
% figure;
for i = 1:length(numSats)
    % select a random permutation of available satellites
    numSVs = numSats(i)
    indices = datasample(SVIDs,numSVs,'Replace',false);
    
    rhoC = measurement10.prange';
    time = 10;

    % Form rho0 and rhoC vectors
    for i = 1:satpos.numSats
        rho0(i,:) = norm([satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]-x0);
        x_sat_matrix(i,:) = [satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]';
    end

    format long
    [x_LS,~] = getPositionLS3D(x_sat_matrix(indices,:), rhoC(indices), rho0(indices), x0);
    x_LS

    % repeat for ranges100.mat
    measurement100 = load('ranges100.mat');
    time = 100;

    % Form rho0 and rhoC vectors
    for i = 1:satpos.numSats
        rho0(i,:) = norm([satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]-x0);
        x_sat_matrix(i,:) = [satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]';
    end

    rhoC = measurement100.prange';
    [x_LS,~] = getPositionLS3D(x_sat_matrix(indices,:), rhoC(indices), rho0(indices), x0);
    x_LS

    % repeat for pranges200.mat
    measurement200 = load('pranges200.mat');
    time = 200;

    % Form rho0 and rhoC vectors
    for i = 1:satpos.numSats
        rho0(i,:) = norm([satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]-x0);
        x_sat_matrix(i,:) = [satpos.x(time,i); satpos.y(time,i); satpos.z(time,i)]';
    end

    rhoC = measurement200.prange';
    [x_LS,b_LS,~,~] = getPositionLS(x_sat_matrix(indices,:), rhoC(indices), rho0(indices), x0,b0);
    x_LS
    b_LS
end

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

function [x_LS,dx_history] = getPositionLS3D(x_sat_matrix, rhoC, rho0, x0)
y0 = x0;
d_rho = rhoC-rho0;

% iterate until convergence to obtain position and bias information
iter = 1;
max_iter = 8;
epsilon = 1e-4; % meters

dx_history = [];
while iter<max_iter && norm(d_rho) > epsilon
    % form matrix
    G = getGeometryMatrix3D(x_sat_matrix, x0);
    
    % update guess    
    for i = 1:size(x_sat_matrix, 1)
        rho0(i) = norm(x_sat_matrix(i,:)'-x0);
    end
    d_rho = rhoC-rho0;
    
    dy = G\d_rho;
    
    y0 = y0 + dy;
    x0 = y0(1:end);


%     norm(d_rho)
    
    % Store history for plotting
    dx = norm(dy(1:end));
    dx_history = [dx_history;dx];

    
    % break criteria
    iter = iter + 1;
end
x_LS = x0;

end


function G = getGeometryMatrix(x_sat_matrix, x0)
% takes in Nx3 satellite matrix and guess at position, returns geometry
% matrix
N = size(x_sat_matrix,1);
G = [];
for i = 1:N
    LoS = getLoSVector(x_sat_matrix(i,:)', x0);
    G = [G; -LoS', 1];
end

end

function G = getGeometryMatrix3D(x_sat_matrix, x0)
% takes in Nx3 satellite matrix and guess at position, returns geometry
% matrix
N = size(x_sat_matrix,1);
G = [];
for i = 1:N
    LoS = getLoSVector(x_sat_matrix(i,:)', x0);
    G = [G; -LoS'];
end

end

function LoS = getLoSVector(x_satellite, x0)
% Use column vectors
LoS = (x_satellite - x0)/norm(x_satellite-x0);
end

function R = getECEF2ENU(x_ref_LLH)
% Takes in reference location in LLH [rad, km] and outputs rotation
% matrix from ECEF to ENU
lat = x_ref_LLH(1);
long = x_ref_LLH(2);

R = [-sin(long), cos(long), 0;...
    -sin(lat)*cos(long), -sin(lat)*sin(long), cos(lat);...
    cos(lat)*cos(long), cos(lat)*sin(long), sin(lat)];
end


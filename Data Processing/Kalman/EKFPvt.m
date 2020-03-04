function [x0,H,P] = EKFPvt(prs,gpsEph,x0,P)
% [xHat,z,svPos,H,Wpr,Wrr] = WlsPvt(prs,gpsEph,x0)
% calculate a weighted least squares PVT solution, xHat
% given pseudoranges, pr rates, and initial state
%
% Inputs:
%  prs: matrix of raw pseudoranges, and pr rates, each row of the form:
%  [trxWeek,trxSeconds,sv,prMeters,prSigmaMeters,prrMps,prrSigmaMps]
%   trxWeek, trxSeconds: Rx time of measurement
%      where trxSeconds = seconds in the current GPS week
%   sv: satellite id number
%   prMeters, prSigmaMeters: pseudorange and standard deviation (meters)
%   prrMps, prrSigmaMps: pseudorange rate and standard deviation (m/s)
%   gpsEph: matching vector of GPS ephemeris struct, defined in ReadRinexNav
%   x0: initial (previous) state, [x,y,z,bc,xDot,yDot,xDot,bcDot]'
%       in ECEF coordinates(meters and m/s)
%
% Outputs: xHat: estimate of state update
%          z = [zPr; zPrr] a-posteriori residuals (measured-calculated)
%          svPos: matrix of calculated sv positions and sv clock error:
%                 [sv prn, x,y,z (ecef m), dtsv (s),xDot,yDot,zDot, dtsvDot]
%          H: H observation matrix corresponding to svs in svPos(:,1)
%          Wpr,Wrr Weights used in WlsPvt = 1/diag(sigma measurements)
%                  use these matrices to compute variances of xHat
%
% The PVT solution = x0 + xHat, in ECEF coordinates
% For unweighted solution, set all sigmas = 1

%Author: Frank van Diggelen
%Open Source code for processing Android GNSS Measurements

% Modified by Zack Miller for Stanford AA 272c Winter 2020 for EFK

jWk=1; jSec=2; jSv=3; jPr=4; jPrSig=5; jPrr=6; jPrrSig=7;%index of columns

[bOk,numVal] = checkInputs(prs, gpsEph, x0);
if ~bOk
    error('inputs not right size, or not properly aligned with each other')
end

xHat=[]; z=[]; H=[]; svPos=[];
xyz0 = x0(1:3);
bc = x0(4);

if numVal<4
    return
end
ttxWeek = prs(:,jWk); %week of tx. Note - we could get a rollover, when ttx_sv
%goes negative, and it is handled in GpsEph2Pvt, where we work with fct
ttxSeconds =  prs(:,jSec) - prs(:,jPr)/GpsConstants.LIGHTSPEED; %ttx by sv clock
% this is accurate satellite time of tx, because we use actual pseudo-ranges
% here, not corrected ranges
% write the equation for pseudorange to see the rx clock error exactly cancel
% to get precise GPS time: we subtract the satellite clock error from sv time,
% as done next:
dtsv = GpsEph2Dtsv(gpsEph,ttxSeconds);
dtsv = dtsv(:); %make into a column for compatibility with other time vectors
ttx = ttxSeconds - dtsv; %subtract dtsv from sv time to get true gps time

%calculate satellite position at ttx
[svXyzTtx,dtsv,svXyzDot,dtsvDot]=GpsEph2Pvt(gpsEph,[ttxWeek,ttx]);
svXyzTrx = svXyzTtx; %initialize svXyz at time of reception

for i=1:length(gpsEph)
    % calculate tflight from, bc and dtsv
    dtflight = (prs(i,jPr)-bc)/GpsConstants.LIGHTSPEED + dtsv(i);
    % Use of bc: bc>0 <=> pr too big <=> tflight too big.
    %   i.e. trx = trxu - bc/GpsConstants.LIGHTSPEED
    % Use of dtsv: dtsv>0 <=> pr too small <=> tflight too small.
    %   i.e ttx = ttxsv - dtsv
    svXyzTrx(i,:) = FlightTimeCorrection(svXyzTtx(i,:), dtflight);
end

%EKF Calculations ---------------------------------------------------
% Noise Parameters --------------------------------------------------

Q = 1 * eye(4);
R = 1e-5 * eye(numVal);

% Jacobians ---------------------------------------------------------

F = eye(4);

v = xyz0(:)*ones(1,numVal,1) - svXyzTrx';%v(:,i) = vector from sv(i) to xyz0
range = sqrt( sum(v.^2) );
v = v./(ones(3,1)*range); % line of sight unit vectors from sv to x0
H = [v', ones(numVal,1)];

% Predict -----------------------------------------------------------

x0(1:4) = x0(1:4);
P = F*P*F' + Q;

% Update ------------------------------------------------------------

% zPr = prs(:,jPr) - range' - x0(4);   %innovation
prHat = range(:) + bc - GpsConstants.LIGHTSPEED*dtsv;
zPr = prs(:,jPr)-prHat; 
K = P*H' * inv(R + H*P*H');

x0 = x0 + K * zPr;
P = (eye(size(P)) - K*H) * P * (eye(size(P)) - K*H)' + K*R*K';

end %end of function EKFPvt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [bOk,numVal] = checkInputs(prs, gpsEph, x0)
%utility function for WlsPvt
jWk=1; jSec=2; jSv=3; jPr=4; jPrSig=5; jPrr=6; jPrrSig=7;%index of columns

bOk=false;
%check inputs
numVal=size(prs,1);
if (max(prs(:,jSec))-min(prs(:,jSec)))> eps
    return
elseif length(gpsEph)~=numVal
    return
elseif any(prs(:,jSv) ~= [gpsEph.PRN]')
    return
elseif  any(size(x0) ~= [4,1])
    return
elseif size(prs,2)~=7
    return
else
    bOk = true;
end

%We insist that gpsEph and prs are aligned first.
%ClosestGpsEph.m does this, and passes back indices for prs - this is the way to
%do it right, so we don't have nested searches for svId

end %end of function checkInputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2016 Google Inc.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.


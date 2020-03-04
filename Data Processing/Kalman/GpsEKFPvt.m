function gpsPvt = GpsEKFPvt(gnssMeas,allGpsEph,bRaw)
%gpsPvt = GpsWlsPvt(gnssMeas,allGpsEph,bRaw)
%compute PVT from gnssMeas
% Input: gnssMeas, structure of pseudoranges, etc. from ProcessGnssMeas
%        allGpsEph, structure with all ephemeris
%        [bRaw], default true, true => use raw pr, false => use smoothed
%
% Output: 
% gpsPvt.FctSeconds    Nx1 time vector, same as gnssMeas.FctSeconds
%       .allLlaDegDegM Nx3 matrix, (i,:) = [lat (deg), lon (deg), alt (m)]
%       .sigmaLlaM     Nx3 standard deviation of [lat,lon,alt] (m)
%       .allBcMeters   Nx1 common bias computed with llaDegDegM
%       .allVelMps     Nx3 (i,:) = velocity in NED coords
%       .sigmaVelMps   Nx3 standard deviation of velocity (m/s)
%       .allBcDotMps   Nx1 common freq bias computed with velocity
%       .numSvs        Nx1 number of satellites used in corresponding llaDegDegM
%       .hdop          Nx1 hdop of corresponding fix
%
%Algorithm: Weighted Least Squares

%Author: Frank van Diggelen
%Open Source code for processing Android GNSS Measurements

% Modified by Zack Miller for Stanford AA 272c Winter 2020 for EFK

if nargin<3
    bRaw = true;
else
    %check that smoothed pr fields exists in input
    if any(~isfield(gnssMeas,{'PrSmM','PrSmSigmaM'}))
       error('If bRaw is false, gnssMeas must have fields gnssMeas.PrSmM and gnssMeas.PrSmSigmaM')
    end
end

x0 =zeros(4,1);%initial state: [center of the Earth, bc=0, velocities = 0]'

weekNum = floor(gnssMeas.FctSeconds/GpsConstants.WEEKSEC);
%TBD check for week rollover here (it is checked in ProcessGnssMeas, but
%this function should stand alone, so we should check again, and adjust 
%tRxSeconds by +- a week if necessary)
%btw, Q. why not just carry around fct and not worry about the hassle of
%weeknumber, and the associated week rollover problems?
% A. because you cannot get better than 1000ns (1 microsecond) precsision
% when you put fct into a double. And that would cause errors of ~800m/s * 1us
% (satellite range rate * time error) ~ 1mm in the range residual computation
% So what? well, if you start processing with carrier phase, these errors
% could accumulate.

N = length(gnssMeas.FctSeconds);
gpsPvt.FctSeconds      = gnssMeas.FctSeconds;
gpsPvt.allLlaDegDegM   = zeros(N,3)+NaN; 
gpsPvt.sigmaLLaM       = zeros(N,3)+NaN;
gpsPvt.allBcMeters     = zeros(N,1)+NaN;
gpsPvt.allVelMps       = zeros(N,3)+NaN;
gpsPvt.sigmaVelMps     = zeros(N,3)+NaN;
gpsPvt.allBcDotMps     = zeros(N,1)+NaN;
gpsPvt.numSvs          = zeros(N,1);
gpsPvt.hdop            = zeros(N,1)+inf;

for i=1:N
    iValid = find(isfinite(gnssMeas.PrM(i,:))); %index into valid svid
    svid    = gnssMeas.Svid(iValid)';
    
    [gpsEph,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas.FctSeconds(i));
    svid = svid(iSv); %svid for which we have ephemeris
    numSvs = length(svid); %number of satellites this epoch
    gpsPvt.numSvs(i) = numSvs;
    if numSvs<4
        continue;%skip to next epoch
    end
    
    %% PVT -----------------------------------------------------------------
    %for those svIds with valid ephemeris, pack prs matrix for WlsNav
    prM     = gnssMeas.PrM(i,iValid(iSv))';
    prSigmaM= gnssMeas.PrSigmaM(i,iValid(iSv))';
    
    prrMps  = gnssMeas.PrrMps(i,iValid(iSv))';
    prrSigmaMps = gnssMeas.PrrSigmaMps(i,iValid(iSv))';
    
    tRx = [ones(numSvs,1)*weekNum(i),gnssMeas.tRxSeconds(i,iValid(iSv))'];
    
    prs = [tRx, svid, prM, prSigmaM, prrMps, prrSigmaMps];
    
    % NEW CODE
    if i==1
        [x0,~,~,H,Wpr,~] = WlsPvt(prs,gpsEph,zeros(8,1));%compute WLS solution for timestep 1
        P = inv(H'*(Wpr'*Wpr)*H); 
    else
        [x0,~,P] = EKFPvt(prs,gpsEph,x0(1:4),P);%compute EKF solution
    end
    % END OF NEW CODE
    
    %extract position states
    llaDegDegM = Xyz2Lla(x0(1:3)');
    gpsPvt.allLlaDegDegM(i,:) = llaDegDegM;
    gpsPvt.allBcMeters(i) = x0(4);
    
    %%end WLS PVT --------------------------------------------------------------
end

end %end of function GpsWlsPvt
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


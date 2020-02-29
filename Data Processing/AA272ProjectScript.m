%AA272ProjectScript.m, script to read GnssLogger output

%Author: Frank van Diggelen
% Modified Winter 2020 Zack Miller
%Open Source code for processing Android GNSS Measurements

clear; clc; close all;
param.llaTrueDegDegM = [];

%% data
% save data from GnssLogger App, and edit dirName and prFileName appropriately
dirName = [pwd,filesep,'..',filesep,'Log Files'];
prFileName = 'Rooftop.txt'; %Rooftop, Walking, Huang
addpath(genpath(pwd));

%% parameters
%param.llaTrueDegDegM = [];
%enter true WGS84 lla, if you know it:
% param.llaTrueDegDegM = [37+25/60+36.85322/3600, -(122+10/60+23.79153/3600), 32];%Durand Roof

%% Set the data filter and Read log file
dataFilter = SetDataFilter;
[gnssRaw,gnssAnalysis] = ReadGnssLogger(dirName,prFileName,dataFilter);
if isempty(gnssRaw), return, end

%% Get online ephemeris from Nasa ftp, first compute UTC Time from gnssRaw:
fctSeconds = 1e-3*double(gnssRaw.allRxMillis(end));
utcTime = Gps2Utc([],fctSeconds);
allGpsEph = GetNasaHourlyEphemeris(utcTime,dirName);
if isempty(allGpsEph), return, end

%% process raw measurements, compute pseudoranges:
[gnssMeas] = ProcessGnssMeas(gnssRaw);

%% plot pseudoranges and pseudorange rates
% h1 = figure;
% [colors] = PlotPseudoranges(gnssMeas,prFileName);
% h2 = figure;
% PlotPseudorangeRates(gnssMeas,prFileName,colors);
% h3 = figure;
% PlotCno(gnssMeas,prFileName,colors);

%% compute WLS position and velocity
tic()
gpsWLSPvt = GpsWlsPvt(gnssMeas,allGpsEph);
toc()
%% plot Pvt results
h4 = figure;
ts = 'Weighted Least Squares solution';
PlotPvt(gpsWLSPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;
ax = gca;

%% compute Huber position and velocity
tic()
gpsHuberPvt = GpsHuberPvt(gnssMeas,allGpsEph);
toc()

%% plot Pvt results
h6 = figure;
ts = 'Raw Pseudoranges, Huber solution';
PlotPvt(gpsHuberPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;
ax = gca;

%% Compute Welsh position
tic()
gpsWelshPvt = GpsWelshPvt(gnssMeas,allGpsEph);
toc()

%% plot Pvt results
h6 = figure;
ts = 'Raw Pseudoranges, Welsh solution';
PlotPvt(gpsWelshPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;
ax = gca;

% Plot Accumulated Delta Range 
if any(any(isfinite(gnssMeas.AdrM) & gnssMeas.AdrM~=0))
    [gnssMeas]= ProcessAdr(gnssMeas);
    h6 = figure;
    PlotAdr(gnssMeas,prFileName,colors);
    [adrResid]= GpsAdrResiduals(gnssMeas,allGpsEph,param.llaTrueDegDegM);drawnow
    h7 = figure;
    PlotAdrResids(adrResid,gnssMeas,prFileName,colors);
end
%% end of ProcessGnssMeasScript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2016 Google Inc.
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

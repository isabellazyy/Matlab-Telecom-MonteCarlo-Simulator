clear all, 
close all,
clc

%% Monte Carlo code - starting point of the thesis
% Here, you can change the number of iterations, cells and users,
% size of the simulation area can also be changed 
numOfIterations = 30;
numOfUsers = 50;
numOfPaths = 10;
carrierFreq = 2600e6; % Hz

%% System Parameters

txPower = 30; %dBm
txAntennaGain = 3; %dBi
cableLoss = 0; %dB
EIRP = txPower + txAntennaGain - cableLoss;

BoltzmanConstant = 1.3806503*10^(-23);
NoisedBm = 10*log10(BoltzmanConstant*290*180000*50/0.001);

RxAntennaGain = 0; %dBi
 
% 10 MHz bandwidth, it has 50 PRB
% 2x2 MIMO configuration
% TP = Nprb * BWprb * BWeff * log2(1 +SINR / SINReff)
BWeff = 0.42;
SINReff = 1.1;
Nprb = 50;
BWprb = 180000;

% criteria for the communication
% Gmax = 150;
% RSSmin = -150;
SINRmin = -10;
Smax = 7.67;  % spectral efficiency
  
%matrix is 2D
pathloss = antennaData();
%matrix is 3D
pathloss = repmat(pathloss,1,1,numOfIterations);
%fast fading is created
%ff = FastFading(numOfPaths,carrierFreq);

rayleigh = load('/home/taylan/Documents/MATLAB/monte/FF');
ff = rayleigh.FF;
% FF = zeros(size(pathloss));

[numOfCells, numOfPix, numOfIterations] = size(pathloss);
%fast fading is added to each cell and pixel relation for each iteration
for i=1:numOfIterations
    FF(:,:,i) = ff(randi(size(ff),numOfCells,numOfUsers));
end

%totalLoss = pathloss + abs(FF);
%users = randi([1,length(totalLoss)],1,numOfUsers,numOfIterations);
users = randi([1,length(pathloss)],1,numOfUsers,numOfIterations);

RSSs = zeros(numOfCells,length(users(1,:,1)),numOfIterations);

for i=1:numOfIterations
    %RSSs(:,:,i) = EIRP + RxAntennaGain - totalLoss(:,users(1,:,i),i);    
    RSSs(:,:,i) = EIRP + RxAntennaGain - pathloss(:,users(1,:,i),i);    

end

RSSmax = max(RSSs); % max RSS in dB

% Which cell serves which users
servingCells = zeros(numOfCells,numOfUsers,numOfIterations); % this matrix can be called map between users and cells
tempRSSs = repmat(RSSmax,numOfCells,1);
servingCells = (RSSs == tempRSSs);

% After cell association, fast fading is also subtracted from received
% signal strength
RSSs = RSSs - abs(FF);
RSSmax = max(RSSs);

% converting to linear scale
RSSs_Linear = 10 .^(RSSs / 10);
RSSmax_Linear =  10.^(RSSmax/10);

% received interference in linear scale
Interference = sum(RSSs_Linear) - RSSmax_Linear;           
Noise = 10 ^ (NoisedBm/10);

 %% SINR calculations

temp_Cell = sum(servingCells);

SINRdB = zeros(1,numOfUsers,numOfIterations);

SINRdB = 10*log10(RSSmax_Linear ./ (Interference + Noise) .* servingCells);
indices = find(SINRdB < -10);
SINRdB(indices) = 0;
SINRdB = sum(SINRdB);
 
roundRobin = zeros(numOfCells,numOfUsers,numOfIterations);    % this matrix gives us the shares between users

howManyUsers = sum(permute(servingCells,[2,1,3]));
howManyUsers = permute(howManyUsers,[2,1,3]);
roundRobin = Nprb * servingCells ./ howManyUsers;
indices = find(isnan(roundRobin));
roundRobin(indices) = 0;

%% Throughput calculations

A = zeros(1,numOfUsers,numOfIterations);
A = BWeff * log2(1 +(10 .^(SINRdB/10)) ./ SINReff);
indices = find(A > Smax);
A(indices) = Smax;
TP = max(roundRobin * BWprb .* A) / (10^6);   % All user throughputs for the iterations

%% Jain's index 

Fairness = zeros(1,numOfIterations);
Fairness = rdivide(sum(TP).^2, sum(TP.^2)) / numOfUsers;
 
 
%% SINR-CDF plot

figure

subplot(3,1,1)
 
SINRdB_vector = reshape(SINRdB,[1,numOfUsers*numOfIterations]);

cdfplot(SINRdB_vector);
xlabel('SINR [dB]')
ylabel('CDF')
title('SINR vs CDF')

%% Throughput-CDF plot
subplot(3,1,2)

TP_vector = reshape(TP,[1,numOfUsers*numOfIterations]);

cdfplot(TP_vector);
xlabel('TP [Mbps]')
ylabel('CDF')
% axis([0 maxTP 0 1]);
title('TP vs CDF');

%% Jain's index
subplot(3,1,3)

cdfplot(Fairness)
xlabel('Fairness')
ylabel('CDF')

% save the SINR
save('/home/taylan/Documents/MATLAB/subfolder/ALL_SINR_Values','SINRdB');


clear all, 
close all,
clc

%% Monte Carlo code - starting point of the thesis
% Here, you can change the number of iterations, cells and users,
% size of the simulation area can also be changed 
numOfIterations = 30;
numOfCells = 50;
numOfUsers = 300;
numOfPaths = 20;
meter = 2000; % meter
carrierFreq = 800e6; % Hz

%% System Parameters

txPower = 30; %dBm
txAntennaGain = 3; %dBi
cableLoss = 0; %dB
EIRP = txPower + txAntennaGain - cableLoss;

BoltzmanConstant = 1.3806503*10^(-23);
NoisedBm = 10*log10(BoltzmanConstant*290*180000*50/0.001);

RxAntennaGain = 0; %dBi
PenetrationLoss = 0; % not considered here since users are in open square

APL = EIRP + RxAntennaGain - PenetrationLoss;

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
   
% distance function brings the distances between users and cells
distance = distance(numOfUsers,numOfCells,numOfIterations,meter);
 
%% matrices are 3D
 
PathLoss = losses(distance,APL,carrierFreq);    % after SF and G
rayleigh = FastFading(numOfPaths,carrierFreq); 
FF = rayleigh(randi(size(rayleigh),numOfCells,numOfUsers,numOfIterations));
RSSs = APL - PathLoss - FF;
RSSmax = max(RSSs);                  % max RSS in dB

% Which cell serves which users
servingCells= zeros(numOfCells,numOfUsers,numOfIterations); % this matrix can be called map between users and cells

for i = 1:numOfIterations
    for k = 1:numOfCells
        for j=1:numOfUsers
            if RSSmax(1,j,i) == RSSs(k,j,i)
                servingCells(k,j,i) = 1;
            else
                servingCells(k,j,i) = 0;
            end
        end
    end
end

% converting to linear scale
RSSs_Linear = 10 .^(RSSs / 10);
RSSmax_Linear =  10.^(RSSmax/10);

% received interference in linear scale
Interference = sum(RSSs_Linear) - RSSmax_Linear;           
Noise = 10 ^ (NoisedBm/10);

 %% SINR calculations

temp_Cell = sum(servingCells);

SINRdB = zeros(1,numOfUsers,numOfIterations);
 
for i=1:numOfIterations
% by using the changed servingCells, we can calculate the SINR 
    for j=1:numOfUsers
    
         if temp_Cell(1,j,i) == 1  % is there any cell serving the user
             SINR = RSSmax_Linear(1,j,i) / (Interference(1,j,i) + Noise); 
             SINRdB(1,j,i) = 10 * log10(SINR);
         
             if SINRdB(1,j,i) < -10
                SINRdB(1,j,i) = 0;
             end
         else
            SINRdB(1,j,i) = 0;
         end 
    end
end

roundRobin = zeros(numOfCells,numOfUsers,numOfIterations);    % this matrix gives us the shares between users

for i = 1:numOfIterations
    
    howManyUsers = sum(servingCells(:,:,i)'); % this vector gives us the number of users for each cell

    temp = howManyUsers';       % transpose 

    for j = 1:numOfUsers
        roundRobin(:,j,i) = Nprb * servingCells(:,j,i) ./ temp;
    end         
end


%% Throughput calculations

A = zeros(1,numOfUsers,numOfIterations);

for i=1:numOfIterations
    for j=1:numOfUsers
        A(1,j,i) = BWeff * log2(1 +(10 ^(SINRdB(1,j,i)/10)) / SINReff);
    
        if A(1,j,i) > Smax
            A(1,j,i) = Smax;
        end
    
    end
end
TP = max(roundRobin * BWprb .* A) / (10^6);   % All user throughputs for the iterations

%% Jain's index 

Fairness = zeros(1,numOfIterations);

for i=1:numOfIterations
    Fairness(i)  =  sum(TP(:,:,i))^2 / (numOfUsers) / sum(TP(:,:,i).^2);
end
 
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
% save all the folder
save('/home/taylan/Documents/MATLAB/subfolder/ALL_SINR_Values','SINRdB');
% save the SINR

clear all, 
close all,
clc

%% Monte Carlo code - starting point of the thesis

numOfIterations = 100;
numOfCells = 20;
numOfUsers = 300;
meter = 2000; % meter

%% Open Square ABG model is used for path loss calculations
% Formula ==> PL(f,d)[dB] = 10 * a * log10(d) + B + 10 * y* log10(f) + SF

a = 4.14;
B = 3.66;
y = 2.43;
SF = 7;
f = 0.8; % GHz
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
% TP = Nprb * BWprb * BWeff * log2(1 +SINR / SINReff)
BWeff = 0.4;
SINReff = 1.5;
Nprb = 50;
BWprb = 180;

%% Random locations
% this section can be done in another M file
% Rectangular is chosen; e.g. width = 150m, height = 150m
% cell locations are chosen for only one time
% C =[C1x, C1y; C2x, C2y]
% U is the matrix for the locations of the users. User 
% locations are randomly chosen for each iteration
% U =[U1x, U1y; U2x, U2y; U3x, U3y]

distance = [zeros(numOfCells,numOfUsers)]; % distance between Users and Cells
   
C = randi([0 meter],numOfCells,2); % cell random locations 
% AllTP = zeros(numOfIterations,numOfUsers); % Throughput list for each iterations
SINR_All = zeros(numOfIterations,numOfUsers);

%% Monte Carlo Starts
for counter = 1:numOfIterations

    % here the user locations are created
    U = randi([0 meter],numOfUsers,2);
   
   for i = 1:numOfCells
      for k = 1:numOfUsers
          temp = C(i,:) - U(k,:);
          distance(i,k) = sqrt(sum((temp.^2)));   % distance between users and cells
      end 
   end
    
% Path Loss matrix ~ G(distance) 

G = zeros(size(distance));
G = 10 * a * log10(distance) + B + 10 * y* log10(f) + SF;

% Received Powers from the cells
RSSs = APL - G ;    % Received Signal Strengths in dB
RSSmax = max(RSSs);

% Which cell serves the user

        
[x y] = size(RSSs);
servingCells= zeros(x,y); % this matrix can be called map

for i = 1:y
    for k = 1:x
        if RSSmax(i) == RSSs(k,i)
            servingCells(k,i) = 1;
        else
            servingCells(k,i) = 0;
        end
    end
end

% if the received RSS is lower than threshold, no communication
for i=1 : length(RSSmax)   
    if RSSmax(i) < -100
       RSSmax(i) = 0; 
    end    
end

% here is for calculations, no special meaning
G = reshape(G,[1, numOfUsers*numOfCells]);
temp_Cell = reshape(servingCells,[1,numOfUsers*numOfCells]);
% here is for calculations, no special meaning


% if the path loss is highher than threshold, neglect the channel
% In order to do this, we need to understand which cell has the best power
% for a user and then we need to see the path loss for this channel
for i = 1:numOfUsers*numOfCells
    if temp_Cell(i) == 1 && G(i) > 70
        temp_Cell(i) = 0;
    end
end

G = reshape(G,[numOfCells, numOfUsers]);
servingCells = reshape(temp_Cell,[numOfCells,numOfUsers]);
% so, servingCells is changed 


% converting to linear scale
RSSs_Linear = (10 .^(RSSs / 10));
RSSmax_Linear =  10.^(RSSmax/10);

Interference = sum(RSSs_Linear) - RSSmax_Linear;           % received interference in linear scale
Noise = (10 ^ (NoisedBm/10));

temp_Cell = sum(servingCells);

SINRdB = zeros(1,numOfUsers);

% by using the changed servingCells, we can calculate the SINR
for i= 1:length(servingCells)
    
    if temp_Cell(i) == 1
        SINR = RSSmax_Linear(i) / (Interference(i) + Noise);
        SINRdB(i) = 10 * log10(SINR);
    else
        SINRdB(i) = 0;
    end     
end

% A cell may serve more than one user so calculate how many users are
% served by one cell, then we can decide the scheduling
% Therefore, transpose of the layout is a need
layoutTrans = servingCells';
howManyUsers = sum(layoutTrans); % this vector gives us the number of users for each cell

temp = howManyUsers';       % transpose 
roundRobin = zeros(x,y);    % this matrix gives us the shares between users

for i = 1:y
    roundRobin(:,i) = servingCells(:,i) ./ temp;
end

roundRobin = roundRobin * Nprb; % how many PRB for each user
              
% U1,U2,U3,...,Un
% TP = max(roundRobin * BWprb * BWeff .* log2(1 +SINR / SINReff));
% 
% AllTP(counter,:) = TP;  % user throughputs for the iterations

SINR_All(counter,:) = SINRdB;
end
 
% AllTP

% reverseAllTP = AllTP';
% 
% 
% figure(1)
% x = 1:numOfIterations;
% 
% for i = 1:numOfUsers
%     userTP = reverseAllTP(i,:);
%     plot(x,userTP);
%     hold on
%     
% end
% hold off
% xlabel('iterations')
% ylabel('throughput')
% 
% figure(2)
% 
% RAllTP = reshape(AllTP,[1,numOfUsers * numOfIterations]);
% cdfplot(RAllTP);
% %sortRALLTP = sort(RAllTP);
% 
% 
SINR_vector = reshape(SINR_All,[1,numOfUsers*numOfIterations]);
maxSINR = max(SINR_vector);
cdfplot(SINR_vector);
xlabel('SINR [dB]')
ylabel('CDF')
axis([-5 maxSINR 0 1])
title('SINR vs CDF')


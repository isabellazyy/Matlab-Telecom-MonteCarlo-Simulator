clear all, 
close all,
clc

%% Monte Carlo code - starting point of the thesis
%Here, you can change the number of iterations, cells and users,
% size of the simulation area can also be changed 
numOfIterations = 100;
numOfCells = 2;
numOfUsers = 3;
meter = 200; % meter

%% Open Square ABG model is used for path loss calculations
% Formula ==> PL(f,d)[dB] = 10 * a * log10(d) + B + 10 * y* log10(f) + SF

a = 4.14;
B = 3.66;
y = 2.43;
SF = 7;  % shadow fading
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
% 2x2 MIMO configuration
% TP = Nprb * BWprb * BWeff * log2(1 +SINR / SINReff)
BWeff = 0.42;
SINReff = 1.1;
Nprb = 50;
BWprb = 180000;

% criteria for the communication
Gmax = 90;
RSSmin = -100;
SINRmin = -10;
Smax = 7.67;

%% Random locations
% this section can be done in another M file
% Rectangular is chosen; e.g. width = 150m, height = 150m
% cell locations are chosen for only one time
% C =[C1x, C1y; C2x, C2y]
% U is the matrix for the locations of the users. User 
% locations are randomly chosen for each iteration
% U =[U1x, U1y; U2x, U2y; U3x, U3y]

% distance between Users and Cells
distance = [zeros(numOfCells,numOfUsers)];

% cell random locations 
C = randi([0 meter],numOfCells,2);

% All the SINR and TP results will be in these vectors at the end of
% simulations
SINR_All = zeros(numOfIterations,numOfUsers);  
TP_All = zeros(numOfIterations,numOfUsers);
Fairness = zeros(1,numOfIterations);

%% Monte Carlo Starts
for counter = 1:numOfIterations

    % here the user locations are created
    U = randi([0 meter],numOfUsers,2);
   
   for i = 1:numOfCells
      for k = 1:numOfUsers
          temp = C(i,:) - U(k,:);
          % distance between users and cells
          distance(i,k) = sqrt(sum((temp.^2)));   
      end 
   end
    
% Path Loss matrix ~ G(distance) 
G = zeros(size(distance));
G = 10 * a * log10(distance) + B + 10 * y* log10(f) + SF;

% Received Powers from the cells
RSSs = APL - G ;    % Received Signal Strengths (RSS) in dB
RSSmax = max(RSSs); % max RSS in dB

% Which cell serves which users
[x y] = size(RSSs);
servingCells= zeros(x,y); % this matrix can be called map between users and cells

for i = 1:y
    for k = 1:x
        if RSSmax(i) == RSSs(k,i)
            servingCells(k,i) = 1;
        else
            servingCells(k,i) = 0;
        end
    end
end

% if the max RSS of a user is lower than threshold, no communication
for i=1 : length(RSSmax)   
    if RSSmax(i) < RSSmin % dB
       RSSmax(i) = 0; 
    end    
end

% here is for calculations, no special meaning
G = reshape(G,[1, numOfUsers*numOfCells]);
temp_Cell = reshape(servingCells,[1,numOfUsers*numOfCells]);
 


% if the path loss is higher than threshold, neglect the channel
% In order to do this, we need to understand which cell has the best power
% for a user and then we need to see the path loss for this channel
for i = 1:numOfUsers * numOfCells
    if temp_Cell(i) == 1 && G(i) > Gmax
        temp_Cell(i) = 0;   % if path loss is higher than threshold, remove 
                            % the channel. It means there is no serving
                            % cell for the user
    end
end
% here is for calculations, no special meaning
G = reshape(G,[numOfCells, numOfUsers]);
servingCells = reshape(temp_Cell,[numOfCells,numOfUsers]);

% so, servingCells is changed (after considering the RSS and path loss
% rules, cell-user map changed. If conditions are not met, relations are
% removed


% converting to linear scale
RSSs_Linear = 10 .^(RSSs / 10);
RSSmax_Linear =  10.^(RSSmax/10);

% received interference in linear scale
Interference = sum(RSSs_Linear) - RSSmax_Linear;           
Noise = 10 ^ (NoisedBm/10);

temp_Cell = sum(servingCells);

%% SINR calculations
SINRdB = zeros(1,numOfUsers);
[row column] = size(servingCells);

% by using the changed servingCells, we can calculate the SINR
for i= 1:column
    
    if temp_Cell(i) == 1  % is there any cell serving the user
        SINR = RSSmax_Linear(i) / (Interference(i) + Noise);
        SINRdB(i) = 10 * log10(SINR);
        
        if SINRdB(i) < SINRmin   % If SINR is less than -10 dB,
                             % there is no communication
            SINRdB(i) = 0;
        end
    else
        SINRdB(i) = 0;
    end   
    %SINRdB(i)
end



% A cell may serve more than one user. Calculating how many users are
% served by one cell, then we can decide the scheduling
% Therefore, transpose of the servingcells is a need
layoutTrans = servingCells';
howManyUsers = sum(layoutTrans); % this vector gives us the number of users for each cell

temp = howManyUsers';       % transpose 
roundRobin = zeros(x,y);    % this matrix gives us the shares between users

for i = 1:y
    roundRobin(:,i) = servingCells(:,i) ./ temp;
end

roundRobin = roundRobin * Nprb; % how many PRB for each user
         

SINR_All(counter,:) = SINRdB;

[x y] = size(SINRdB);
     

%% Throughput calculations

for i=1:y
    A(i) = BWeff * log2(1 +(10 ^(SINRdB(i)/10)) / SINReff);
    
    if A(i) > Smax
        A(i) = Smax;
    end
    
end
TP = max(roundRobin * BWprb .* A);
TP_All(counter,:) = TP / (10^6);  % All user throughputs for the iterations


%% Jain's index
Fairness(counter) =  sum(TP)^2 / length(TP) / sum(TP.^2);



end
 
%% SINR-CDF plot
subplot(3,1,1)
 
SINR_vector = reshape(SINR_All,[1,numOfUsers*numOfIterations]);
maxSINR = max(SINR_vector);
cdfplot(SINR_vector);
xlabel('SINR [dB]')
ylabel('CDF')
axis([0 maxSINR 0 1])
title('SINR vs CDF')

%% Throughput-CDF plot
subplot(3,1,2)

TP_vector = reshape(TP_All,[1,numOfUsers*numOfIterations]);
maxTP = max(TP_vector);
cdfplot(TP_vector);
xlabel('TP [Mbps]')
ylabel('CDF')
% axis([0 maxTP 0 1]);
title('TP vs CDF');

%% Jain's index
subplot(3,1,3)

cdfplot(Fairness)
xlabel('Fairness')

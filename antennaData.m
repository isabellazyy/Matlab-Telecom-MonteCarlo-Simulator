% this function provides the received power
% between the cells and pixel
function [C] = antennaData()

    a = load('lossDataSmallCell2600MHz.mat');
    numOfAntenna = length(a.lossData);
    % C=zeros(numOfAntenna,13974);

    for i=1:numOfAntenna % 3
        C(i,:) = a.lossData(i).pathloss'; %(1:10);
    end
end
% 
% [x y] = size(C);
% 
% ff = FastFading(10,800e6);
% FF = zeros(size(C));
% 
% FF = ff(randi(size(ff),x,y));
% 
% ffData = C + FF;
% users = randi([1,1000],1,50);
% 
% received_power = min(C(:,users))
% 

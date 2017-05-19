function [RSSs] = losses(distance,APL,numOfPaths)

N = numOfPaths;                                    %Number of paths
% Open Square ABG model is used for path loss calculations
% Formula ==> PL(f,d)[dB] = 10 * a * log10(d) + B + 10 * y* log10(f) + SF

a = 4.14;
B = 3.66;
y = 2.43;
SF = 7;  % shadow fading
f = 0.8; % GHz

[numOfCells,numOfUsers] = size(distance);

M = zeros(size(distance));
X = zeros(1,N);
Y = zeros(1,N);

for k=1:numOfCells
   
    for j = 1:numOfUsers
        d = distance(k,j);
    
        for i=1:N
             
             G = 10 * a * log10(randi([floor(d+1) floor((d+1)*sqrt(2))],1,1)) + B + 10 * y* log10(f) + SF;
             P = APL - G;
             P = 10 ^ (P/10);
             angle = rand(1);
             x = P * cos(2 * pi * angle);
             y = P * sin(2 * pi * angle);
             X(i) = x;
             Y(i) = y;
             
        end
        
               
    M=abs(sum(X)+j*sum(Y));
    RSSs(k,j)=10*log10(M);
    
    end
end

end
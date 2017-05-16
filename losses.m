function [G] = losses(distance,APL,carrierFreq)

% In this function, shadow fading and path loss are considered

% Open Square ABG model is used for path loss calculations
% Formula ==> PL(f,d)[dB] = 10 * a * log10(d) + B + 10 * y* log10(f) + SF

a = 4.14;
B = 3.66;
y = 2.43;
SF = 7;  % shadow fading
f = carrierFreq/1e9; % GHz
        
G = 10 * a * log10(distance) + B + 10 * y* log10(f) + SF;

         
end
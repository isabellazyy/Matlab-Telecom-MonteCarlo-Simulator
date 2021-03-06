 
function [FF] = FastFading(numberOfPaths,carrierFreq)
% SUM OF SINUSOIDS FADING SIMULATORS
% fc - carrier frequency
% fd - Doppler frequency  
% fs - Sampling frequency
% ts - Sampling period
% N - Number of sinusoids

v=30 / 3600;			% vehicular speed	 
c=300 * 10 ^ 3;         % speed of light
fc=carrierFreq;         % carrier freq
fd=fc * v / c;          % doppler freq
fs=100000;             % sampling freq
ts=1 / fs;              
t=0 : ts : 1;           
N=numberOfPaths;                   % num of paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sum of sinusoids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=zeros(1,length(t));
y=zeros(1,length(t));

for n=1:N
    alpha=(rand-0.5)*2*pi;
    phi=(rand-0.5)*2*pi;       
    x=x+randn*cos(2*pi*fd*t*cos(alpha)+phi);
    y=y+randn*sin(2*pi*fd*t*cos(alpha)+phi);
end
z=(1/sqrt(N))*(x+1i*y);
r1=abs(z);
FF = 10*log10(r1);

plot(t,10*log10(r1))
title('Rayleigh')
 
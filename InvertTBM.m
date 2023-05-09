
function [betaanswer] = InvertTBM(theta,M,g,strength)
% programmer: David Willy, NAU
% Notes:
%   1. Simple function to find beta by inverting the theta-beta-mach equat.
%   2. Uses a guess to find the root. 
%   3. Strength - finds the weak or strong shock solution
%           a. if strenth =1, weak shocks (40)
%           b. if strenth =2, strong shocks (80)

f= @(beta,theta,M) theta - atand(2*cotd(beta)* ...
    (((M^2)*((sind(beta))^2))-1)/(((g+(cosd(2*beta)))*M^2)+2));

if strength==1
    betaguess=40;      % feel free to change if things get buggy
elseif strength==2
    betaguess=80;
end

options = optimset('Display','off');
betaanswer = fsolve(@(beta)f(beta,theta,M), betaguess, options);



function [T2, T02, P2, P02, Rho2, M2] = ...
    obliqueShocker(theta,beta,M,gamma,T01,P01,Rho01)
%obliqueShocker Calculates oblique shock properties
%   Inputs: Theta - Turning Angle
%           
Mnorm = M*sind(beta);
To_T = 1 + Mnorm.^2 .* ((gamma - 1) / 2);
rho_0_rho = (To_T).^(1 / (gamma - 1));
po_p = (To_T).^(gamma / (gamma - 1));

M2norm = sqrt((Mnorm^2+(2/(gamma-1)))/...
        ((2*gamma/(gamma-1))*Mnorm^2-1));
    
p2_p1 = 1 + (2*gamma/(gamma+1)) * (M2norm^2 - 1);

rho2_rho1 = ((gamma+1)*M2norm^2)/(2+(gamma-1)*M2norm^2);

T2_T1 = (1+(2*gamma/(gamma+1))*(M2norm^2-1))*((2+(gamma-1)*M2norm^2)/...
    ((gamma+1)*M2norm^2));

% diffs2s1 = Cp * log((1+(2*gam/(gam+1))*(Mnorm^2-1)) * ...
%     ((2+(gam-1)*Mnorm^2)/((gam+1)*Mnorm^2))) - ...
%     R * log(1+(2*gam/(gam+1))*(Mnorm^2-1));

% p02_p01 = exp(-1*diffs2s1/R);

p02_p1 = (((gamma+1)^2*M2norm^2)/(4*gamma*M2norm^2-2*(gamma-1)))^(gamma/ ...
    (gamma-1)) * ((1 - gamma + 2 * gamma * M2norm^2)/(gamma + 1));


T1 = To_T / T01;
T2 = T2_T1 * T1;
T02 = T01;

P1 = po_p / P01;
P2 = P1 * p2_p1;
P02 = P1 * p02_p1;

Rho1 = rho_0_rho / Rho01;
Rho2 = Rho1 * rho2_rho1;

M2 = M2norm / sind(beta - theta);

end
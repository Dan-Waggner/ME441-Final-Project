function [vectorCoordinates,vectorLength] = vectorGenerator(x1, x2, y1, y2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Generates a vector from 2 points' x- and y- values
%
%   Inputs:
%       x1: x-coordinate of point 1
%       x2: x-coordinate of point 2
%       y1: y-coordinate of point 1
%       y2: y-coordinate of point 2
%
%   Outputs:
%       vectorCoordinates:
%           xValue: x-value of the vector modeled from (0,0)
%           yValue: y-value of the vector modeled from (0,0)
%       vectorLength: The length of the vector between points 1 and 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Created by: Daniel Waggner
%   Date:       07 May 2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine x- and y- values
vectorCoordinates.xValue = x2 - x1;
vectorCoordinates.yValue = y2 - y1; 
% Determine Length of Vector
vectorLength = sqrt((vectorCoordinates.xValue)^2 + ...
    (vectorCoordinates.yValue)^2);
end
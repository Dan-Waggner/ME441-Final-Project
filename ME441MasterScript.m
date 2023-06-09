%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ME441-001 Final Project
%   11 May 2023
%   Daniel Waggner + Corbin Strycker
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial CMDs
clear; 
clc;
close all;
format short;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Givens %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Angle of Attack
alpha = 0;
% Requested Mach Number(s)
M1Init = 3.0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initial Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variable DFNs
%Flow Vector
flowVectorMag = 0.2;
flowVectorPosx = -1 * flowVectorMag * cosd(alpha);
flowVectorPosy = flowVectorMag * sind(alpha); 
%Coordinates
%  Top Coordinates
xvalstop = [flowVectorPosx, 0, 0.3, 0.7, 1.0];
yvalstop = [flowVectorPosy, 0, 0.1, 0.075, 0];
% Bottom Coordinates
xvalsbottom = [flowVectorPosx, 0, 1.0];
yvalsbottom = [flowVectorPosy, 0, 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Assumptions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. Acting Fluid is Air modeled as a calorically perfect gas
%   2. Frictionless Flow (Inviscid)
%   3. Fluid is Air @ STP
%   4. Flow is Adiabatic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Via Assumptions
gamma = 1.4;
R = 287; % J/kg-K
T_flow_abs = 20+273; % k
P_flow = 101325; % Pa
rho_flow = P_flow / (R * T_flow_abs); % kg/m^3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Various Calcs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Speed of Sound
a = sqrt(gamma * R * T_flow_abs);
% Flow Speed
v_flow = M1Init * a; % m/s
% Dynamic Pressure
P_dyn = 0.5 * rho_flow * v_flow^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Animation Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1
figure(1);
hold on;
[Flift, Fdrag] = ...
    airfoilPropertiesPlotter(P_flow,T_flow_abs,...
    M1Init, xvalstop, yvalstop, xvalsbottom, yvalsbottom, alpha,...
    rho_flow, gamma, flowVectorPosx, flowVectorPosy, R);
fprintf("Lift and Drag  %4.2f kN/m and %4.2f kN/m\n", ...
        Flift, Fdrag);
hold off;
pause(1);
%% Part 2
%Varying Alpha
alpha = -10:1:10; 
%Generating Animation Frames
% Working Directory DFN
workingDir = pwd;
mkdir(workingDir,'images');
% Generating Frames via For Loop
for k = 1:length(alpha)
    fig = figure("Resize","off",'visible','off');
%   Display Status of Generation
    fprintf(join(["Generating frame ", num2str(floor(k)),...
        " of ", num2str(floor(length(alpha))), "...\n"]));
%   Generate Graph
    hold on;
    flowVectorMag = 0.2;
    flowVectorPosx = -1 * flowVectorMag * cosd(alpha(k));
    flowVectorPosy = flowVectorMag * sind(alpha(k));
%   Plotting Airfoil Top/Bottom
    xvalstop = [flowVectorPosx, 0, 0.3, 0.7, 1.0];
    yvalstop = [flowVectorPosy, 0, 0.1, 0.075, 0];

    xvalsbottom = [flowVectorPosx, 0, 1.0];
    yvalsbottom = [flowVectorPosy, 0, 0];
%   Determining Properties and Graphing Shocks/Expansions
    [Flift, Fdrag] = ...
        airfoilPropertiesPlotter(P_flow,T_flow_abs,...
        M1Init, xvalstop, yvalstop, xvalsbottom,...
        yvalsbottom, alpha(k),rho_flow, gamma,...
        flowVectorPosx, flowVectorPosy, R);
%   Recording Lift/Drag forces
    FL(k) = Flift;
    FD(k) = Fdrag;
%   Pause for 1s to Enhance Performance (hopefully)
    pause(1);
%   Generate File Name
    saveName = join(["Animation Frame ",...
        num2str(floor(k)),".jpg"]);
%   Generate Exact Location to Save Image
    fullFileName = fullfile(workingDir, 'images', saveName);
%   Save Image
    saveas(fig,fullFileName);
%   Save File Name (in Order of Generation)
    imageNames(k) = saveName;
%   Finish with Current Frame
    hold off;
    close(fig);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Recording/Graphing Forces / Force Coefficients
%   Chord Length
    L = 1;
%   Width
    W = 1; % Unitz
%   Planform Area
    Ap = L*W; % m^2
% Coefficients
%   Lift
Cl = 2 * FL / (v_flow * Ap);
%   Drag
Cd = 2 * FD / (v_flow * Ap);
% Lift/Drag/Coefficients vs Alpha
tiledlayout(2,2);
figure(100); 
nexttile;
grid on;
hold on;
plot(alpha, FL/1000, 'DisplayName','Lift vs AoA');
title('Lift vs AoA');
xlabel('Angle of Attack, Degrees');
ylabel('Force of Lift, kN/m');
hold off;

nexttile;
grid on;
hold on
plot(alpha, FD/1000, 'DisplayName','Drag vs AoA');
title('Drag vs AoA');
xlabel('Angle of Attack, Degrees');
ylabel('Force of Drag, kN/m');
hold off;

nexttile;
grid on;
hold on;
plot(alpha, Cl);
title(['Coefficient of Lift',' vs AoA']);
xlabel('Angle of Attack, Degrees');
ylabel('Coefficient of Drag');
hold off;

nexttile;
grid on;
hold on;
plot(alpha, Cd);
title(['Coefficient of Drag', ' vs AoA']);
xlabel('Angle of Attack, Degrees');
ylabel('Coefficient of Drag');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Animations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating Animation
% User Notification
fprintf("Generating Animation...\n");
% Creating Output Video File
videoFileName = 'Airfoil_Animation.avi';
outputVideo = VideoWriter(fullfile( ...
    workingDir,videoFileName));
% Open Output Video File to Edit
open(outputVideo);
% Writing Images to Video
for ii = 1:length(imageNames)
%   User Notification that Images are Being Read
    fprintf(join(["Reading ", num2str(floor(ii)), ...
        " of ", num2str(floor(length(imageNames))),...
        " images...\n"]));
%   Read Images
   img = imread(fullfile(workingDir,'images',imageNames{ii}));
%   Send Images to File
   writeVideo(outputVideo,img);
%   Pause for 0.5s to (Hopefully) Enhance Performance
   pause(0.5);
end
%   Close Edited Video Fi,e
close(outputVideo);

%Generate Video to MatLab Video File
airfoilAnimation = VideoReader(fullfile(...
    workingDir,videoFileName));
%Assign Images to MatLab Video Frames
ii = 1;
while hasFrame(airfoilAnimation)
   mov(ii) = im2frame(readFrame(airfoilAnimation));
   ii = ii+1;
end
fprintf("Running Animation...\n");
pause(1);
%Play Video
figure;
%   Set Video Format
imshow(mov(1).cdata, 'Border', 'tight');
%   Play Animation
movie(mov, 20, 5);


%   salesartillery.com/fs/top-100-aerospace-companies
%   usajobs

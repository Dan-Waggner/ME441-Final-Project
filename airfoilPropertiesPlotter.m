function [Flift, Fdrag] = airfoilPropertiesPlotter(P_flow,T_flow_abs,...
    M1Init, xvalstop, yvalstop, xvalsbottom, yvalsbottom, alpha,...
    rho_flow, gamma, flowVectorPosx, flowVectorPosy, R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Plotting Airfoil
fill(xvalstop(2:length(xvalstop)),...
    yvalstop(2:length(yvalstop)), 'r');
xlims = [-0.4 1.4];
ylims = [-0.4, 0.4];
xlim(xlims);
ylim(ylims);
xticks(-0.2:0.2:1.2);

title("Asymmetric Airfoil Diagram");
xlabel("X-direction (m)");
ylabel("Y-direction (m)");
grid on; 
box on;

%Plotting Velocity Vector Line
plot([flowVectorPosx, 0], [flowVectorPosy, 0], 'b');

%Velocity Vector Labeling 
velStr1 = "AirFlow:";
velStr2 = append("M", num2str(round(M1Init,2)));
velStr3 = append("\alpha = ",...
    num2str(round(alpha, 2)));

%Incoming Flow Properties
[~, T, P, ~, ~] = flowisentropic(gamma, M1Init);
T1 = T * T_flow_abs;
Pstatic = P_flow * P;
Rho01 = rho_flow;

%Incoming Flow Labeling
propStr1 = join(["", num2str(round(Pstatic,2))," kPa"]);
propStr2 = join(["",num2str(round(T1,2))," K"]);
propStr3 = join(["",num2str(round(Rho01))," kg/m^3"]);

%Plotting/Labeling Regions
for j = 1:2
    if j == 1 % top
        for i = 2:length(xvalstop)
            plot([xvalstop(i),xvalstop(i)],...
                [yvalstop(i),ylims(2)],'g-');
            text(xvalstop(i-1)+0.1, ...
                ylims(2) - 0.1*ylims(2),...
                append([num2str(i-1)]));
        end
    elseif j == 2 % bottom
        for i = 2:length(xvalsbottom)
            plot([xvalsbottom(i),xvalsbottom(i)],...
                [yvalsbottom(i),ylims(1)],'g-');
            if i >= 3
                text(xvalsbottom(i-1)+0.2, ...
                    ylims(1) - 0.1*ylims(1),...
                    append([num2str(length(xvalstop) ...
                    + i-3)]));
            end
        end
    else
        fprintf(append(["How did this happen,\n",...
            "we're smarter than this!"]));
        break;
    end
end
text((xvalstop(length(xvalstop)) + ...
        0.1*xvalstop(length(xvalstop))),...
        ylims(2) - 0.25*ylims(2),...
        ["Slip", "Line", "Region"]);

%Initialized Variables
Flift = 0;
Fdrag = 0;
%Evaluating airfoil as 2 halves
%   j = 1: top half
%   j = 2: bottom half
    for j = 1:2
        %Resetting the Current Mach No. back to the
        %   oncoming flow's velocity
        %   + properties
        M1 = M1Init;
        P01 = P_flow;
        T01 = T_flow_abs;
    
        [~, T, P, ~, ~] = flowisentropic(gamma, M1);
        T1 = T * T01;
        Pstatic = P01 * P;
        Rho01 = rho_flow;
    
        %Handling Top vs Bottom Cases
        if j == 1
           valueInterestLength = length(xvalstop); 
        else
           valueInterestLength = length(xvalsbottom);
        end
        %Running Flow Analyses
        for i = 3:valueInterestLength
            %DFN relative vectors
            %   i.e. vectors that overlap the
            %   DFN's airfoil's edges
            if j == 1 % top half of foil points
                %x-values
                x1 = xvalstop(i-2);
                x2 = xvalstop(i-1);
                x3 = xvalstop(i);
                %y-values
                y1 = yvalstop(i-2);
                y2 = yvalstop(i-1);
                y3 = yvalstop(i);
            else % bottom half of foil points
                %x-values
                x1 = xvalsbottom(i-2);
                x2 = xvalsbottom(i-1);
                x3 = xvalsbottom(i);
                %y-values
                y1 = yvalsbottom(i-2);
                y2 = yvalsbottom(i-1);
                y3 = yvalsbottom(i);
            end
    
            %Handle Origin
            if x2==0 && y2==0
                linfit = polyfit([x1, x2], [y1, y2], 1);
    
                vec1coords = [x2, 0.2;...
                    y2, linfit(1)*0.2 + linfit(2)];
                vec2coords = [x2, x3;...
                    y2, y3];
                vec3coords = [0.2, x3;...
                    linfit(1)*0.2 + linfit(2), y3];
                entranceedge = [x2,x3;y2,y3];
                vecAxisCoords = [x2,0.2; y2,y2];
            else
                vec1coords = [x2, 2*x2; y2, 2*y2];
                vec2coords = [x2, x3; y2, y3];
                vec3coords = [2*x2, x3; 2*y2, y3]; 
                entranceedge = [x2,x3;y2,y3];
                vecAxisCoords = [x2,x3; y2,y2];
            end
            
            %Prepping Law of Cosines
            vec1length = sqrt((vec1coords(1,2) ...
                - vec1coords(1,1))^2 + (vec1coords(2,2)...
                - vec1coords(2,1))^2);
            
            vec2length = sqrt((vec2coords(1,2) ...
                - vec2coords(1,1))^2 + (vec2coords(2,2)...
                - vec2coords(2,1))^2);
        
            vec3length = sqrt((vec3coords(1,2) ...
                - vec3coords(1,1))^2 + (vec3coords(2,2)... 
                - vec3coords(2,1))^2);
        
            vec1 = [vec1coords(1,2) - vec1coords(1,1), ...
                vec1coords(2,2) - vec1coords(2,1)];
        
            vec2 = [vec2coords(1,2) - vec2coords(1,1), ...
                vec2coords(2,2) - vec2coords(2,1)];
    
            vecAxis = [vecAxisCoords(1,2) - vecAxisCoords(1,1),...
                vecAxisCoords(2,2) - vecAxisCoords(2,1)];
    
            entranceEdgevec = [entranceedge(1,2) ...
                - entranceedge(1,1), entranceedge(2,2) ...
                - entranceedge(2,1)];
            
            %Compute Turning Angles via Vectors
            %   Law of Cosines
            theta = acosd((vec1length^2 + vec2length^2 ...
                - vec3length^2) / (2 * vec1length ...
                * vec2length));     % degrees
    
            nu = alpha + rad2deg(AngleIn2D(vecAxis, ...
                entranceEdgevec));    % degrees
            %   Evaluating Direction
            %       (-) Clockwise:
            %           Expansion Fan
            %       (+) Counter-Clockwise:
            %           Oblique Shockwave
            
            if j == 2 % handling bottom side case
                        %   via flipping angle along 
                        %   x-axis
                thetaDirection = -1 * AngleIn2D(vec1,...
                    vec2);
                                            % radians
            else
                thetaDirection = AngleIn2D(vec1, vec2);
                                            % radians
            end
            
            %Handling Cases
            if thetaDirection < 0 % Expansion Fan
    
                [~, nu1, ~] = flowprandtlmeyer(...
                    gamma, M1, 'mach');
                nu2 = nu1 - theta;
                [M1, nu2, ~] = flowprandtlmeyer(...
                    gamma, abs(nu2), 'nu');
    
                [~, T, P, ~, ~] = ...
                    flowisentropic(gamma, M1);
                Tstatic = T01 * T;
                Pstatic = P01 * P;
    
                Flift = Pstatic * vec2length * cosd(nu) + Flift;
                Fdrag = Pstatic * vec2length * sind(nu) + Fdrag;
    
                dataStr = [join(["M ",...
                    num2str(round(M1,2))]),...
                    join(["",...
                    num2str(round(Pstatic/1000,2)),...
                    " kPa"]),...
                    join(["",...
                    num2str(round(Tstatic,2)),...
                    " K"])];
    
                %handling bottom half of foil
                %   via flipping angle along
                %   x-axis
                if j == 2
                    nu1 = nu1 - 90;
                    nu2 = nu2 - 90;
                end
                
                %Plotting Expansion Fans
                plot([x2, x2+.3*cosd(nu1)],[y2, ...
                    y2+.3*sind(nu1)], 'r');
                plot([x2, x2+.3*cosd(nu2)],[y2, ...
                    y2+.3*sind(nu2)], 'r');
    
                %Placing Data
                if j == 1
                    ypt = 0.75*ylims(2);
                    text(x2+0.1, ypt, dataStr); 
                else
                    ypt = 0.75*ylims(1);
                    text(x2+0.2, ypt, dataStr); 
                end           
        
            elseif thetaDirection > 0 % Oblique Shock
                beta = InvertTBM(theta, M1, gamma, 1);
    
                %Flow Conditions
                state.mach = M1;
                state.pressure = P01;
                state.temp = T1;
                state.rho = Rho01;
                
                shock.machNorm = M1*sind(beta);
                shock.theta = theta;
                shock.beta = beta;
    
                fluid.gam = gamma;
                fluid.R = R;

                %Oblique Shock Properties
                [state, shock, ~] = normal_shock(state,...
                    shock, fluid);
                Tstatic = state.temp(2);
                Pstatic = state.pressure(2);
                M1 = state.mach(2);
                
                %Lift and Drag
                Flift = Pstatic * vec2length * cosd(nu) + Flift;
                Fdrag = Pstatic * vec2length * sind(nu) + Fdrag;
                
                dataStr = [join(["M ",...
                    num2str(round(M1,2))]),...
                    join(["",...
                    num2str(round(Pstatic/1000,2)),...
                    " kPa"]),...
                    join(["",...
                    num2str(round(Tstatic,2)),...
                    " K"])];
    
                %Placing Data
                if j == 1
                    ypt = 0.75*ylims(2);
                    text(x2+0.1, ypt, dataStr);
                else
                    ypt = 0.75*ylims(1);
                    text(x2+0.2, ypt, dataStr);
                end
                
                if j == 2
                    beta = beta - 90;
                end
    
                %handling bottom half of foil
                %   via flipping angle along
                %   x-axis
                plot([x2, x2+.3*cosd(beta)], ...
                    [y2, y2+.3*sind(beta)], 'b--');
                
            else    
                %Case Where an Inputted Airfoil 
                %   Has Points Defining a Straight Edge
                continue;
            end
        end
    end
    liftStr = join(["F_L = ", num2str(round(abs(Flift/1000),2)), "kN/m"]);
    dragStr = join(["F_D = ", num2str(round(abs(Fdrag/1000),2)), "kN/m"]);
    text(-0.35, -0.2,...
    [velStr1, velStr2, velStr3, propStr1, propStr2, propStr3, liftStr,...
    dragStr]);
    Flift = abs(Flift/1000); % kN/m
    Fdrag = abs(Fdrag/1000); % kN/m
end
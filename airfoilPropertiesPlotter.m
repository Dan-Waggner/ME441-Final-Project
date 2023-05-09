function [Flift, Fdrag] = airfoilPropertiesPlotter(P_flow,T_flow_abs,...
    M1Init, xvalstop, yvalstop, xvalsbottom, yvalsbottom, alpha,...
    rho_flow, gamma, flowVectorPosx, flowVectorPosy, R)
%airfoilPropertiesPlotter Generates the plot of a given airfoil and places
%   all of its properties onto the graph

%%%%%%%%%%%%%%%%%%%%%%%%% Text Offset Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

xoffset = 0.025;
yoffset = 0.125;
lineLengthCoeff = .15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting Airfoil %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%% Plotting Airflow Vector %%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting Velocity Vector Line
plot([flowVectorPosx, 0], [flowVectorPosy, 0], 'b');
% Velocity Vector Labeling 
velStr2 = join(["M", num2str(round(M1Init,2))]);
velStr3 = join(["\alpha = ",...
    num2str(round(alpha, 2))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Airflow Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Incoming Flow Properties
[~, T, P, rhoi, ~] = flowisentropic(gamma, M1Init);
T1 = T * T_flow_abs;
Pstatic = P_flow * P;
Rho01 = rho_flow * rhoi;
% Incoming Flow Labeling
propStr1 = join(["P=", join([num2str(round(Pstatic,2)),"Pa"])]);
propStr2 = join(["T=",join([num2str(round(T1,2)),"K"])]);
propStr3 = join(["Rho=",join([num2str(round(Rho01,4)),"kg/m^3"])]);

%%%%%%%%%%%%%%%% Plotting Expansion Fans / Oblique Shocks %%%%%%%%%%%%%%%%%

% Plotting/Labeling Regions
for j = 1:2
    if j == 1 % top
        for i = 2:length(xvalstop)
            plot([xvalstop(i),xvalstop(i)],...
                [yvalstop(i),ylims(2)],'g-');

            % Handling Case for "Region 1"
            %   I.E. the Incoming Flow
            if i == 2
                text(-0.38, 0.0,...
                    join(["R",num2str(i-1)]));
            else
                text(xvalstop(i-1)+xoffset, ...
                    ylims(2)-0.0125,...
                    join(["R",num2str(i-1)]));
            end
            
        end
    elseif j == 2 % bottom
        for i = 2:length(xvalsbottom)
            plot([xvalsbottom(i),xvalsbottom(i)],...
                [yvalsbottom(i),ylims(1)],'g-');
            if i >= 3
                text(xvalsbottom(i-1)+4*xoffset, ...
                    ylims(1)+2*yoffset,...
                    join(["R",num2str(length(xvalstop) ...
                    + i-3)]));
            end
        end
    else
        fprintf(append(["How did this happen?\n",...
            "We're smarter than this!\n"]));
        break;
    end
end
text((xvalstop(length(xvalstop)) + ...
        8*xoffset),...
        0,...
        ["Slip", "Line", "Region"]);

% Initialized Variables
Flift = 0;
Fdrag = 0;
phi = linspace(-20,20,10);
Pstatics = zeros(2,length(phi));
thetaDirections = zeros(2,length(phi));
thetas = zeros(2,length(phi));
Machs = zeros(1,4);
Temps = zeros(1,4);
TempStatics = zeros(1,4);
% Evaluating airfoil as 2 halves
%   j = 1: top half
%   j = 2: top slip
%   j = 3: bottom half
%   j = 4: bottom slip
    for j = 1:4
        %Resetting the Current Mach No. back to the
        %   oncoming flow's velocity
        %   + properties
        if j == 1 || j == 3
            M1(1) = M1Init;
            P01 = P_flow;
            T01 = T_flow_abs;
        
            [~, T, P, ~, ~] = flowisentropic(gamma, M1(1));
            T1 = T * T01;
            Pstatic = P01 * P;
            Rho01 = rho_flow * rhoi; 
        end
         
        % Handling Top vs Bottom Cases
        if j == 1
           valueInterestLength = length(xvalstop); 
        elseif j ==3
           valueInterestLength = length(xvalsbottom);
        else
            valueInterestLength = 3;
        end
        % Running Flow Analyses
        for i = 3:valueInterestLength
            % DFN relative vectors
            %   i.e. vectors that overlap the
            %   DFN's airfoil's edges
            if j == 1 % top half of foil points
                % x-values
                x1 = xvalstop(i-2);
                x2 = xvalstop(i-1);
                x3 = xvalstop(i);
                % y-values
                y1 = yvalstop(i-2);
                y2 = yvalstop(i-1);
                y3 = yvalstop(i);
                % x-plot coordinate for txt
                xplotval = x2+xoffset;
                % y-plot coordinate for txt
                yplotval = ylims(2) - yoffset;
            elseif j == 3 % bottom half of foil points
                % x-values
                x1 = xvalsbottom(i-2);
                x2 = xvalsbottom(i-1);
                x3 = xvalsbottom(i);
                % y-values
                y1 = yvalsbottom(i-2);
                y2 = yvalsbottom(i-1);
                y3 = yvalsbottom(i);
                % x-plot coordinate for txt
                xplotval = x2+4*xoffset;
                % y-plot coordinate for txt
                yplotval = ylims(1)+yoffset;
            elseif j == 2   % top slip
                % x-values
                x1 = xvalstop(i-1);
                x2 = xvalstop(i);
                x3 = x2 + xoffset;
                % y-values
                y1 = yvalstop(i-1);
                y2 = yvalstop(i);
                y3 = y2 + yoffset;
                % x-plot coordinate for txt
                xplotval = x2+4*xoffset;
                % y-plot coordinate for txt
                yplotval = ylims(2)+yoffset;
            else    % bottom slip
                % x-values
                x1 = xvalsbottom(i-1);
                x2 = xvalsbottom(i);
                x3 = x2 + xoffset;
                % y-values
                y1 = yvalsbottom(i-1);
                y2 = yvalsbottom(i);
                y3 = y2 + yoffset;
                % x-plot coordinate for txt
                xplotval = x2+4*xoffset;
                % y-plot coordinate for txt
                yplotval = ylims(1)+yoffset;
            end
    
            % Prepping Law of Cosines
            if x2==0 && y2==0 % Handle Origin
                linfit = polyfit([x1, x2], [y1, y2], 1);
    
                [vec1coords, vec1length] = vectorGenerator(x2, 0.2, y2,...
                    linfit(1)*0.2 + linfit(2));
                vec1 = [vec1coords.xValue, vec1coords.yValue];
                
                [vec2coords, vec2length] = vectorGenerator(x2, x3, y2, y3);
                vec2 = [vec2coords.xValue, vec2coords.yValue];

                [~, vec3length] = vectorGenerator(0.2, x3,...
                    linfit(1)*0.2 + linfit(2), y3);

                [entranceedge, ~] = vectorGenerator(x2,x3,y2,y3);
                entranceEdgevec = [entranceedge.xValue,...
                    entranceedge.yValue];

                [vecAxisCoords, ~] = vectorGenerator(x2,0.2,y2,y2);
                vecAxis = [vecAxisCoords.xValue, vecAxisCoords.yValue];

            else
                [vec1coords, vec1length] = vectorGenerator(x2, 2*x2, ...
                    y2, 2*y2);
                vec1 = [vec1coords.xValue, vec1coords.yValue];

                [vec2coords, vec2length] = vectorGenerator(x2, x3, y2, y3);
                vec2 = [vec2coords.xValue, vec2coords.yValue];

                [~, vec3length] = vectorGenerator(2*x2, x3, ...
                    2*y2, y3);

                [entranceedge, ~] = vectorGenerator(x2,x3,y2,y3);
                entranceEdgevec = [entranceedge.xValue,...
                    entranceedge.yValue];

                [vecAxisCoords, ~] = vectorGenerator(x2,x3,y2,y2);
                vecAxis = [vecAxisCoords.xValue, vecAxisCoords.yValue];
                
            end
            
            % Compute Turning Angles via Vectors
            %   Law of Cosines
            theta = acosd((vec1length^2 + vec2length^2 ...
                - vec3length^2) / (2 * vec1length ...
                * vec2length));     % degrees
            if j == 2
                theta = theta + phi;
                thetas(1,:) = theta;
            elseif j == 4
                theta = theta - phi;
                thetas(2,:) = theta;
            end
            
            %   Angle for Finding Lift/Drag Forces
            nu = alpha + rad2deg(AngleIn2D(vecAxis, ...
                entranceEdgevec));  % degrees
            % Evaluating Direction
            %       (-) Clockwise:
            %           Expansion Fan
            %       (+) Counter-Clockwise:
            %           Oblique Shockwave
            
            if j == 3 || j == 4   % handling bottom side case
                        %   via flipping angle along 
                        %   x-axis
                thetaDirection = rad2deg(-1 * AngleIn2D(vec1,...
                    vec2));
                if j == 4
                thetaDirection = (thetaDirection + -phi);
                end

                thetaDirections(2,:) = -thetaDirection;
                                            % degrees
            else
                thetaDirection = rad2deg(AngleIn2D(vec1, vec2));
                                            % degrees
                if j == 2
                thetaDirection = thetaDirection + phi;
                end

                thetaDirections(1,:) = thetaDirection;
                
            end
            

            % Handling Cases
            for k = 1:length(theta)
                if thetaDirection(k) < 0 % Expansion Fan
                    % Evaluating Properties
                    %   Before
                    if M1 < 1
                        continue
                    else
                        [~, nu1, ~] = flowprandtlmeyer(...
                        gamma, M1, 'mach');
                        nu2 = nu1 - abs(theta(k));
                        [M1, nu2, ~] = flowprandtlmeyer(...
                            gamma, abs(nu2), 'nu');
                        %   After
                        [~, T, P, rho, ~] = ...
                            flowisentropic(gamma, M1);
                        Tstatic = T01 * T;
                        Pstatic = P01 * P;
                        T1 = Tstatic;
                        Rho01 = Rho01 * rho;
                    end
                    
                    % Lift/Drag Forces
                    if j == 1 || j == 3
                        Flift = Pstatic * vec2length * cosd(nu) + Flift;
                        Fdrag = Pstatic * vec2length * sind(nu) + Fdrag;
                    end
        
                    % Prepping Data for Plotting
                    dataStr = [join(["M=",...
                        num2str(round(M1,2))]),...
                        append([join(["P=",...
                        num2str(round(Pstatic))]), "Pa"]),...
                        append([join(["T=",...
                        num2str(round(Tstatic,2))]), "Kelvin"])];
        
                    % handling bottom half of foil
                    %   via flipping angle along
                    %   x-axis
                    thetalog = theta(k);
                    if j == 3 || j == 4
                        nu1 = nu1 - 90;
                        nu2 = nu2 - 90;
                        thetalog = -theta(k);
                    end
                    thetas(j) = thetalog;
                    
                    % Get all static Pressures for varying thetas
                    if j == 2 || j == 4
                        if j == 2
                            m = 1;
                        elseif j == 4
                            m = 2;
                        end
                        Pstatics(m,k) = Pstatic;
                    end
                    
    
    %                 if i == valueInterestLength
    %                     Pend(j) = Pstatic;
    %                     PendStag(j) = P01;
    %                     thetaEnd(j) = theta(k);
    %                     Mend(j) = M1;
    %                     Tend(j) = Tstatic;
    %                     rhoEnd(j) = Rho01;
    %                 end
                    if j == 1 || j == 3
                        % Plotting Expansion Fans
                        plot([x2, x2+lineLengthCoeff*cosd(nu1)],[y2, ...
                            y2+lineLengthCoeff*sind(nu1)], 'r');
                        plot([x2, x2+lineLengthCoeff*cosd(nu2)],[y2, ...
                            y2+lineLengthCoeff*sind(nu2)], 'r');
            
                        % Placing Data onto Plot     
                        text(xplotval, yplotval, dataStr); 
                    end

                    Machs(j) = M1;
                    Temps(j) = T01(1);
                    TempStatics(j) = T1(1);

                elseif thetaDirection(k) > 0 % Oblique Shock
                    beta = InvertTBM(theta(k), M1, gamma, 1);
        
                    % Flow Conditions
                    state.mach = M1;
                    state.pressure = Pstatic;
                    state.temp = T1;
                    state.rho = Rho01;
                    
                    shock.machNorm = M1*sind(beta);
                    shock.theta = theta(k);
                    shock.beta = beta;
        
                    fluid.gam = gamma;
                    fluid.R = R;
    
                    % Oblique Shock Properties
                    [state, shock, ~] = oblique_shock(state,...
                        shock, fluid);
                    Tstatic = state.temp(2);
                    Pstatic = state.pressure(2);
                    M1 = state.mach(2);
                    Rho01 = state.rho(2);
                    P01 = state.stagPressure(2);
                    T1 = Tstatic;
    
                    % Lift/Drag Forces
                    if j == 1 || j == 3
                        Flift = Pstatic * vec2length * cosd(nu) + Flift;
                        Fdrag = Pstatic * vec2length * sind(nu) + Fdrag;
                    end
                    
                    % Prepping Data for Plotting
                    dataStr = [join(["M=",...
                        num2str(round(M1,2))]),...
                        append([join(["P=",...
                        num2str(round(Pstatic))]), "Pa"]),...
                        append([join(["T=",...
                        num2str(round(Tstatic,2))]), "Kelvin"])];
                    
                    % handling bottom half of foil
                    %   via flipping angle along
                    %   x-axis
                    thetalog = theta(k);
                    if j == 3 || j == 4
                        beta = beta - 90;
                        thetalog = -theta(k); 
                    end
                    thetas(j) = thetalog;
    
                    % Get all static Pressures for varying thetas
                    if j == 2 || j == 4
                        if j == 2
                            m = 1;
                        elseif j == 4
                            m = 2;
                        end
                        Pstatics(m,k) = Pstatic;
                    end
    
                    % Plotting
                    if j == 1 || j == 3
                        % Placing Data onto Plot
                        text(xplotval, yplotval, dataStr); 
        
                        % Plotting Oblique Shockwave
                        plot([x2, x2+lineLengthCoeff*cosd(beta)], ...
                            [y2, y2+lineLengthCoeff*sind(beta)], 'b--');
                    else
                        continue;
                    end
                    
                    Machs(j) = M1(1);
                    Temps(j) = T01(1);
                    TempStatics(j) = T1(1);
                    
                else    
                    % Case Where an Inputted Airfoil 
                    %   Has Points Defining a Straight Edge
                    
                    % Since no new information is computed, prior info
                    %   must be stored
                    Machs(j) = M1(1);
                    Temps(j) = T01(1);
                    TempStatics(j) = T1(1);
                    thetas(j) = 0;

                    % Generating Data String
                    dataStr = [join(["M=",...
                        num2str(round(M1,2))]),...
                        append([join(["P=",...
                        num2str(round(Pstatic))]), "Pa"]),...
                        append([join(["T=",...
                        num2str(round(T1(1),2))]), "Kelvin"])];

                    % Placing Data onto Plot
                        text(xplotval, yplotval, dataStr); 
                    continue;
                end
            end
        end
    end

%%%%%%%%%%%%%%%%%%%% Placing Incoming Flow Properties %%%%%%%%%%%%%%%%%%%%%

    % Prepping String to Display Lifting Force
    liftStr = join(["F_L=", num2str(round(abs(Flift/1000),2)), "kN/m"]);
    % Prepping String to Display Drag force
    dragStr = join(["F_D=", num2str(round(abs(Fdrag/1000),2)), "kN/m"]);
    % Placing Lift/Drag Forces onto Plot
    text(-0.38, -0.18,...
    [velStr2, velStr3, propStr1, propStr2, propStr3, liftStr,...
    dragStr]);
    % Converting Lift/Drag into kN per Meter
    Flift = abs(Flift/1000); % kN/m
    Fdrag = abs(Fdrag/1000); % kN/m

%%%%%%%%%%%%%%%%%%%%%%%% Handling Slip Line Region %%%%%%%%%%%%%%%%%%%%%%%%

    try
    % Finding Correct Phi
        [phi,Ps] = intersections(phi, Pstatics(1,:), phi, Pstatics(2,:));
        [~, pos] = min(abs(phi));
        SlipLineStaticP = Ps(pos);
        SlipLineAngle = phi(pos);
    
    % Prepping for Stat Calcs
        thetaDirection = [thetaDirections(1,pos),thetaDirections(2,pos)];
        theta = thetaDirection; 
    
        M1 = [Machs(1), Machs(3)];
        T01 = [Temps(1), Temps(3)];
        T1 = [TempStatics(1),TempStatics(3)];
    
        % top x-values
        x1(1) = xvalstop(length(xvalstop)-1);
        x2(1) = xvalstop(length(xvalstop));
        x3(1) = x2(1) + lineLengthCoeff*cosd(SlipLineAngle);
        % top y-values
        y1(1) = yvalstop(length(yvalstop)-1);
        y2(1) = yvalstop(length(yvalstop));
        y3(1) = y2(1) + lineLengthCoeff*sind(SlipLineAngle);
        % bottom x-values
        x1(2) = xvalsbottom(length(xvalsbottom)-1);
        x2(2) = x2(1);
        x3(2) = x3(1);
        % bottom y-values
        y1(2) = xvalsbottom(length(yvalsbottom)-1);
        y2(2) = y2(1);
        y3(2) = y3(1);
    
    % Plotting Slip Line
        plot([x2(2),x3(2)], [y2(2),y3(2)], 'g--');
        
        for j = 1:2
    
            if j == 1   % top slip
                % x-plot coordinate for txt
                xplotval = x2(j)+xoffset;
                % y-plot coordinate for txt
                yplotval = ylims(2) - yoffset;
            else    % bottom slip
                % x-plot coordinate for txt
                xplotval = x2(j)+xoffset;
                % y-plot coordinate for txt
                yplotval = ylims(1)+yoffset;
            end
    
            if thetaDirection(j) < 0 % Expansion Fan
                % Evaluating Properties
                %   Before
                [~, nu1, ~] = flowprandtlmeyer(...
                    gamma, M1(j), 'mach');
                nu2 = nu1 - abs(theta(j));
                [M2, nu2, ~] = flowprandtlmeyer(...
                    gamma, abs(nu2), 'nu');
                %   After
                [~, T, ~, ~, ~] = ...
                    flowisentropic(gamma, M2);
                Tstatic = T01(j) * T;
                Pstatic = SlipLineStaticP;
                % Lift/Drag Forces
                Flift = Pstatic * vec2length * cosd(nu) + Flift;
                Fdrag = Pstatic * vec2length * sind(nu) + Fdrag;
                
    
                % Prepping Data for Plotting
                dataStr = [join(["M=",...
                    num2str(round(M2,2))]),...
                    append([join(["P=",...
                    num2str(round(Pstatic))]), "Pa"]),...
                    append([join(["T=",...
                    num2str(round(Tstatic,2))]), "Kelvin"])];
    
                % handling bottom half of foil
                %   via flipping angle along
                %   x-axis
                if j == 2
                    nu1 = nu1 - 90;
                    nu2 = nu2 - 90;
                end
    
                % Plotting Expansion Fans
                plot([x2(j), x2(j)+lineLengthCoeff*cosd(nu1)],[y2(j), ...
                    y2(j)+lineLengthCoeff*sind(nu1)], 'r');
                plot([x2(j), x2(j)+lineLengthCoeff*cosd(nu2)],[y2(j), ...
                    y2(j)+lineLengthCoeff*sind(nu2)], 'r');
    
                % Placing Data onto Plot     
                text(xplotval, yplotval, dataStr); 
        
            elseif thetaDirection(j) > 0 % Oblique Shock
                beta = InvertTBM(theta(j), M1(j), gamma, 1);
    
                % Flow Conditions
                Pstatic = SlipLineStaticP;
                
                state.mach = M1(j);
                state.pressure = Pstatic;
                state.temp = T1(j);
                state.rho = Rho01;
                
                shock.machNorm = M1(j)*sind(beta);
                shock.theta = theta(j);
                shock.beta = beta;
    
                fluid.gam = gamma;
                fluid.R = R;
    
                % Oblique Shock Properties
                [state, shock, ~] = oblique_shock(state,...
                    shock, fluid);
                Tstatic = state.temp(2);
                Pstatic = SlipLineStaticP;
                M2 = state.mach(2);
                Rho01 = state.rho(2);
    
                % Lift/Drag Forces
                Flift = Pstatic * vec2length * cosd(nu) + Flift;
                Fdrag = Pstatic * vec2length * sind(nu) + Fdrag;
                
                % Prepping Data for Plotting
                dataStr = [join(["M=",...
                    num2str(round(M2,2))]),...
                    append([join(["P=",...
                    num2str(round(Pstatic))]), "Pa"]),...
                    append([join(["T=",...
                    num2str(round(Tstatic,2))]), "Kelvin"])];
                
                % handling bottom half of foil
                %   via flipping angle along
                %   x-axis
                if j == 2
                    beta = beta - 90;
                end
    
                % Plotting
                    % Placing Data onto Plot
                    text(xplotval, yplotval, dataStr); 
    
                    % Plotting Oblique Shockwave
                    plot([x2(j), x2(j)+lineLengthCoeff*cosd(beta)], ...
                        [y2(j), y2(j)+lineLengthCoeff*sind(beta)], 'b--');
                
            else    
                % Case Where an Inputted Airfoil 
                %   Has Points Defining a Straight Edge
                continue;
            end
            
        end
    catch
        text(xplotval, yplotval, ["No Valid ", "Slip Line Found"]); 
        
    end

end
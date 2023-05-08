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
propStr3 = join(["Rho=",join([num2str(round(Rho01)),"kg/m^3"])]);

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
        xoffset),...
        ylims(2) - yoffset,...
        ["Slip", "Line", "Region"]);

% Initialized Variables
Flift = 0;
Fdrag = 0;
% Evaluating airfoil as 2 halves
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
        Rho01 = rho_flow * rhoi;  
        % Handling Top vs Bottom Cases
        if j == 1
           valueInterestLength = length(xvalstop); 
        else
           valueInterestLength = length(xvalsbottom);
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
            else % bottom half of foil points
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
            
            %   Angle for Finding Lift/Drag Forces
            nu = alpha + rad2deg(AngleIn2D(vecAxis, ...
                entranceEdgevec));  % degrees
            % Evaluating Direction
            %       (-) Clockwise:
            %           Expansion Fan
            %       (+) Counter-Clockwise:
            %           Oblique Shockwave
            
            if j == 2   % handling bottom side case
                        %   via flipping angle along 
                        %   x-axis
                thetaDirection = -1 * AngleIn2D(vec1,...
                    vec2);
                                            % radians
            else
                thetaDirection = AngleIn2D(vec1, vec2);
                                            % radians
            end
            
            % Handling Cases
            if thetaDirection < 0 % Expansion Fan
                % Evaluating Properties
                %   Before
                if M1 < 1
                    continue;
                end
                [~, nu1, ~] = flowprandtlmeyer(...
                    gamma, M1, 'mach');
                nu2 = nu1 - theta;
                [M1, nu2, ~] = flowprandtlmeyer(...
                    gamma, abs(nu2), 'nu');
                %   After
                [~, T, P, rho, ~] = ...
                    flowisentropic(gamma, M1);
                Tstatic = T01 * T;
                Pstatic = P01 / P;
                Rho01 = Rho01 * rho;
                % Lift/Drag Forces
                Flift = Pstatic * vec2length * cosd(nu) + Flift;
                Fdrag = Pstatic * vec2length * sind(nu) + Fdrag;
    
                % Prepping Data for Plotting
                dataStr = [join(["M=",...
                    num2str(round(M1,2))]),...
                    append([join(["P=",...
                    num2str(round(Pstatic,2))]), "Pa"]),...
                    append([join(["T=",...
                    num2str(round(Tstatic,2))]), "Kelvin"])];
    
                % handling bottom half of foil
                %   via flipping angle along
                %   x-axis
                if j == 2
                    nu1 = nu1 - 90;
                    nu2 = nu2 - 90;
                end

                % Get last region's pressure value
                if i == valueInterestLength
                    Pend(j) = Pstatic;
                    PendStag(j) = P01;
                    thetaEnd(j) = theta;
                    Mend(j) = M1;
                    Tend(j) = Tstatic;
                    rhoEnd(j) = Rho01;
                end
                
                % Plotting Expansion Fans
                plot([x2, x2+lineLengthCoeff*cosd(nu1)],[y2, ...
                    y2+lineLengthCoeff*sind(nu1)], 'r');
                plot([x2, x2+lineLengthCoeff*cosd(nu2)],[y2, ...
                    y2+lineLengthCoeff*sind(nu2)], 'r');
    
                % Placing Data onto Plot     
                text(xplotval, yplotval, dataStr); 
        
            elseif thetaDirection > 0 % Oblique Shock
                beta = InvertTBM(theta, M1, gamma, 1);
    
                % Flow Conditions
                state.mach = M1;
                state.pressure = Pstatic;
                state.temp = T1;
                state.rho = Rho01;
                
                shock.machNorm = M1*sind(beta);
                shock.theta = theta;
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

                % Lift and Drag
                Flift = Pstatic * vec2length * cosd(nu) + Flift;
                Fdrag = Pstatic * vec2length * sind(nu) + Fdrag;
                
                % Prepping Data for Plotting
                dataStr = [join(["M=",...
                    num2str(round(M1,2))]),...
                    append([join(["P=",...
                    num2str(round(Pstatic,2))]), "Pa"]),...
                    append([join(["T=",...
                    num2str(round(Tstatic,2))]), "Kelvin"])];
                
                % handling bottom half of foil
                %   via flipping angle along
                %   x-axis
                if j == 2
                    beta = beta - 90;
                end

                % Get last region's pressure value
                if i == valueInterestLength
                    Pend(j) = Pstatic;
                    PendStag(j) = P01;
                    thetaEnd(j) = theta;
                    Mend(j) = M1;
                    Tend(j) = Tstatic;
                    rhoEnd(j) = Rho01;
                end

                % Placing Data onto Plot     
                text(xplotval, yplotval, dataStr); 
    
                % Plotting Oblique Shockwave
                plot([x2, x2+lineLengthCoeff*cosd(beta)], ...
                    [y2, y2+lineLengthCoeff*sind(beta)], 'b--');
                
            else    
                % Case Where an Inputted Airfoil 
                %   Has Points Defining a Straight Edge
                continue;
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
%     % Phi Values
%     phi = linspace(-90,90,180); % degrees
%     
%     % last 2 top x-values
%     x1 = xvalstop(length(xvalstop)-1);
%     x2 = xvalstop(length(xvalstop));
% 
%     % last bottom x-value (not end pt)
%     x3 = xvalsbottom(length(xvalsbottom)-1);
%     
%     % last 2 top y-values
%     y1 = yvalstop(length(yvalstop)-1);
%     y2 = yvalstop(length(yvalstop));
% 
%     % last bottom y-value (not end pt)
%     y3 = yvalsbottom(length(yvalsbottom)-1);
%     
%     % Top line vector
%     [vtopendcoords, vtopendlength] = vectorGenerator(x1,x2,y1,y2);
%     vtopend = [vtopendcoords.xValue, vtopendcoords.yValue];
% 
%     % Extension of top line vector
%     linfittopext = polyfit([x1, x2], [y1, y2], 1);
%     [vtopendcoordsext, vtopendlengthext] = ...
%         vectorGenerator(x2, 0.2+x2, y2,...
%                     linfittopext(1)*(0.2+x2) + linfittopext(2));
%     vtopendext = [vtopendcoordsext.xValue, vtopendcoordsext.yValue];
% 
%     % Bottom line vector
%     [vbottomendcoords, vbottomendlength] = vectorGenerator(x3,x2,y3,y2);
%     vbottomend = [vbottomendcoords.xValue, ...
%         vbottomendcoords.yValue];
% 
%     % Extension of bottom line vector
%     linfitbottomext = polyfit([x3, x2], [y3, y2], 1);
%     [vbottomendcoordsext, vbottomendlengthext] = ...
%         vectorGenerator(x2, 0.2+x2, y2,...
%                     linfitbottomext(1)*(0.2+x2) + linfitbottomext(2));
%     vbottomendext = [vbottomendcoordsext.xValue,...
%         vbottomendcoordsext.yValue];
% 
%     for i = 1:length(phi)
%         
%         % x-value of imaginary slip line
%         x4 = x2 + lineLengthCoeff;
%     
%         % y-value of imaginary slip line
%         y4 = y2+lineLengthCoeff*sind(phi(i));
% 
%         % imaginary slip line vector
%         [sliplinecoords, sliplinelength] = vectorGenerator(x3,x4,y3,y4);
%         vslipline = [sliplinecoords.xValue, sliplinecoords.yValue];
%     
%         % "side c" for law of cosines - top
%         [sidectopcoords, sidectoplength] = ...
%             vectorGenerator(vtopendext(1), vslipline(1),...
%             vtopendext(2), vslipline(2));
%         vsidectop = [sidectopcoords.xValue, sidectopcoords.yValue];
%     
%         % "side c" for law of cosines - bottom
%         [sidecbottomcoords, sidecbottomlength] = ...
%             vectorGenerator(vbottomendext(1), vslipline(1),...
%             vbottomendext(2), vslipline(2));
%         vsidecbottom = [sidecbottomcoords.xValue,...
%             sidecbottomcoords.yValue];
%     
%         % Vectors from
%         
%         % Compute Top Turning Angle via Vectors
%         %   Law of Cosines
%         thetaTop = acosd((vtopendlengthext^2 + sliplinelength^2 ...
%             - sidectoplength^2) / (2 * vtopendlengthext ...
%             * sliplinelength));     % degrees
% 
%         thetaDirectionTop = AngleIn2D(vtopendext,...
%                     vslipline);
%     
%         % Compute Top Turning Angle via Vectors
%         %   Law of Cosines
%         thetaBottom = acosd((vbottomendlengthext^2 + sliplinelength^2 ...
%             - sidecbottomlength^2) / (2 * vbottomendlengthext ...
%             * sliplinelength));     % degrees
% 
%         thetaDirectionBottom = -AngleIn2D(vbottomendext,...
%                     vslipline);
%         
%         if thetaDirectionTop < 0
%             [~, nu1, ~] = flowprandtlmeyer(...
%                 gamma, Mend(1), 'mach');
%             nu2 = nu1 - thetaTop;
%             [Mslip(1), nu2, ~] = flowprandtlmeyer(...
%                 gamma, abs(nu2), 'nu');
%             %   After
%             [~, ~, P, rho, ~] = ...
%                 flowisentropic(gamma, Mslip(1));
%             PstaticTop = PendStag(1) / P;
%             Rho01 = Rho01 * rho;
%             
%         else
%             betaEnd(1) = InvertTBM(thetaTop, Mend(1), gamma, 1);
%     
%             % Flow Conditions
%             state.mach = Mend(1);
%             state.pressure = Pend(1);
%             state.temp = Tend(1);
%             state.rho = rhoEnd(1);
%             
%             shock.machNorm = Mend(1)*sind(betaEnd(1));
%             shock.theta = thetaTop;
%             shock.beta = betaEnd(1);
% 
%             fluid.gam = gamma;
%             fluid.R = R;
% 
%             % Oblique Shock Properties
%             [state, shock, ~] = oblique_shock(state,...
%                 shock, fluid);
%             PstaticTop = state.pressure(2);
%             Rho01 = state.rho(2);
%             Mslip(1) = state.mach(2);
%         end
% 
%         if thetaDirectionBottom < 0
%             [~, nu1, ~] = flowprandtlmeyer(...
%                 gamma, Mend(1), 'mach');
%             nu2 = nu1 - thetaBottom;
%             [Mslip(2), nu2, ~] = flowprandtlmeyer(...
%                 gamma, abs(nu2), 'nu');
%             %   After
%             [~, ~, P, rho, ~] = ...
%                 flowisentropic(gamma, Mslip(2));
%             PstaticBottom = PendStag(2) / P;
%             Rho01 = Rho01 * rho;
%             
%         else
%             betaEnd(2) = InvertTBM(thetaBottom, Mend(2), gamma, 1);
%     
%             % Flow Conditions
%             state.mach = Mend(2);
%             state.pressure = Pend(2);
%             state.temp = Tend(2);
%             state.rho = rhoEnd(2);
%             
%             shock.machNorm = Mend(2)*sind(betaEnd(2));
%             shock.theta = thetaBottom;
%             shock.beta = betaEnd(2);
% 
%             fluid.gam = gamma;
%             fluid.R = R;
% 
%             % Oblique Shock Properties
%             [state, shock, ~] = oblique_shock(state,...
%                 shock, fluid);
%             PstaticBottom = state.pressure(2);
%             Rho01 = state.rho(2);
%         end
%         
%         PStaticTopCoords(i) = PstaticTop;
%         PStaticBottomCoords(i) = PstaticBottom;
% 
%     end
%     
%     [PStaticSlipline, sliplineAngle] = ...
%         intersections(PStaticTopCoords, phi, PStaticBottomCoords, phi);
% 
%     phi = sliplineAngle(1);
% 
%     % x-value of imaginary slip line
%     x4 = x2 + lineLengthCoeff;
% 
%     % y-value of imaginary slip line
%     y4 = y2+lineLengthCoeff*sind(phi(1));
% 
%     plot([x2,x4],[y2,y4], 'm--');
end
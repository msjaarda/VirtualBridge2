function [Max,pd,x_values,y_values] = qInvestInitial(BM,ClassType,DistTypes,MaxEvents,ILData)
%qINVESTINITIAL Steps 2, 3, and 4 of former qInvestigation

% --- Step 3: Build Structure with Block Maxima ---
% Convert ClassT to m (number) form m = 1 is All, m = 2 is ClassOW, m = 3 is Class
% Reminder, m = 1 is ClassT 'All', m = 2 is 'ClassOW', and m = 3 is 'Class'

% For each Influence case
for r = 1:length(ILData)
    
    % Reset MaxEvents and select r
    MaxEventsSub = MaxEvents(MaxEvents.InfCase == r,:);
    % Transform AxTandem into Array (necessary for splitapply)
    Z = MaxEventsSub;
    Z.DTS = datenum(Z.DTS);
    Z = table2array(Z);
    
    for i = 1:length(ClassType)
        Class = ClassType{i};
        
        % Filter based on Class - MaxEvents is compromised after this (deleting entries)
        if strcmp(Class,'ClassOW')
            MaxEventsSub(MaxEventsSub.m == 1,:) = []; %#ok<*SAGROW>
            Z(Z(:,4) == 1,:) = [];
        elseif strcmp(Class,'Class')
            MaxEventsSub(MaxEventsSub.m == 1,:) = [];
            Z(Z(:,4) == 1,:) = [];
            MaxEventsSub(MaxEventsSub.m == 2,:) = [];
            Z(Z(:,4) == 2,:) = [];
        end
        
        for j = 1:length(BM)
            BlockM = BM{j};
            
            % Initialize
            Max(r).(Class).(BlockM) = [];
            
            if strcmp(BlockM,'Daily')
                
                % Make groups out of unique locations and days
                [Gr, ~, ~, ~] = findgroups(dateshift(MaxEventsSub.DTS,'start','day'),MaxEventsSub.SITE,MaxEventsSub.InfCase);
               
            elseif strcmp(BlockM,'Weekly')
                [Gr, ~, ~, ~] = findgroups(dateshift(MaxEventsSub.DTS,'start','week'),MaxEventsSub.SITE,MaxEventsSub.InfCase);
            else
                [Gr, ~, ~, ~] = findgroups(year(MaxEventsSub.DTS),MaxEventsSub.SITE,MaxEventsSub.InfCase);
            end
            
            % Perform splitapply (see function at end... not just Max as we want whole rows involving maxes)
            Max(r).(Class).(BlockM) = splitapply(@(Z)maxIndex(Z),Z,Gr);
            % Transform back into table form
            Max(r).(Class).(BlockM) = array2table(Max(r).(Class).(BlockM));
            %Max(r).(Class).(BlockM).Properties.VariableNames = {'R', 'SITE', 'Max', 'InfCase', 'DayRank', 'L1Veh', 'L2Veh', 'L1Load', 'L2Load', 'L1Ax', 'L2Ax', 'L1Sp', 'L2Sp', 'DTS', 'm'};
            Max(r).(Class).(BlockM).Properties.VariableNames = {'R', 'SITE', 'Max', 'InfCase', 'm', 'DayRank', 'BrStInd', 'DTS'};
            Max(r).(Class).(BlockM).DTS = datetime(Max(r).(Class).(BlockM).DTS,'ConvertFrom',"datenum"); Max(r).(Class).(BlockM).R = [];
            
        end
    end
end


% --- Step 4: Curve Fitting ---

% Set CDF Scaling Factors for estimates
D2WFactor = 5; W2YFactor = 50; D2YFactor = D2WFactor*W2YFactor;

% For each Influence Case
for r = 1:length(ILData)

    % Fit Block Maxima to Normal Curve
    for i = 1:length(ClassType)
        Class = ClassType{i};
        
        for j = 1:length(BM)
            BlockM = BM{j};
            
            for k = 1:length(DistTypes)
                Dist = DistTypes{k};
                
                [pd(r).(Class).(BlockM).(Dist),x_values(r).(Class).(BlockM),y_values(r).(Class).(BlockM).(Dist).PDF_Fit,...
                        y_values(r).(Class).(BlockM).(Dist).CDF_Fit] = GetBlockMaxFit(Max(r).(Class).(BlockM).Max,Dist,false);
         
                if strcmp(BlockM,'Daily')
                    y_values(r).(Class).('Weekly').(Dist).CDF_D2W = y_values(r).(Class).(BlockM).(Dist).CDF_Fit.^D2WFactor;
                    y_values(r).(Class).('Weekly').(Dist).PDF_D2W = [0 diff(y_values(r).(Class).('Weekly').(Dist).CDF_D2W)];
                    y_values(r).(Class).('Yearly').(Dist).CDF_D2Y = y_values(r).(Class).(BlockM).(Dist).CDF_Fit.^D2YFactor;
                    y_values(r).(Class).('Yearly').(Dist).PDF_D2Y = [0 diff(y_values(r).(Class).('Yearly').(Dist).CDF_D2Y)];
                    
                    if strcmp(Dist,'Normal')
                        mu = trapz(y_values(r).(Class).('Weekly').(Dist).CDF_D2W,x_values(r).(Class).(BlockM));
                        sigma = sqrt(trapz(y_values(r).(Class).('Weekly').(Dist).CDF_D2W,(x_values(r).(Class).(BlockM)-mu).^2));
                        pd(r).(Class).D2W.(Dist) = makedist(Dist,"mu",mu,"sigma",sigma);
                        
                        mu = trapz(y_values(r).(Class).('Yearly').(Dist).CDF_D2Y,x_values(r).(Class).(BlockM));
                        sigma = sqrt(trapz(y_values(r).(Class).('Yearly').(Dist).CDF_D2Y,(x_values(r).(Class).(BlockM)-mu).^2));
                        pd(r).(Class).D2Y.(Dist) = makedist(Dist,"mu",mu,"sigma",sigma);
                    elseif strcmp(Dist,'Lognormal')
                        mu = trapz(y_values(r).(Class).('Weekly').(Dist).CDF_D2W,x_values(r).(Class).(BlockM));
                        sigma = sqrt(trapz(y_values(r).(Class).('Weekly').(Dist).CDF_D2W,(x_values(r).(Class).(BlockM)-mu).^2));
                        zeta = sqrt(log(1+(sigma/mu)^2));
                        lambda = log(mu)-0.5*zeta^2;
                        pd(r).(Class).D2W.(Dist) = makedist(Dist,"mu",lambda,"sigma",zeta);
                        
                        mu = trapz(y_values(r).(Class).('Yearly').(Dist).CDF_D2Y,x_values(r).(Class).(BlockM));
                        sigma = sqrt(trapz(y_values(r).(Class).('Yearly').(Dist).CDF_D2Y,(x_values(r).(Class).(BlockM)-mu).^2));
                        zeta = sqrt(log(1+(sigma/mu)^2));
                        lambda = log(mu)-0.5*zeta^2;
                        pd(r).(Class).D2Y.(Dist) = makedist(Dist,"mu",lambda,"sigma",zeta);
                    end
                    
                elseif strcmp(BlockM,'Weekly')
                    y_values(r).(Class).('Yearly').(Dist).CDF_W2Y = y_values(r).(Class).(BlockM).(Dist).CDF_Fit.^W2YFactor;
                    y_values(r).(Class).('Yearly').(Dist).PDF_W2Y = [0 diff(y_values(r).(Class).('Yearly').(Dist).CDF_W2Y)];
                    
                    if strcmp(Dist,'Normal')
                        mu = trapz(y_values(r).(Class).('Yearly').(Dist).CDF_W2Y,x_values(r).(Class).(BlockM));
                        sigma = sqrt(trapz(y_values(r).(Class).('Yearly').(Dist).CDF_W2Y,(x_values(r).(Class).(BlockM)-mu).^2));
                        pd(r).(Class).W2Y.(Dist) = makedist(Dist,"mu",mu,"sigma",sigma);
                        
                    elseif strcmp(Dist,'Lognormal')
                        mu = trapz(y_values(r).(Class).('Yearly').(Dist).CDF_W2Y,x_values(r).(Class).(BlockM));
                        sigma = sqrt(trapz(y_values(r).(Class).('Yearly').(Dist).CDF_W2Y,(x_values(r).(Class).(BlockM)-mu).^2));
                        zeta = sqrt(log(1+(sigma/mu)^2));
                        lambda = log(mu)-0.5*zeta^2;
                        pd(r).(Class).W2Y.(Dist) = makedist(Dist,"mu",lambda,"sigma",zeta);
                        
                    end
                end
            end
        end
    end
end


end

function out = maxIndex(Z)
    [ymax, loc]=max(Z(:,2));
    out=[ymax, Z(loc,:)];
end
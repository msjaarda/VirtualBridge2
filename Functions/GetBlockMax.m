function [Max] = GetBlockMax(MaxEvents,ClassType,BM)
%GETBLOCKMAX Gives back Max, BlockMaximumData
% MaxEvents  - table including the maximum daily events (either axles  or load effects) w/ cols
% DTS  - datetime
% SITE - location
% Max OR MaxLE  - column we are 'block' finding for
% InfCase - influence line case (optional)
% m - 1 (Contains unclassified) 2 (Not 1, contains classified +) 3 (Not 1, 2, contains classified)
% Other cols (to be retained)
% Typical values of other vars:
% ClassType = {'All', 'ClassOW', 'Class'};            % i
% BM = {'Daily', 'Weekly', 'Monthly', 'Yearly'};      % j

% See if there is a column called ILData
try 
    AA = unique(MaxEvents.InfCase); NoInf = 0;
catch
    % If not, we do them all!
    AA = 1; MaxEvents.InfCase = ones(height(MaxEvents),1); NoInf = 1;
end

if ~isdatetime(MaxEvents.DTS(1))
    MaxEvents.DTS = datetime(MaxEvents.DTS,'ConvertFrom','datenum');
end

% Get Variable Names (for reconstructing table)
VN = string(MaxEvents.Properties.VariableNames);
% If we have a 'MaxLE', switch it to 'Max'
VN(strcmp("MaxLE",VN)) = "Max";

% Detect the column number of: DTS, SITE, and Max OR MaxLE, and m
DTSID = find(VN == "DTS");
mID = find(VN == "m");
SITEID = find(VN == "SITE");
InfCaseID = find(VN == "InfCase");
MaxID = find(VN == "Max");
if isempty(MaxID)
    MaxID = find(VN == "MaxLE");
end

% --- Build Structure with Block Maxima ---

% For each Influence case
for r = 1:length(AA)
    
    % Reset MaxEvents and select r
    %MaxEventsSub = MaxEvents(MaxEvents.InfCase == r,:);
    % Transform AxTandem into Array (necessary for splitapply)
    %Z = MaxEventsSub;
    
    Z = MaxEvents(MaxEvents.InfCase == r,:);
    Z.DTS = datenum(Z.DTS);
    Z = table2array(Z);
    
    for i = 1:length(ClassType)
        Class = ClassType{i};
        
        % Filter based on Class - MaxEvents is compromised after this (deleting entries)
        if strcmp(Class,'ClassOW')
            Z(Z(:,mID) == 1,:) = [];
        elseif strcmp(Class,'Class')
            Z(Z(:,mID) == 1 | Z(:,mID) == 2,:) = [];
        end
        
        for j = 1:length(BM)
            BlockM = BM{j};
            
            % Initialize
            Max(r).(Class).(BlockM) = [];
            
            if strcmp(BlockM,'Daily')
                % Make groups out of unique locations and days
                [Gr, ~, ~, ~] = findgroups(dateshift(datetime(Z(:,DTSID),'ConvertFrom','datenum'),'start','day'),Z(:,SITEID),Z(:,InfCaseID));
            elseif strcmp(BlockM,'Weekly')
                [Gr, ~, ~, ~] = findgroups(dateshift(datetime(Z(:,DTSID),'ConvertFrom','datenum'),'start','week'),Z(:,SITEID),Z(:,InfCaseID));
            elseif strcmp(BlockM,'Monthly')
                [Gr, ~, ~, ~] = findgroups(dateshift(datetime(Z(:,DTSID),'ConvertFrom','datenum'),'start','month'),Z(:,SITEID),Z(:,InfCaseID));
            else % Yearly
                [Gr, ~, ~, ~] = findgroups(year(datetime(Z(:,DTSID),'ConvertFrom','datenum')),Z(:,SITEID),Z(:,InfCaseID));
            end
            
            % Perform splitapply (see function at end... not just Max as we want whole rows involving maxes)
            Max(r).(Class).(BlockM) = splitapply(@(Z)maxIndex(Z,BlockM,DTSID,MaxID),Z,Gr);
            % Delete -1 values
            Max(r).(Class).(BlockM)(Max(r).(Class).(BlockM)(:,1) == -1,:) = []; Max(r).(Class).(BlockM)(:,1) = [];
            % Transform back into table form
            Max(r).(Class).(BlockM) = array2table(Max(r).(Class).(BlockM));
            Max(r).(Class).(BlockM).Properties.VariableNames = VN;
            Max(r).(Class).(BlockM).DTS = datetime(Max(r).(Class).(BlockM).DTS,'ConvertFrom',"datenum");
            if NoInf % Return the structure without the InfCase var
                Max(r).(Class).(BlockM).InfCase = [];
            end
        end
    end
end
end

function out = maxIndex(Z,BlockM,DTSID,MaxID)
    if strcmp(BlockM,'Yearly')
        LimH = 0.6*365*24;
    elseif strcmp(BlockM,'Weekly')
        LimH = 4*24;
    elseif strcmp(BlockM,'Monthly')
        LimH = 18*24;
    else
        LimH = 12;
    end
    
    if hours(max(datetime(Z(:,DTSID),'ConvertFrom','datenum')) - min(datetime(Z(:,DTSID),'ConvertFrom','datenum'))) < LimH
        % We cannot have the daily limit of 12 if we only have 1 daily value
        if size(Z,2) < 100
            [ymax, loc] = max(Z(:,MaxID));
            out = [ymax, Z(loc,:)];
        else
            out = [-1, Z(1,:)];
        end
    else
        [ymax, loc] = max(Z(:,MaxID));
        out = [ymax, Z(loc,:)];
    end
    
end

% Raph Explanation

% Anonymous functions will look for variables in the workspace... that is
% why BlockM still goes.

% Two functions are being called by splitapply kinda.

% function out = function_handle(Z)
%     out = maxIndexTridem(Z)
% end
% 
% function out = @(Z)
%     BlockM = "Yearly";
%     out = maxIndexTridem(Z,BlockM);
% end
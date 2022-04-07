% ------------------------------------------------------------------------
%                             VBAxleStatsBasic
% ------------------------------------------------------------------------
% Assemble Simple, Tandem, and Tridem Axles, using geometry, to gain
% information on maximum axle loads (especially Q1)
% Successor to AxleStatsBasic... which was much slower using WIMtoAllTr...
% Differs from function AxleStats (func) because it finds axle groups in
% non-classified vehicles. Optional variable save at end (large vars)
% Recommended instead to only save the blockmax result, in structure

 % Initial commands
clear, clc, tic, close all

% Inputs
Grav = 9.81;
TAxSingle = false;
TAxTandem = false;
TAxTridem = true;

Country = 'SWISS';
%Country = 'GER';
%Country = 'AUST';
%Country = 'USA';
% Can put another name here instead if a country name, then set Sitesx

% Former
%Sites = [402 405:409 415:416];
% Replaced by ALL
%Sites = [402 405:423 431:434 441 456 489];
% Trial
%Sitesx = [177];
% ALL from that Country
Sitesx = 0;

load('C:\Users\mjsja\Desktop\SwissTraffic2\Misc\Sites.mat')
% Consider, for Swiss, incorporating Calibration info

% Initialize
VNSingle = {'AWT1kN','SITE','CLASS','DTS'};
AxSingle = cell2table(cell(0,4),'VariableNames',VNSingle);
VNTandem = {'AWT1kN','AWT2kN','W1_2','SITE','CLASS','DTS'};
AxTandem = cell2table(cell(0,6),'VariableNames',VNTandem);
VNTridem = {'AWT1kN','AWT2kN','AWT3kN','W1_2','W2_3','SITE','CLASS','DTS'};
AxTridem = cell2table(cell(0,8),'VariableNames',VNTridem);

% Get List of Sites
if Sitesx == 0
    % Take all Core for that Country
    Sitesx = Sites.SITE(Sites.(Country) == 1 & Sites.Core == 1);
else
    Country = '';
end

% Get Axle Groups
for i = 1:length(Sitesx)
    % Get individual SITE
    SITE = Sitesx(i);
    % Load WIM File
    load(['C:\Users\mjsja\Desktop\SwissTraffic2\WIM\' num2str(SITE) '.mat'])
    %PDs.DTS = datenum(PDs.DTS);
    
    % Perform Stage2Prune
    PDs = Stage2Prune(PDs);
    % DON'T NEED TO REMOVE DUPS (Save Time)... Makes no difference
    % Find and delete duplicates if possible... check if required
%     if sum(SITE == Sites.SITE(Sites.Uni2L == 1)) > 0
%         % Get Duplicates
%         PDs = FindDup2(PDs,0,0);
%         % Delete Duplicates - from L1
%         PDs(PDs.Dup & PDs.LANE == 1,:) = [];
%     end
            
    % Find all tandem axles
    IndW = find(string(PDs.Properties.VariableNames) == "W1_2");
    NumW = sum(cell2mat(regexp(string(PDs.Properties.VariableNames), 'W\d_*')));
    WB = PDs{:,IndW:IndW+NumW-1}; WB(PDs{:,IndW:IndW+NumW-1} == 0) = NaN;
    IndA = find(string(PDs.Properties.VariableNames) == "AWT01");
    NumA = sum(cell2mat(regexp(string(PDs.Properties.VariableNames), 'AWT')));
    AXkN = PDs{:,IndA:IndA+NumA-1}; AXkN(PDs{:,IndA:IndA+NumA-1} == 0) = NaN;
    AXkN = AXkN.*(Grav/1000);
    
    % Create the table of axeles. 0: >2m ; 1: <2m, NaN: no axle
    axleTable = double(PDs{:,IndW:IndW+NumW-1} < 200);
    axleTable(PDs{:,IndW:IndW+NumW-1} == 0) = NaN;
    axleTable(:,size(axleTable,2)+1) = NaN;
    
    % Initilize the empty count of axles and count of axles per group.
    count = ones(size(axleTable,1),1);
    countGroup = nan(size(axleTable));
    Finished = zeros(size(axleTable,1),1);

    for j = 1:size(axleTable,2)
        
        % Find vehicles that are finished (first NaN)
        idFinish = isnan(axleTable(:,j));
        % Find vehicules with current axle seperated from the previous group
        idAlone = axleTable(:,j) == 0;
        % Find vehicules with current axle being part of the previous group
        idTogether = axleTable(:,j) == 1;
        % If it is part of the previous group, then add 1 to the axle count
        count(idTogether) = count(idTogether) + 1;
        % If new group, first write down the count of the previous group
        countGroup(idAlone,j) = count(idAlone);
        % Also add those that are finished
        countGroup(idFinish & ~Finished,j) = count(idFinish & ~Finished);
        % Finish these
        Finished(idFinish) = 1;
        % And reinitilize the count of axle in the new group to 1
        count(idAlone) = 1;
        
    end
    
    if TAxSingle
        % First do Single
        Y = countGroup == 1;
        Ind = find(Y);
        [R, ~] = ind2sub([size(countGroup,1) size(countGroup,2)],Ind);
        AxSingle = [AxSingle; array2table([AXkN(Y) PDs.SITE(R) PDs.CLASS(R) datenum(PDs.DTS(R))],'VariableNames',VNSingle)];
    else
        AxSingle = [];
    end
    
    if TAxTandem
        % Next Tandem
        % Now get axles themselves using countGroup and AXkN
        % Grab all 1s (Single), 2s (Tandem), and 3s (Tridem) in countGroup.
        % For Tandem, use these indices in AXkN. These are all AWT2kN.
        % Then go 1 index to the left and these are AWT1kN. Same index is W1_2.
        Y = countGroup == 2;
        X = circshift(Y,-1,2);
        
        Ind = find(Y);
        [R, ~] = ind2sub([size(countGroup,1) size(countGroup,2)],Ind);
        AxTandem = [AxTandem; array2table([AXkN(X) AXkN(Y) WB(X)./100 PDs.SITE(R) PDs.CLASS(R) datenum(PDs.DTS(R))],'VariableNames',VNTandem)];
    else
        AxTandem = [];
    end
    
    if TAxTridem
        % Finally Tridem
        Y = countGroup == 3;
        X = circshift(Y,-1,2);
        X2 = circshift(Y,-2,2);
        
        Ind = find(Y);
        [R, ~] = ind2sub([size(countGroup,1) size(countGroup,2)],Ind);
        AxTridem = [AxTridem; array2table([AXkN(X2) AXkN(X) AXkN(Y) WB(X2)./100 WB(X)./100 PDs.SITE(R) PDs.CLASS(R) datenum(PDs.DTS(R))],'VariableNames',VNTridem)];
    else
        AxTridem = [];
    end
    
    fprintf('\nSITE: %i (%i of %i)',SITE, i, length(Sitesx))
    
end

% Let's do block max right here before saving the much smaller BM version...

fprintf('\nTotal time: %.2f seconds\n\n',toc)
clearvars -except AxSingle AxTandem AxTridem TAxSingle TAxTandem TAxTridem Country

% NEXT STEP: BLOCK MAXIMA

% Block Maxima (always use j)
BM = {'Daily', 'Weekly', 'Yearly'};
% Select which Classification you want to analyze (always use i)
ClassType = {'All', 'ClassOW', 'Class'};
% Must be in the order 'All' 'Class+' 'Class' due to deletions
ClassT = {'All', 'Classified+', 'Classified'};

% Variable Names
VNSingle = {'Max','AWT1kN','SITE','CLASS','DTS'};
VNTandem = {'Max','AWT1kN','AWT2kN','W1_2M','SITE','CLASS','DTS'};
VNTridem = {'Max','AWT1kN','AWT2kN','AWT3kN','W1_2M','W2_3M','SITE','CLASS','DTS'};

if TAxSingle
    Z = table2array(AxSingle);
    AxSingle.DTS = datetime(AxSingle.DTS,'ConvertFrom','datenum');
    % Delete Cols for Speed
    AxSingle.AWT1kN = [];
    
    for i = 1:length(ClassType)
        Class = ClassType{i};
        
        % Filter based on Class - AxSingle is compromised after this (deleting entries)
        if strcmp(Class,'ClassOW')
            AxSingle(AxSingle.CLASS == 0,:) = [];
            Z(Z(:,3) == 0,:) = [];
        elseif strcmp(Class,'Class')
            AxSingle(AxSingle.CLASS == 0,:) = [];
            Z(Z(:,3) == 0,:) = [];
            AxSingle(AxSingle.CLASS > 39 & AxSingle.CLASS < 50,:) = [];
            Z(Z(:,3) > 39 & Z(:,3) < 50,:) = [];
        end
        
        for j = 1:length(BM)
            BlockM = BM{j};
            
            % Initialize
            MaxAx.Single.(Class).(BlockM) = [];
            
            if strcmp(BlockM,'Daily')
                % Make groups out of unique locations and days
                [Gr, GrIDDay, GrIDZST] = findgroups(dateshift(AxSingle.DTS,'start','day'),AxSingle.SITE);
            elseif strcmp(BlockM,'Weekly')
                [Gr, GrIDWeek, GrIDZST] = findgroups(dateshift(AxSingle.DTS,'start','week'),AxSingle.SITE);
            else
                [Gr, GrIDYear, GrIDZST] = findgroups(year(AxSingle.DTS),AxSingle.SITE);
            end
            
            % Perform splitapply (see function at end... not just Max as we want whole rows involving maxes)
            MaxAx.Single.(Class).(BlockM) = splitapply(@(Z)maxIndexSingle(Z),Z,Gr);
            % Transform back into table form
            MaxAx.Single.(Class).(BlockM) = array2table(MaxAx.Single.(Class).(BlockM));
            MaxAx.Single.(Class).(BlockM).Properties.VariableNames = VNSingle;
            MaxAx.Single.(Class).(BlockM).DTS = datetime(MaxAx.Single.(Class).(BlockM).DTS,'ConvertFrom',"datenum");
            
        end
    end
    % Add to overall struct
end

if TAxTandem
    Z = table2array(AxTandem);
    AxTandem.DTS = datetime(AxTandem.DTS,'ConvertFrom','datenum');
    % Delete Cols for Speed
    AxTandem.AWT1kN = []; AxTandem.AWT2kN = []; AxTandem.W1_2 = [];
    
    for i = 1:length(ClassType)
        Class = ClassType{i};
        
        % Filter based on Class - AxTandem is compromised after this (deleting entries)
        if strcmp(Class,'ClassOW')
            AxTandem(AxTandem.CLASS == 0,:) = [];
            Z(Z(:,5) == 0,:) = [];
        elseif strcmp(Class,'Class')
            AxTandem(AxTandem.CLASS == 0,:) = [];
            Z(Z(:,5) == 0,:) = [];
            AxTandem(AxTandem.CLASS > 39 & AxTandem.CLASS < 50,:) = [];
            Z(Z(:,5) > 39 & Z(:,5) < 50,:) = [];
        end
        
        for j = 1:length(BM)
            BlockM = BM{j};
            
            % Initialize
            MaxAx.Tandem.(Class).(BlockM) = [];
            
            if strcmp(BlockM,'Daily')
                % Make groups out of unique locations and days
                [Gr, GrIDDay, GrIDZST] = findgroups(dateshift(AxTandem.DTS,'start','day'),AxTandem.SITE);
            elseif strcmp(BlockM,'Weekly')
                [Gr, GrIDWeek, GrIDZST] = findgroups(dateshift(AxTandem.DTS,'start','week'),AxTandem.SITE);
            else
                [Gr, GrIDYear, GrIDZST] = findgroups(year(AxTandem.DTS),AxTandem.SITE);
            end
            
            % Perform splitapply (see function at end... not just Max as we want whole rows involving maxes)
            MaxAx.Tandem.(Class).(BlockM) = splitapply(@(Z)maxIndexTandem(Z),Z,Gr);
            % Transform back into table form
            MaxAx.Tandem.(Class).(BlockM) = array2table(MaxAx.Tandem.(Class).(BlockM));
            MaxAx.Tandem.(Class).(BlockM).Properties.VariableNames = VNTandem;
            MaxAx.Tandem.(Class).(BlockM).W1_2M = MaxAx.Tandem.(Class).(BlockM).W1_2M;
            MaxAx.Tandem.(Class).(BlockM).DTS = datetime(MaxAx.Tandem.(Class).(BlockM).DTS,'ConvertFrom',"datenum");
            
        end
    end
    % Add to overall struct
end

if TAxTridem
    Z = table2array(AxTridem);
    AxTridem.DTS = datetime(AxTridem.DTS,'ConvertFrom','datenum');
    % Delete Cols for Speed
    AxTridem.AWT1kN = []; AxTridem.AWT2kN = []; AxTridem.AWT3kN = []; AxTridem.W1_2 = []; AxTridem.W2_3 = [];
    
    for i = 1:length(ClassType)
        Class = ClassType{i};
        
        % Filter based on Class - AxTridem is compromised after this (deleting entries)
        if strcmp(Class,'ClassOW')
            AxTridem(AxTridem.CLASS == 0,:) = [];
            Z(Z(:,7) == 0,:) = [];
        elseif strcmp(Class,'Class')
            AxTridem(AxTridem.CLASS == 0,:) = [];
            Z(Z(:,7) == 0,:) = [];
            AxTridem(AxTridem.CLASS > 39 & AxTridem.CLASS < 50,:) = [];
            Z(Z(:,7) > 39 & Z(:,7) < 50,:) = [];
        end
        
        for j = 1:length(BM)
            BlockM = BM{j};
            
            % Initialize
            MaxAx.Tridem.(Class).(BlockM) = [];
            
            if strcmp(BlockM,'Daily')
                % Make groups out of unique locations and days
                [Gr, GrIDDay, GrIDZST] = findgroups(dateshift(AxTridem.DTS,'start','day'),AxTridem.SITE);
            elseif strcmp(BlockM,'Weekly')
                [Gr, GrIDWeek, GrIDZST] = findgroups(dateshift(AxTridem.DTS,'start','week'),AxTridem.SITE);
            else
                [Gr, GrIDYear, GrIDZST] = findgroups(year(AxTridem.DTS),AxTridem.SITE);
            end
            
            % Perform splitapply (see function at end... not just Max as we want whole rows involving maxes)
            MaxAx.Tridem.(Class).(BlockM) = splitapply(@(Z)maxIndexTridem(Z),Z,Gr);
            % Transform back into table form
            MaxAx.Tridem.(Class).(BlockM) = array2table(MaxAx.Tridem.(Class).(BlockM));
            MaxAx.Tridem.(Class).(BlockM).Properties.VariableNames = VNTridem;
            MaxAx.Tridem.(Class).(BlockM).W1_2M = MaxAx.Tridem.(Class).(BlockM).W1_2M;
            MaxAx.Tridem.(Class).(BlockM).W2_3M = MaxAx.Tridem.(Class).(BlockM).W2_3M;
            MaxAx.Tridem.(Class).(BlockM).DTS = datetime(MaxAx.Tridem.(Class).(BlockM).DTS,'ConvertFrom',"datenum");
            
        end
    end
    % Add to overall struct
end

% MIGHT NEED TO ADD A FILTER TO REMVOE OUTLIERS ON THE LOW SIDE...
% HAPPENS WHEN YEARLY BLOCKMAX COMES FROM STATION WITH ONLY 2 DAYS FOR EX.
% NOW I HAVE GONE AND DELETED ONE BY HAND (80 for tridem SWISS yearly)


% Optional Saving
% Load MaxAx
load('Misc/BMAxles.mat')
%BMAxles.(Country).MaxAx = MaxAx;
%BMAxles.(Country).MaxAx.Tridem = MaxAx.Tridem;

%save('Misc/BMAxles.mat','BMAxles')

fprintf('\nTotal time: %.2f seconds\n\n',toc)

function out = maxIndexSingle(Z)
    [ymax, loc]=max(sum([Z(:,1)],2));
    out=[ymax, Z(loc,:)];
end

function out = maxIndexTandem(Z)
    [ymax, loc]=max(sum([Z(:,1) Z(:,2)],2));
    out=[ymax, Z(loc,:)];
end

function out = maxIndexTridem(Z)
    [ymax, loc]=max(sum([Z(:,1) Z(:,2) Z(:,3)],2));
    out=[ymax, Z(loc,:)];
end
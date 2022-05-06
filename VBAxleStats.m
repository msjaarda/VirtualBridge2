% ------------------------------------------------------------------------
%                             VBAxleStats
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
TAxSingle = true;
TAxTandem = true;
TAxTridem = true;

Country = 'CH';
%Country = 'DE';
%Country = 'AT';
%Country = 'US';
%Country = '';
% Can put another name here instead if a country name, then set Sites

% Former
%Sites = [402 405:409 415:416];
% Replaced by ALL
%Sites = [402 405:423 431:434 441 456 489];
% Trial
%Sitesx = [177];
% ALL from that Country
Sites = 0;

% Get SiteSet
Sites = VBGetSiteSet(Sites,1,1,Country);

% Consider, for Swiss, incorporating Calibration info

Grav = 9.81;

% Initialize
VNSingle = {'AWT1kN','SITE','CLASS','DTS'};
AxSingle = cell2table(cell(0,4),'VariableNames',VNSingle);
VNTandem = {'AWT1kN','AWT2kN','W1_2','SITE','CLASS','DTS'};
AxTandem = cell2table(cell(0,6),'VariableNames',VNTandem);
VNTridem = {'AWT1kN','AWT2kN','AWT3kN','W1_2','W2_3','SITE','CLASS','DTS'};
AxTridem = cell2table(cell(0,8),'VariableNames',VNTridem);
AxSingle = []; AxTandem = []; AxTridem = [];

% Get Axle Groups
for i = 1:length(Sites)
    % Get individual SITE
    SITE = Sites(i);
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
    end
    
    if TAxTridem
        % Finally Tridem
        Y = countGroup == 3;
        X = circshift(Y,-1,2);
        X2 = circshift(Y,-2,2);
        
        Ind = find(Y);
        [R, ~] = ind2sub([size(countGroup,1) size(countGroup,2)],Ind);
        AxTridem = [AxTridem; array2table([AXkN(X2) AXkN(X) AXkN(Y) WB(X2)./100 WB(X)./100 PDs.SITE(R) PDs.CLASS(R) datenum(PDs.DTS(R))],'VariableNames',VNTridem)];
    end
    
    fprintf('\nSITE: %i (%i of %i)',SITE, i, length(Sites))
    
end

% Let's do block max right here before saving the much smaller BM version...
fprintf('\nTotal time: %.2f seconds\n\n',toc)
clearvars -except AxSingle AxTandem AxTridem TAxSingle TAxTandem TAxTridem Country

% NEXT STEP: BLOCK MAXIMA

% Block Maxima (always use j)
BM = {'Daily', 'Weekly', 'Monthly', 'Yearly'};
% Select which Classification you want to analyze (always use i)
ClassType = {'All', 'ClassOW', 'Class'};


% % Get Variable Names (for reconstructing table)
% VN = string(MaxEvents.Properties.VariableNames);
% 
% % Detect the column number of: DTS, SITE, and Max OR MaxLE, and m
% DTSID = find(VN == "DTS");
% mID = find(VN == "m");
% SITEID = find(VN == "SITE");
% InfCaseID = find(VN == "InfCase");
% MaxID = find(VN == "Max");
% if isempty(MaxID)
%     MaxID = find(VN == "MaxLE");
% end
% 
% if ~isdatetime(MaxEvents.DTS(1))
%     MaxEvents.DTS = datetime(MaxEvents.DTS,'ConvertFrom','datenum');
% end
% 



if TAxSingle
    % Let's trim in advance based on experience...
    AxSingle.Max = AxSingle.AWT1kN;
    AxSingle(AxSingle.Max < 50,:) = [];
    AxSingle.m = ones(height(AxSingle),1); AxSingle.m(AxSingle.CLASS > 0) = 3;
    AxSingle.m(AxSingle.CLASS > 39 & AxSingle.CLASS < 50) = 2;
    MaxAx.Single = GetBlockMax(AxSingle,ClassType,BM);
end
if TAxTandem
    AxTandem.Max = AxTandem.AWT1kN + AxTandem.AWT2kN;
    AxTandem(AxTandem.Max < 90,:) = [];
    AxTandem.m = ones(height(AxTandem),1); AxTandem.m(AxTandem.CLASS > 0) = 3;
    AxTandem.m(AxTandem.CLASS > 39 & AxTandem.CLASS < 50) = 2;
    MaxAx.Tandem = GetBlockMax(AxTandem,ClassType,BM);
end
if TAxTridem
    AxTridem.Max = AxTridem.AWT1kN + AxTridem.AWT2kN + AxTridem.AWT3kN;
    AxTridem(AxTridem.Max < 120,:) = [];
    AxTridem.m = ones(height(AxTridem),1); AxTridem.m(AxTridem.CLASS > 0) = 3;
    AxTridem.m(AxTridem.CLASS > 39 & AxTridem.CLASS < 50) = 2;
    MaxAx.Tridem = GetBlockMax(AxTridem,ClassType,BM);
end

% Load MaxAx
load('Misc/BMAxles.mat')
BMAxles.(Country).MaxAx = MaxAx;
% Optional Saving
save('Misc/BMAxles.mat','BMAxles')

fprintf('\nTotal time: %.2f seconds\n\n',toc)


% 
% function out = maxIndex(Z,BlockM,DTSID,MaxID)
%     if strcmp(BlockM,'Yearly')
%         LimH = 0.6*365*24;
%     elseif strcmp(BlockM,'Weekly')
%         LimH = 4*24;
%     elseif strcmp(BlockM,'Monthly')
%         LimH = 18*24;
%     else
%         LimH = 12;
%     end
%     
%     if hours(max(datetime(Z(:,DTSID),'ConvertFrom','datenum')) - min(datetime(Z(:,DTSID),'ConvertFrom','datenum'))) < LimH
%         out = [-1, Z(1,:)];
%     else
%         [ymax, loc] = max(Z(:,MaxID));
%         out = [ymax, Z(loc,:)];
%     end
%     
% end
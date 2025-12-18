% PRUNE CLASSIFY SAVE COMPLETE
% Combine all years for each stations into one variable
% Make fixes along the way

clear, clc, tic

% INPUT
Folder = 'Raw WIM\';
Files = dir(Folder);
Files(1:2) = [];
SNames = {Files(:).name};

load('Sites.mat')
load('SiteLanes.mat')

% Initialize
% Sites.StartDate = repmat(datetime('now'),height(Sites),1);
% Sites.EndDate = repmat(datetime('now'),height(Sites),1);
% Sites.TrRateEst = zeros(height(Sites),1);

for i = 1:length(SNames)
    
    SName = SNames{i};
    
    % Get File List
    Folder = ['Raw WIM\' SName];
    Files = dir(Folder);
    Files(1:2) = [];
    
    % Initialize
    PD = [];
    
    % Go through and add up RD
    for p = 1:length(Files)
        % Method to not override (poof) RD
        RDi = load(strcat(Files(p).folder,'\',Files(p).name));
        PD = [PD; RDi.RD];
    end
    
    % Get all nans
    A = sum(isnan(PD{:,3:end}));
    
    % Loop through and look for other nans
    for p = 3:29
        if A(p-2) > 0
            PD(isnan(PD{:,p}),:) = [];
        end
    end
    PD(isnat(PD.DTS),:) = [];
    PD(PD.LANE == 0,:) = [];
    PD.DIR = logical(PD.DIR);
    
    % ---------- START PRUNE -------------------
    
    % Use Method from ASTRA Reports
    
    % 1) Remove vehicles under 3.5 tonnes
    
    UnderW = PD.GW_TOT<3500;
    TotUnderW = sum(UnderW);
    PD(UnderW,:) = [];
    
    % Step 2 in ASTRA Reports is not a real step
    
    % 3) Remove vehicles with no lengths
    
    NoL = PD.LENTH<1;
    TotNoL = sum(NoL);
    PD(NoL,:) = [];
    
    % 4) Remove vehicles with lengths over 26m
    
    OverL = PD.LENTH>2600;
    TotOverL = sum(OverL);
    PD(OverL,:) = [];
    
    % 5) Remove vehicles with missing axle weights
    
    % Get column names starting with AWT
    InAxs = contains(PD.Properties.VariableNames, 'AWT');
    
    % Ugly but fast
    NoAW = PD.AX == 1 |...
        PD.AX == 2 & (PD.AWT01<1 | PD.AWT02<1 | PD.AWT03>1 | PD.AWT04>1 | PD.AWT05>1 | PD.AWT06>1 | PD.AWT07>1 | PD.AWT08>1 | PD.AWT09>1) |...
        PD.AX == 3 & (PD.AWT01<1 | PD.AWT02<1 | PD.AWT03<1 | PD.AWT04>1 | PD.AWT05>1 | PD.AWT06>1 | PD.AWT07>1 | PD.AWT08>1 | PD.AWT09>1) |...
        PD.AX == 4 & (PD.AWT01<1 | PD.AWT02<1 | PD.AWT03<1 | PD.AWT04<1 | PD.AWT05>1 | PD.AWT06>1 | PD.AWT07>1 | PD.AWT08>1 | PD.AWT09>1) |...
        PD.AX == 5 & (PD.AWT01<1 | PD.AWT02<1 | PD.AWT03<1 | PD.AWT04<1 | PD.AWT05<1 | PD.AWT06>1 | PD.AWT07>1 | PD.AWT08>1 | PD.AWT09>1) |...
        PD.AX == 6 & (PD.AWT01<1 | PD.AWT02<1 | PD.AWT03<1 | PD.AWT04<1 | PD.AWT05<1 | PD.AWT06<1 | PD.AWT07>1 | PD.AWT08>1 | PD.AWT09>1) |...
        PD.AX == 7 & (PD.AWT01<1 | PD.AWT02<1 | PD.AWT03<1 | PD.AWT04<1 | PD.AWT05<1 | PD.AWT06<1 | PD.AWT07<1 | PD.AWT08>1 | PD.AWT09>1) |...
        PD.AX == 8 & (PD.AWT01<1 | PD.AWT02<1 | PD.AWT03<1 | PD.AWT04<1 | PD.AWT05<1 | PD.AWT06<1 | PD.AWT07<1 | PD.AWT08<1 | PD.AWT09>1) |...
        PD.AX == 9 & (PD.AWT01<1 | PD.AWT02<1 | PD.AWT03<1 | PD.AWT04<1 | PD.AWT05<1 | PD.AWT06<1 | PD.AWT07<1 | PD.AWT08<1 | PD.AWT09<1);
    TotNoAW = sum(NoAW);
    PD(NoAW,:) = [];
    
    % 6) Remove vehicles with axle/wheelbase distances less than 60 cm
    
    % Get column names starting with AWT
    InWbs = logical(contains(PD.Properties.VariableNames, 'W').*contains(PD.Properties.VariableNames, '_'));
    InWbs(find(string(PD.Properties.VariableNames) == 'GW_TOT')) = 0;
    
    % Ugly but fast
    UnderWB = PD.AX == 1 |...
        PD.AX == 2 & (PD.W1_2<60 | PD.W2_3>1 | PD.W3_4>1 | PD.W4_5>1 | PD.W5_6>1 | PD.W6_7>1 | PD.W7_8>1 | PD.W8_9>1) |...
        PD.AX == 3 & (PD.W1_2<60 | PD.W2_3<60 | PD.W3_4>1 | PD.W4_5>1 | PD.W5_6>1 | PD.W6_7>1 | PD.W7_8>1 | PD.W8_9>1) |...
        PD.AX == 4 & (PD.W1_2<60 | PD.W2_3<60 | PD.W3_4<60 | PD.W4_5>1 | PD.W5_6>1 | PD.W6_7>1 | PD.W7_8>1 | PD.W8_9>1) |...
        PD.AX == 5 & (PD.W1_2<60 | PD.W2_3<60 | PD.W3_4<60 | PD.W4_5<60 | PD.W5_6>1 | PD.W6_7>1 | PD.W7_8>1 | PD.W8_9>1) |...
        PD.AX == 6 & (PD.W1_2<60 | PD.W2_3<60 | PD.W3_4<60 | PD.W4_5<60 | PD.W5_6<60 | PD.W6_7>1 | PD.W7_8>1 | PD.W8_9>1) |...
        PD.AX == 7 & (PD.W1_2<60 | PD.W2_3<60 | PD.W3_4<60 | PD.W4_5<60 | PD.W5_6<60 | PD.W6_7<60 | PD.W7_8>1 | PD.W8_9>1) |...
        PD.AX == 8 & (PD.W1_2<60 | PD.W2_3<60 | PD.W3_4<60 | PD.W4_5<60 | PD.W5_6<60 | PD.W6_7<60 | PD.W7_8<60 | PD.W8_9>1) |...
        PD.AX == 9 & (PD.W1_2<60 | PD.W2_3<60 | PD.W3_4<60 | PD.W4_5<60 | PD.W5_6<60 | PD.W6_7<60 | PD.W7_8<60 | PD.W8_9<60);
    TotUnderWB = sum(UnderWB);
    PD(UnderWB,:) = [];
    
    % 7) Eliminate Max Total Weight greater than 100 tonnes (ASTRA does 65)
    
    OverW = PD.GW_TOT>100000;
    TotOverW = sum(OverW);
    PD(OverW,:) = [];
    
    % 8) Remove vehicles with axle weights greater than 25 tonnes (ASTRA does 18)
    
    OverAW = sum(PD{:,InAxs} > 25000,2)>0;
    TotOverAW = sum(OverAW);
    PD(OverAW,:) = [];
    
    % 9) Remove vehicles with total lengths < 4m
    
    UnderL = PD.LENTH<400;
    TotUnderL = sum(UnderL);
    PD(UnderL,:) = [];
    
    % We perform some additional steps of our own ----------------------------
    
    % A # Axles should be 9 or less
    
    OverAx = PD.AX>9;
    TotOverAx = sum(OverAx);
    PD(OverAx,:) = [];
    
    % B Individual wheelbases should be less than 15 m
    
    OverWB = sum(PD{:,InWbs} > 1500,2)>0;
    TotOverWB = sum(OverWB);
    PD(OverWB,:) = [];
    
    % ---------- PRUNE DONE -------------------
    
    PD = VBClassify(PD);    

    % Fix Ceneri - 408 is used for both directions before 23-Nov-2006
    % Call it 489
    if strcmp(SName,'Ceneri')
        % Delete the error spots
        PD(PD.SITE == 408 & PD.LANE < 3 & PD.DTS > datetime(2014,1,1),:) = [];
        % Find date where new stations were installed
        DChange = max(PD.DTS(PD.SITE == 408 & PD.LANE < 3));
        % Set Old station to 489 instead of 408
        PD.SITE(PD.SITE == 408 & PD.DTS <= DChange) = 489;
        % Rename 408 Lane 3 as 1 and 4 as 1
        PD.LANE(PD.SITE == 408 & PD.LANE == 3) = 2;
        PD.LANE(PD.SITE == 408 & PD.LANE == 4) = 1;
    end
    % Fix Denges - 405 is used before 2010 for 1 lane in each direction
    % Call it 456
    if strcmp(SName,'Denges')
        % Find date where new stations were installed
        DChange = datetime(2010,6,6);
        % Rename 405 Lane 1 as 1 and 406 Lane 1 as 2
        PD.LANE(PD.DTS <= DChange & PD.SITE == 406) = 2;
        % Set Old station to 456 instead of 405 and 406
        PD.SITE(PD.DTS <= DChange) = 456;
    end
    
    
    % Fix SanBernardino 439/440 and Simplon, 441/442 to be at same SITE
    if strcmp(SName,'SanBernardino')
        % Delete sites 439 and 440 as they are very small
        PD(PD.SITE > 423,:) = [];
    end
    if strcmp(SName,'Simplon')
        % Change lane of 442 to 2
        PD.LANE(PD.SITE == 442) = 2;
        PD.SITE(PD.SITE == 442) = 441;
    end
    
    % Go through and correct directions, COUNTIDs, and SITE LANE Errors
    % Separate into Stations... or even lanes?
    SiteDetails = Sites(strcmp(Sites.SName,SName),:);
    
    for j = 1:height(SiteDetails)
        SITE = SiteDetails.SITE(j);

        LaneDetails = SiteLanes(SiteLanes.SITE == SITE,:);
        % True for DIR is North and East
        if length(unique(LaneDetails.DIR)) > 1
            for k = 1:height(LaneDetails)
                if LaneDetails.NSEW(k) == 1 || LaneDetails.NSEW(k) == 3
                    PD.DIR(PD.SITE == SITE & PD.LANE == LaneDetails.LANE(k)) = true;
                else
                    PD.DIR(PD.SITE == SITE & PD.LANE == LaneDetails.LANE(k)) = false;
                end
            end            
        end
        
        PDs = PD(PD.SITE == SITE,:);
        % Fix problem specific to 405 and 406 with PUN  NEVER BEEN DONE AFTER
        % IMPORT YET!
        if SITE == 405 || SITE == 406
            PDs(PDs.LANE == 1 & PDs.GAPT > 99.8,:) = [];
        end
        SiteDetails.StartDate(j) = min(PDs.DTS);
        SiteDetails.EndDate(j) = max(PDs.DTS);
        Diffs = diff(PDs.COUNTID);
        AvgDiff = mean(Diffs(Diffs > 0 & Diffs < 10*mean(Diffs(Diffs>0))));
        SiteDetails.TrRateEst(j) = 1/AvgDiff;
        
        % Save 1) see below
        save(strcat('WIM\',num2str(SITE)),'PDs','-v7.3')
        fprintf('\n%i Completed\n',SITE);
        
    end
    
    Sites(strcmp(Sites.SName,SName),:) = SiteDetails;
    
    % Save in 3 forms: 
    % 1) SName based, PD
    % 2) SName based without Axles and Wheelbases, PD(:,[1:12 30 31])
    % 3) Station based, PDs
    
    % Save 2)
    %save(strcat('WIM\',SName),'PD','-v7.3') 
    %PD = PD(:,[1:12 30 31]);
    % Save 3)
    %save(strcat('WIM\',SName,'R'),'PD','-v7.3')
    
end

save('Misc\Sites','Sites')
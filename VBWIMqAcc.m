% ------------------------------------------------------------------------
%                            VBWIMqAcc
% ------------------------------------------------------------------------
% See what accompanies other vehicles when they are in the worst IL spot

% Initializing commands
clear, clc, tic, format long g, load('Sites.mat'), rng('shuffle'), close all;

% Specified Accompaniment Target
Tar = [23]; LaneOps = [1]; Mult = 100; SampleSize = 10; % Must be 1 if Tar = [41,48];

% Read Input File
FName = 'Input/VBWIMqInput60t.xlsx';
BaseData = VBReadInputFile(FName);

% Initialize parpool if necessary and initialize progress bar
if BaseData.Parallel(1) > 0, gcp; clc; end

% Each row of BaseData represents one analysis OR analysis 'set'/'group'
for g = 1:height(BaseData)
    
    % Initialize variables and start row counter
    MaxEvents = []; RamUsed = []; LenPrint = [];
    LostTrucks(g) = 0;
    
    % Recognize if BaseData.SITE(g) is actually a 'set'
    Sitesx = VBGetSiteSet(BaseData.SITE(g),BaseData.LightVehs(g),0,BaseData.Country(g));
        
    for w = 1:length(Sitesx)
        
        % Modify BaseData.SITE(g) based on SiteSet
        BaseData.SITE(g) = Sitesx(w);
        % Update analysis data for current row of BaseData
        [Num,Lane,ILData,~,~] = VBUpdateData(BaseData(g,:));
        ESIA = E.ESPTR;
        
        % Get MaxLength for Spacesave
        MaxLength = (max(arrayfun(@(x) size(x.v,1),ILData))-1)*BaseData.ILRes(g);
        
        % Load File
        load(['WIM/',num2str(BaseData.SITE(g)),'.mat']);
        
        if BaseData.Stage2P(g)
            PDs = Stage2Prune(PDs);
        end

        if Sites.Layout(Sites.SITE == BaseData.SITE(g)) == 11
            % Get Duplicates
            PDs = FindDup2(PDs,0,0);
            % Delete Duplicates - from L1
            PDs(PDs.Dup & PDs.LANE == 1,:) = [];
        end
        
        % Get Only the ClassType Specified
        if strcmp(BaseData.ClassType(g),'Class')
            PDs = PDs(PDs.CLASS > 0 & (PDs.CLASS > 90 | PDs.CLASS < 40),:);
        elseif strcmp(BaseData.ClassType(g),'ClassOW')
            PDs = PDs(PDs.CLASS > 0,:);
        end
                          
        % Separate for each year...
        if ismember('Year',BaseData.Properties.VariableNames)
            UYears = BaseData.Year(g);
        else
            UYears = unique(year(PDs.DTS));
        end
        
        % Start Progress Bar
        u = StartProgBar(length(UYears), 1, str2num([num2str(g) '.' num2str(w)]), str2num([num2str(height(BaseData)) '.' num2str(length(Sitesx))])); tic; st = now;
        
        for r = 1:length(UYears)
            
            PDsy = PDs(year(PDs.DTS) == UYears(r),:);
            
            % For 23 analyse
            % Sample size is how many 23 compare to 41 and 48 we want to keep for analysis. 1 ==> same sample size as 41+48 *samplesize
            Numb23 = sum(sum(PDsy.CLASS == Tar & PDsy.LANE == LaneOps,2)); % Number of 23 actually inside the traffic sample
            NumbAna23 = sum(sum(PDsy.CLASS == [41,48],2)).*SampleSize; % Number of 23 that should be kept for the analysis.
            if Numb23 < NumbAna23                
                LostTrucks(g) = LostTrucks(g) + NumbAna23-Numb23;
                NumbAna23 = Numb23;
            end
            TargIndsFull = find(PDsy.CLASS == Tar & PDsy.LANE == LaneOps);
            TargInds = TargIndsFull(randperm(length(TargIndsFull),NumbAna23));
            Selection = zeros(height(PDsy),1); Selection(TargInds) = 1;
             
            % Modify to only include desired and surrounding 20 vehicles + and - 
            PDsy = PDsy(logical(conv(Selection,[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1],'same')),:);
            
            % For 23 analyse
            % We need to amplify the axle loads of the 23 trucks to be sure that they will be always determinant. To do this we multiply the loads with "Mult"
            
            AXD = strncmp(PDsy.Properties.VariableNames,'AW',2);
            PDsy{PDsy.CLASS == Tar & PDsy.LANE == LaneOps,AXD} = PDsy{PDsy.CLASS == Tar & PDsy.LANE == LaneOps,AXD}*Mult;
            % Don't think we need to modify the GW_TOT... do we?
            PDsy.GW_TOT(PDsy.CLASS == Tar & PDsy.LANE == LaneOps) = PDsy.GW_TOT(PDsy.CLASS == Tar & PDsy.LANE == LaneOps)*Mult;
            
            if isempty(PDsy)
                user = memory;
                RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
                LenPrint = VBUpProgBar(st,g,length(UYears),RamUsed(end),r,LenPrint);
                continue
            end
            
            % Get TrLineUp, AllTrAx, Starti and Endi in sliced form
            [TrLineUpGr,PDsy] = VBGetSlicedPDs2AllTrAx(PDsy,MaxLength,Lane,BaseData.ILRes(g),'day');

            % Perform search for maximums for each day
            %parfor (z = 1:max(PDsy.Group), BaseData.Parallel(g)*100)
            for z = 1:max(PDsy.Group)
                
                % Initialize
                MaxEvents1 = [];
                
                % Get TrLineUpSub and AllTrAxSub
                TrLineUpSub = TrLineUpGr{z};
                % Must sort due to direction switch and accumarray not having negatives for first few
                TrLineUpSub = sortrows(TrLineUpSub);
                
                % Get Lanes
                Lanes = unique(PDsy.LANE);
                
                AllTrAxGr = zeros(max(TrLineUpSub.ATAIndex)-TrLineUpSub.ATAIndex(1)+1,length(Lanes));
                for i = 1:length(Lanes)
                    A = accumarray(TrLineUpSub.ATAIndex(TrLineUpSub.LaneID == Lanes(i))-TrLineUpSub.ATAIndex(1)+1,TrLineUpSub.AxleValue(TrLineUpSub.LaneID == Lanes(i)));
                    AllTrAxGr(1:length(A),i) = A;
                end
                
                % Don't bother running if the segment is too small
                if length(AllTrAxGr) < 450/BaseData.ILRes(g)
                    continue
                    disp('skipped one')
                end % modified temp by lucas true value 2000
                
                % For each InfCase
                for t = 1:Num.InfCases
                    
                    % Get length of bridge in number of indices
                    BrLengthInd = size(ILData(t).v,1);
                    
                    % Reset for each t
                    AllTrAxSub = AllTrAxGr;
                    
                    % Eliminate the need for padding or BrStInd index issues
                    AllTrAxSub(1:BrLengthInd,:) = 0; AllTrAxSub(end-BrLengthInd:end,:) = 0;
                        
                    % Now the issue is that I have to figure out how
                    % much of MaxLE is from the ST... and remove
                    
                    % Subject Influence Line to Truck Axle Stream
                    [MaxLE,DLF,BrStInd,R] = VBGetMaxLE(AllTrAxSub,ILData(t).v(:,1:length(Lanes)),BaseData.RunDyn(g));
                    
                    % Next find a way to save T for every week
                    
                    if MaxLE == 0
                        continue
                    end
                    
                    % Adjust BrStInd for Starti [now TrLineUpSub.ATAIndex(1)]
                    BrStIndx = BrStInd + TrLineUpSub.ATAIndex(1) - 1;
                    % Get BrEndIndx
                    BrEndIndx = BrStIndx + BrLengthInd - 1;
                    % Get Bridge Indices
                    BrIndsx = [BrStIndx:BrEndIndx]';
                    BrInds = [BrStInd:BrStInd+BrLengthInd - 1]';
                    %AxOnBr = sum(AllTrAxSub(BrInds,:),2);
                    %sum(sum(AllTrAxSub(BrInds,:).*flip(ILData(t).v(:,1:length(Lanes)))))
                    
                    % Get Key Info to Save
                    TrNums = TrLineUpSub.TrNum(TrLineUpSub.ATAIndex >= min(BrIndsx) & TrLineUpSub.ATAIndex <= max(BrIndsx));
                    TrLineUpOnBr = TrLineUpSub(TrLineUpSub.ATAIndex >= min(BrIndsx) & TrLineUpSub.ATAIndex <= max(BrIndsx),:);
                       
                    [MaxM, MaxI] = max(TrLineUpOnBr.AxleValue);
                    TrIdMax = TrLineUpOnBr.TrNum(MaxI);
                    
                    TrNumsU = unique(TrNums);
                    
                    % We use TrNums because they don't depend on Starti shift
                    MaxLETime = PDsy.DTS(TrNums(1));
                    Vehs = PDsy.CLASS(TrNumsU);

                    % If there are now 23 vehicles for some reason... it got deleted... skip
                    if sum(sum(Vehs == Tar)) == 0
                        continue
                    end
                    % Get Modified MaxLE. first set alxes from 60t to 0
                    % Find indices of AllTrAxSub... must use TrLineUpOnBr cols 1 and 4 (but w/ x vs no x offset? hmm
                    B4 = sum(sum(AllTrAxSub(BrInds,:).*flip(ILData(t).v(:,1:length(Lanes)))));
                    % Set Vehs == 41 axles to 0
                    MaxLEContr = [];
                    
                    for v = 1:length(TrNumsU)
                        AllTrAxSubTemp = AllTrAxSub;
                        
                        vInds = TrLineUpOnBr.TrNum == TrNumsU(v)';
                        ATXInds = TrLineUpOnBr.ATAIndex(logical(sum(vInds,2)))-(TrLineUpSub.ATAIndex(1)-1);
                        AllTrAxSubTemp(ATXInds,sum((Lanes == (TrLineUpOnBr.LaneID(logical(sum(vInds,2))))').*[1:length(Lanes)]',1)') = 0;
                        
                        % Calculate MaxLE again...
                        MaxLEContr(v) = sum(sum(AllTrAxSubTemp(BrInds,:).*flip(ILData(t).v(:,1:length(Lanes)))));
                    end
                    % MaxLEContr is actually the INVERSE contribution
                    
                    v = find(sum(Vehs == Tar,2));
                    %if max(B4-MaxLEContr(v)) == max(B4-MaxLEContr)
                    if max(B4-MaxLEContr(v)) == max(B4-MaxLEContr)
                        [~,vposi] = max(B4-MaxLEContr(v)); % Check to see if we have more than one target vehicle - if so, choose the one with the larger contribution
                        % values for the 23 must be rescaled with Mult
                        AllTrAxSub(TrLineUpOnBr.ATAIndex(logical(sum(TrLineUpOnBr.TrNum == (TrNumsU(v))',2)))-(TrLineUpSub.ATAIndex(1)-1),sum((Lanes == (TrLineUpOnBr.LaneID(logical(sum((TrLineUpOnBr.TrNum == TrNumsU(v)'),2))))').*[1:length(Lanes)]',1)') = AllTrAxSub(TrLineUpOnBr.ATAIndex(logical(sum(TrLineUpOnBr.TrNum == (TrNumsU(v))',2)))-(TrLineUpSub.ATAIndex(1)-1),sum((Lanes == (TrLineUpOnBr.LaneID(logical(sum((TrLineUpOnBr.TrNum == TrNumsU(v)'),2))))').*[1:length(Lanes)]',1)')/Mult;
                        v = v(vposi); % NEW for more than 41
                        AllTrAxSub(TrLineUpOnBr.ATAIndex(logical(sum(TrLineUpOnBr.TrNum == (TrNumsU(v))',2)))-(TrLineUpSub.ATAIndex(1)-1),sum((Lanes == (TrLineUpOnBr.LaneID(logical(sum((TrLineUpOnBr.TrNum == TrNumsU(v)'),2))))').*[1:length(Lanes)]',1)') = 0;
                        MaxLE = sum(sum(AllTrAxSub(BrInds,:).*flip(ILData(t).v(:,1:length(Lanes)))));
                    else
                        continue
                    end
                    
                    % Are we sure about this?
                    if MaxLE < 0
                        MaxLE = 0;
                    end
                    
                    % Get ClassT (in m form for now)
                    if min(Vehs) == 0
                        m = 1;
                    elseif sum(sum(Vehs == Tar,2)) > 0
                        m = 2;
                    else
                        m = 3;
                    end
                    
                    if BaseData.Apercu(g) == 1
                        %if MaxLE > 9000 && m == 3
                        if MaxLE > 0
                            TrLineUpSubTemp = TrLineUpSub;
                            %PosiRemTruck = TrNumsU.*(Vehs==41); PosiRemTruck(PosiRemTruck == 0) = []; PosiRemTruck = PosiRemTruck(1);
                            PosiRemTruck = TrNumsU(v);
                            % return 23 loads to original using Mult
                            TrLineUpSub(TrLineUpSub.TrNum == PosiRemTruck,2) = TrLineUpSub(TrLineUpSub.TrNum == PosiRemTruck,2)/Mult;
                            TrLineUpSubTemp(TrLineUpSub.TrNum == PosiRemTruck,:) = [];
                            T = VBApercuv2(PDsy,'',ILData(t),BrStIndx,table2array(TrLineUpSubTemp),1,Lane,BaseData.ILRes(g));
                            
                            AXD = strncmp(PDsy.Properties.VariableNames,'AW',2);
                            PDsy{PDsy.CLASS == Tar & PDsy.LANE == LaneOps,AXD} = PDsy{PDsy.CLASS == Tar & PDsy.LANE == LaneOps,AXD}/Mult;
                            PDsy.GW_TOT(PDsy.CLASS == Tar & PDsy.LANE == LaneOps) = PDsy.GW_TOT(PDsy.CLASS == Tar & PDsy.LANE == LaneOps)/Mult;
                            
                            Tx = VBApercuv2(PDsy,'',ILData(t),BrStIndx,table2array(TrLineUpSub),1,Lane,BaseData.ILRes(g));
                            
                            %return original values
                            TrLineUpSub(TrLineUpSub.TrNum == PosiRemTruck,2) = TrLineUpSub(TrLineUpSub.TrNum == PosiRemTruck,2)*Mult;
                            PDsy{PDsy.CLASS == Tar & PDsy.LANE == LaneOps,AXD} = PDsy{PDsy.CLASS == Tar & PDsy.LANE == LaneOps,AXD}*Mult;
                            PDsy.GW_TOT(PDsy.CLASS == Tar & PDsy.LANE == LaneOps) = PDsy.GW_TOT(PDsy.CLASS == Tar & PDsy.LANE == LaneOps)*Mult;

                        end
                        %exportgraphics(gcf,"Max"  + ".jpg",'Resolution',600)
                    end
                    
                    % Save MaxEvents... save Times and Datenums and then convert
                    MaxEvents1 = [MaxEvents1; datenum(MaxLETime), BaseData.SITE(g), MaxLE, t, m, 1, BrStIndx];
                        
                end % t, InfCases
                
                MaxEvents = [MaxEvents; MaxEvents1];
            end % z, groups
            
            % Update progress bar
            user = memory;
            RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
            LenPrint = VBUpProgBar(st,RamUsed(end),r,LenPrint);
                        
        end % r, years
    end % w, SiteGroups
  
    % If there was no determinent vehicules, initialise
    if height(MaxEvents) == 0
        MaxEvents = [datenum(MaxLETime), BaseData.SITE(g), MaxLE, t, m, k, BrStIndx];
    end
    
    % Convert back to datetime
    MaxEvents = array2table(MaxEvents,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'BrStInd'});clc
    MaxEvents.Datenum = datetime(MaxEvents.Datenum,'ConvertFrom','datenum');
    MaxEvents = renamevars(MaxEvents,"Datenum","DTS");
           
    TName = datestr(now,'mmmdd-yy HHMMSS');
    % Need to go back to original BaseData... no SITE switch
    BaseData = VBReadInputFile(FName);
    OutInfo.Name = TName; OutInfo.BaseData = BaseData(g,:);
    OutInfo.ILData = ILData;
    OutInfo.MaxEvents = MaxEvents;
    
    % Create folders where there are none
    CreateFolders(BaseData.Folder{g},BaseData.VWIM(g),BaseData.Apercu(g),BaseData.Save(g))
    
    if BaseData.Save(g) == 1
        save(['Output' BaseData.Folder{g} '/' OutInfo.Name], 'OutInfo')
    end    
end % g, BaseData

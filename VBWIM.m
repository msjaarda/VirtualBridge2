% ------------------------------------------------------------------------
%                             VBWIM
% ------------------------------------------------------------------------
% Generate a summary of maximum effects on bridges from real WIM
% Equipped with:
%    - StopSim
%    - Fatigue
%    - Platoon Swapping

% Initializing commands
clear, clc, tic, format long g, load('Sites.mat'), rng('shuffle'), close all;

% Read Input File
FName = 'Input/VBWIMInputF325.xlsx';
BaseData = VBReadInputFile(FName);

% Initialize parpool if necessary and initialize progress bar
if BaseData.Parallel(1) > 0, gcp; clc; end

% Each row of BaseData represents one analysis OR analysis 'set'/'group'
for g = 1:height(BaseData)
    
    % Initialize variables
    MaxEvents = []; RamUsed = []; LenPrint = []; MaxEventsStop = [];
    
    % Recognize if BaseData.SITE(g) is actually a set
    SitesList = VBGetSiteSet(BaseData.SITE(g),BaseData.LightVehs(g),0,BaseData.Country(g));
    
    % For each Site    
    for w = 1:length(SitesList)
        
        % Modify BaseData.SITE(g) based on SiteSet
        BaseData.SITE(g) = SitesList(w);
        % Update analysis data for current row of BaseData
        [Num,Lane,ILData,~,~] = VBUpdateData(BaseData(g,:));
        
        % Get MaxLength for Spacesave
        MaxLength = (max(arrayfun(@(x) size(x.v,1),ILData))-1)*BaseData.ILRes(g);
        
        % Load File
        load(['WIM/',num2str(BaseData.SITE(g)),'.mat']);
        
        % Stage 2 Prune
        if BaseData.Stage2P(g); PDs = Stage2Prune(PDs); end
        
        % Find and remove duplicates
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
        
        % Platooning Swapping
        if BaseData.Plat(g)
            VType = VBVTypes; TrTyps = VType.CLASS(VType.Class == 1);
            PlatPct = BaseData.PlatRate(g)*ones(length(TrTyps),1); % PlatPct(TrTyps == 22 | TrTyps == 23) = 0;
            PDs = SwapforPlatoonsWIM(PDs,BaseData.PlatSize(g),BaseData.PlatFolDist(g),TrTyps,PlatPct);
        end
        
        % Get unique years
        UYears = unique(year(PDs.DTS));
        if ismember('Year',BaseData.Properties.VariableNames); UYears = BaseData.Year(g); end
        
        % Clear variable to save space if light vehicles are included
        if Sites.Light(Sites.SITE == BaseData.SITE(g)) && ~BaseData.Stage2P(g); clear PDs; end
        
        % Start Progress Bar
        u = StartProgBar(length(UYears), 1, str2num([num2str(g) '.' num2str(w)]), str2num([num2str(height(BaseData)) '.' num2str(length(SitesList))])); tic; st = now;
        
        % For each unique year
        for r = 1:length(UYears)
            
            if Sites.Light(Sites.SITE == BaseData.SITE(g)) && ~BaseData.Stage2P(g)
                % Save memory by clearing PDs and re-loading for Light Vehs
                [PDsy] = LoadPDYear(['WIM/',num2str(BaseData.SITE(g)),'.mat'],UYears(r));
            else
                PDsy = PDs(year(PDs.DTS) == UYears(r),:);
            end
                       
            % Get TrLineUp, AllTrAx, Starti and Endi in sliced form
            [TrLineUpGr,PDsy] = VBGetSlicedPDs2AllTrAx(PDsy,MaxLength,Lane,BaseData.ILRes(g),'day');

            % Perform search for maximums for each day
            parfor (z = 1:max(PDsy.Group), BaseData.Parallel(g)*100)
            %for z = 1:max(PDsy.Group)
                
                % Initialize (required for parallelization)
                MaxEvents1 = []; MaxEvents1Stop = [];
                
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
                if length(AllTrAxGr) < 2000/BaseData.ILRes(g), continue, end
                
                % For each InfCase
                for t = 1:Num.InfCases
                    
                    % Get length of bridge in number of indices
                    BrLengthInd = size(ILData(t).v,1);
                    
                    % Reset for each t
                    AllTrAxSub = AllTrAxGr;
                    
                    % Eliminate the need for padding or BrStInd index issues
                    AllTrAxSub(1:BrLengthInd,:) = 0; AllTrAxSub(end-BrLengthInd:end,:) = 0;
                    
                    % For each analysis
                    k = 0; % Initialize k
                    while k < BaseData.NumAnalyses(g) && sum(AllTrAxSub,'all') > 0
                        
                        % Subject Influence Line to Truck Axle Stream
                        [MaxLE,DLF,BrStInd,R] = VBGetMaxLE(AllTrAxSub,ILData(t).v(:,1:length(Lanes)),BaseData.RunDyn(g));
                        
                        % Get rainflow for fatigue
                        if BaseData.RunFat(g) == 1
                            [c,hist,edges,rmm,idx] = rainflow(R);
                            T = array2table(c,'VariableNames',{'Count','Range','Mean','Start','End'});
                        end
                        % Next find a way to save T for every week
                        
                        if MaxLE == 0, k = k+1; continue, end
                        k = k+1; % Add to k
                        
                        % Adjust BrStInd for Starti [now TrLineUpSub.ATAIndex(1)]
                        BrStIndx = BrStInd + TrLineUpSub.ATAIndex(1) - 1;
                        % Get BrEndIndx
                        BrEndIndx = BrStIndx + BrLengthInd - 1;
                        % Get Bridge Indices
                        BrIndsx = [BrStIndx:BrEndIndx]';
                        BrInds = [BrStInd:BrStInd+BrLengthInd - 1]';
                        %AxOnBr = sum(AllTrAxt(StripInds,:),2);
                        
                        % Get Key Info to Save
                        TrNums = TrLineUpSub.TrNum(TrLineUpSub.ATAIndex >= min(BrIndsx) & TrLineUpSub.ATAIndex <= max(BrIndsx));
                        TrLineUpOnBr = TrLineUpSub(TrLineUpSub.ATAIndex >= min(BrIndsx) & TrLineUpSub.ATAIndex <= max(BrIndsx),:);
                        [MaxM, MaxI] = max(TrLineUpOnBr.AxleValue);
                        TrIdMax = TrLineUpOnBr.TrNum(MaxI);
                        
                        % Get Unique Truck Identifiers
                        TrNumsU = unique(TrNums);
                        
                        % Perform StopSim if necessary
                        if BaseData.StopSim(g)
                            NumExtra = 20;
                            TrNumsUE = [max(1,TrNumsU(1)-NumExtra):min(TrNumsU(end)+NumExtra,height(PDsy))]';
                            
                            PDe = PDsy(TrNumsUE,:);
                            
                            % Call VBWIMtoAllTrAx w/ mods... must give stationary point or truck
                            [PDe, AllTrAxStop, TrLineUpStop] = VBWIMtoAllTrAxStop(PDe,MaxLength,Lane,BaseData.ILRes(g),find(TrNumsUE == TrIdMax));
                            TrLineUpStop = array2table(TrLineUpStop,'VariableNames',{'ATAIndex','AxleValue','TrNum','LaneID'});
                            
                            % Round TrLineUp to ILRes
                            TrLineUpStop.mStation = TrLineUpStop.ATAIndex;
                            TrLineUpStop.ATAIndex = round(TrLineUpStop.ATAIndex/BaseData.ILRes(g));
                            
                            [MaxLEe,DLFe,BrStInde,Re] = VBGetMaxLE(AllTrAxStop,ILData(t).v,BaseData.RunDyn(g));
                        end
                        
                        % We use TrNums because they don't depend on Starti shift
                        MaxLETime = PDsy.DTS(TrNums(1));
                        Vehs = PDsy.CLASS(TrNumsU);
                        
                        % Get ClassT (in m form for now)
                        m = 3;
                        if min(Vehs) == 0
                            m = 1;
                        elseif sum(Vehs > 39 & Vehs < 90) > 0
                            m = 2;
                        end
                        
                        if BaseData.Apercu(g) == 1
                            T = VBApercuv2(PDsy,'',ILData(t),BrStIndx,table2array(TrLineUpSub),1,Lane,BaseData.ILRes(g)); %exportgraphics(gcf,"Max"  + ".jpg",'Resolution',600);
                            if BaseData.StopSim(g)
                                TStop = VBApercuv2(PDe,'',ILData(t),BrStInde,table2array(TrLineUpStop),1,Lane,BaseData.ILRes(g));
                            end
                        end
                        
                        % Save MaxEvents... save Times and Datenums and then convert
                        if BaseData.Plat(g) % Also return platoon type
                            if max(PDsy.Plat(TrNumsU)) > 0
                                MaxEvents1 = [MaxEvents1; datenum(MaxLETime), BaseData.SITE(g), MaxLE, t, m, k, BrStIndx mode(Vehs)];
                            else
                                MaxEvents1 = [MaxEvents1; datenum(MaxLETime), BaseData.SITE(g), MaxLE, t, m, k, BrStIndx 0];
                            end
                        else
                            MaxEvents1 = [MaxEvents1; datenum(MaxLETime), BaseData.SITE(g), MaxLE, t, m, k, BrStIndx];
                        end
                        
                        % Save for StopSim
                        if BaseData.StopSim(g)
                            MaxEvents1Stop = [MaxEvents1Stop; datenum(MaxLETime), BaseData.SITE(g), MaxLEe, t, m, k, BrStInde];
                        end
                        
                        % Bump k up so that analysis doesn't continue!
                        if m == 3; break; end
                        
                        % Prepare for next run - Set Axles to zero in AllTrAx (can't delete because indices are locations)
                        AllTrAxSub(BrInds,:) = 0;
                        
                    end % k, analyses
                end % t, InfCases
                
                % Save MaxEvetns to master
                MaxEvents = [MaxEvents; MaxEvents1];
                if BaseData.StopSim(g)
                    MaxEventsStop = [MaxEventsStop; MaxEvents1Stop];
                end
            end % z, groups
            
            % Update progress bar
            user = memory;
            RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
            LenPrint = VBUpProgBar(st,RamUsed(end),r,LenPrint);
            
        end % r, years
    end % w, SiteGroups
  
    % Convert back to datetime (building the array required datenum)
    if BaseData.Plat(g)
        MaxEvents = array2table(MaxEvents,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'BrStInd', 'PlatType'});
    else
        MaxEvents = array2table(MaxEvents,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'BrStInd'});
    end
    MaxEvents.Datenum = datetime(MaxEvents.Datenum,'ConvertFrom','datenum'); MaxEvents = renamevars(MaxEvents,"Datenum","DTS");
    if BaseData.StopSim(g)
        MaxEventsStop = array2table(MaxEventsStop,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'BrStInd'});
        MaxEventsStop.Datenum = datetime(MaxEventsStop.Datenum,'ConvertFrom','datenum'); MaxEventsStop = renamevars(MaxEventsStop,"Datenum","DTS");
    end
           
    TName = datestr(now,'mmmdd-yy HHMMSS');
    % Need to go back to original BaseData... no SITE switch
    BaseData = VBReadInputFile(FName);
    OutInfo.Name = TName;
    OutInfo.BaseData = BaseData(g,:);
    OutInfo.BaseData.StopSim = 0;
    % We do not need ESIA... can load this from outside later on
    %OutInfo.ESIA = ESIA;
    OutInfo.ILData = ILData;
    OutInfo.MaxEvents = MaxEvents;
    
    % Create folders where there are none
    CreateFolders(BaseData.Folder{g},BaseData.VWIM(g),BaseData.Apercu(g),BaseData.Save(g))
    
    if BaseData.Save(g) == 1
        save(['Output' BaseData.Folder{g} '/' OutInfo.Name], 'OutInfo','-v7.3')
    end
    
    if BaseData.StopSim(g)
        OutInfo.BaseData.StopSim = 1;
        % Hmmm should this line be there? I don't think so... comment 4 now
        %MaxEventsStop(MaxEventsStop.MaxLE <= 0,:) = [];
        
        % Add some time to not have the same name as non StopSim
        TName = datestr(now+1/86400,'mmmdd-yy HHMMSS');
        OutInfo.Name = TName;
        OutInfo.MaxEvents = MaxEventsStop;
        
        if BaseData.Save(g) == 1
            save(['Output' BaseData.Folder{g} '/' OutInfo.Name], 'OutInfo','-v7.3')
        end
    end
    
end % g, BaseData

% ------------------------------------------------------------------------
%                            VBWIMq
% ------------------------------------------------------------------------
% Generate a summary of maximum effects on bridges from real WIM

% AxleStatsBasic >> AxTandem      >> Q1   Investigation
% VBWIMQ1Q2      >> Q1Q2MaxEvents >> Q1Q2 Investigation
% VBWIMq         >> qMaxEvents    >> q    Investigation

% Initializing commands
clear, clc, tic, format long g, rng('shuffle'), close all;

% Read Input File
FName = 'Input/VBWIMqInputxxx.xlsx';
BaseData = VBReadInputFile(FName);

% Initialize parpool if necessary and initialize progress bar
if BaseData.Parallel(1) > 0, gcp; clc; end

% Each row of BaseData represents one analysis OR analysis 'set'/'group'
for g = 1:height(BaseData)
    
    % Initialize variables and start row counter
    MaxEvents = []; RamUsed = []; LenPrint = []; MaxEventsStop = [];
    
    % Recognize if BaseData.SITE(g) is actually a 'set'
    load('SiteGroups.mat')
    if BaseData.SITE(g) == 11, Sites = SiteGroups.('Uni2L');
    elseif BaseData.SITE(g) == 111, Sites = SiteGroups.('Uni3L');
    elseif BaseData.SITE(g) == 12, Sites = SiteGroups.('Bi2L');
    elseif BaseData.SITE(g) == 1122, Sites = SiteGroups.('Bi4L');
    elseif BaseData.SITE == 110, Traffic = SiteGroups.('LSVAUni2L');
    else Sites = BaseData.SITE(g);
    end
        
    for w = 1:length(Sites)
        
        % Modify BaseData.SITE(g) based on SiteSet
        BaseData.SITE(g) = Sites(w);
        % Update analysis data for current row of BaseData
        [Num,Lane,ILData,~,~,ESIA] = VBUpdateData(BaseData(g,:));
        
        % Get MaxLength for Spacesave
        MaxLength = (max(arrayfun(@(x) size(x.v,1),ILData))-1)*BaseData.ILRes(g);
        
        % Load File
        if BaseData.LSVA(g)
            load(['WIMLSVA/',num2str(BaseData.SITE(g)),'.mat']);
        else
            load(['WIM/',num2str(BaseData.SITE(g)),'.mat']);
        end
        
        if BaseData.Stage2P(g)
            PDs = Stage2Prune(PDs);
        end
        
        % Get Only the ClassType Specified
        try
            if strcmp(BaseData.ClassType(g),'Class')
                PDs = PDs(PDs.CLASS > 0 & (PDs.CLASS > 50 | PDs.CLASS < 40),:);
            elseif strcmp(BaseData.ClassType(g),'ClassOW')
                PDs = PDs(PDs.CLASS > 0,:);
            end
        catch 
        end
        
        % Separate for each year...
        if ismember('Year',BaseData.Properties.VariableNames)
            UYears = BaseData.Year(g);
        else
            UYears = unique(year(PDs.DTS));
        end
        
        if BaseData.LSVA(g) && ~BaseData.Stage2P(g)
            clear PDs
        end
        
        % Start Progress Bar
         u = StartProgBar(length(UYears), 1, g, height(BaseData)); tic; st = now;
        
        for r = 1:length(UYears)
            
            if BaseData.LSVA(g) && ~BaseData.Stage2P(g)
                % Try to save mem by clearing PDs and re-loading
                [PDsy] = LoadPDYear(['WIMLSVA/',num2str(BaseData.SITE(g)),'.mat'],UYears(r));
            else
                PDsy = PDs(year(PDs.DTS) == UYears(r),:);
            end
                       
            % Get TrLineUp, AllTrAx, Starti and Endi in sliced form
            [TrLineUpGr,PDsy] = GetSlicedPDs2AllTrAx(PDsy,MaxLength,Lane,BaseData.ILRes(g));

            % Perform search for maximums for each day
            parfor (z = 1:max(PDsy.Group), BaseData.Parallel(g)*100)
            %for z = 1:max(PDsy.Group)
                
                % Initialize
                MaxEvents1 = []; MaxEvents1Stop = [];
                
                % Get TrLineUpSub and AllTrAxSub
                TrLineUpSub = TrLineUpGr{z};
                % Must sort due to direction switch and accumarray not having negatives for first few
                TrLineUpSub = sortrows(TrLineUpSub);
                
                % Get Lanes
                Lanes = unique(PDsy.LANE);
                AllTrAxGr = zeros(max(TrLineUpSub(:,1))-TrLineUpSub(1,1)+1,length(Lanes));
                for i = 1:length(Lanes)
                    A = accumarray(TrLineUpSub(TrLineUpSub(:,4)==Lanes(i),1)-TrLineUpSub(1,1)+1,TrLineUpSub(TrLineUpSub(:,4)==Lanes(i),2));
                    AllTrAxGr(1:length(A),i) = A;
                end
                
                % Don't bother running if the segment is too small
                if length(AllTrAxGr) < 20000/BaseData.ILRes(g), continue, end
                
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
                        [MaxLE,DLF,BrStInd,~] = VBGetMaxLE(AllTrAxSub,ILData(t).v(:,1:length(Lanes)),BaseData.RunDyn(g));
                        if MaxLE == 0, k = k+1; continue, end
                        k = k+1; % Add to k
                        
                        % Adjust BrStInd for Starti
                        BrStIndx = BrStInd + TrLineUpSub(1,1) -1;
                        % Get BrEndIndx
                        BrEndIndx = BrStIndx + BrLengthInd - 1;
                        % Get Bridge Indices
                        BrIndsx = [BrStIndx:BrEndIndx]';
                        BrInds = [BrStInd:BrStInd+BrLengthInd - 1]';
                        %AxOnBr = sum(AllTrAxt(StripInds,:),2);
                        
                        % Get Key Info to Save
                        TrNums = TrLineUpSub(TrLineUpSub(:,1) >= min(BrIndsx) & TrLineUpSub(:,1) <= max(BrIndsx),3);
                        TrLineUpOnBr = TrLineUpSub(TrLineUpSub(:,1) >= min(BrIndsx) & TrLineUpSub(:,1) <= max(BrIndsx),:);
                        [MaxM, MaxI] = max(TrLineUpOnBr(:,2));
                        TrIdMax = TrLineUpOnBr(MaxI,3);
                        
                        TrNumsU = unique(TrNums);
                        
                        if BaseData.StopSim(g)
                            NumExtra = 20;
                            TrNumsUE = [max(1,TrNumsU(1)-NumExtra):min(TrNumsU(end)+NumExtra,height(PDsy))]';
                            
                            PDe = PDsy(TrNumsUE,:);
                            
                            % Call VBWIMtoAllTrAx w/ mods... must give stationary point or truck
                            [PDe, AllTrAxStop, TrLineUpStop] = VBWIMtoAllTrAxStop(PDe,MaxLength,Lane,BaseData.ILRes(g),find(TrNumsUE == TrIdMax));
                            
                            % Round TrLineUp first row, move unrounded to fifth row
                            TrLineUpStop(:,5) = TrLineUpStop(:,1); TrLineUpStop(:,1) = round(TrLineUpStop(:,1)/BaseData.ILRes(g));
                            % TrLineUpStop [ 1: AllTrAxIndex  2: AxleValue  3: Truck#  4: LaneID  5: Station(m) ]
                            
                            [MaxLEe,DLFe,BrStInde,Re] = VBGetMaxLE(AllTrAxStop,ILData(t).v,BaseData.RunDyn(g));
                        end
                        
                        % We use TrNums because they don't depend on Starti shift
                        MaxLETime = PDsy.DTS(TrNums(1));
                        Vehs = PDsy.CLASS(TrNumsU);
                        
                        if BaseData.Apercu(g) == 1
                            T = VBApercu(PDsy,'',ILData(t),BrStIndx,TrLineUpSub,MaxLE/ESIA.Total(t),1,Lane,BaseData.ILRes(g));
                            %exportgraphics(gcf,"Apercu" + BaseData.Folder(g) + "/" + ApercuTitle + ".jpg",'Resolution',600)
                            if BaseData.StopSim(g)
                                TStop = VBApercu(PDe,'',ILData(t),BrStInde,TrLineUpStop,MaxLEe/ESIA.Total(t),1,Lane,BaseData.ILRes(g));
                            end
                        end
                        
                        % Only collect detailed info if desired... function
                        % not written right now
                        %[L1Veh,L2Veh,L1Spd,L2Spd,L1Load,L2Load,L1Ax,L2Ax] = DetailedVBWIM(PDsy,TrNumsU,Vehs,AllTrAxSub,BrInds,Starti);
                        
                        % Get ClassT (in m form for now)
                        if min(Vehs) == 0
                            m = 1;
                        elseif sum(Vehs > 39 & Vehs < 50) > 0
                            m = 2;
                        else
                            m = 3;
                        end
                        
                        % Save MaxEvents... save Times and Datenums and then convert
                        MaxEvents1 = [MaxEvents1; datenum(MaxLETime), BaseData.SITE(g), MaxLE, t, m, k, BrStIndx];
                        % Rewrite line if DetailedVBWIM Desired
                        if BaseData.StopSim(g)
                            MaxEvents1Stop = [MaxEvents1Stop; datenum(MaxLETime), BaseData.SITE(g), MaxLEe, t, m, k, BrStInde];
                        end
                        
                        if m == 3
                            k = 100; % Bump k up so that analysis doesn't continue!
                        end
                        
                        % Prepare for next run - Set Axles to zero in AllTrAx (can't delete because indices are locations)
                        AllTrAxSub(BrInds,:) = 0;
                        
                    end % k, analyses
                end % t, InfCases
                
                MaxEvents = [MaxEvents; MaxEvents1];
                if BaseData.StopSim(g)
                    MaxEventsStop = [MaxEventsStop; MaxEvents1Stop];
                end
            end % z, groups
            
            % Update progress bar
            user = memory;
            RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
            LenPrint = VBUpProgBar(u,st,g,1,length(UYears),1,RamUsed(end),r,LenPrint);
            
        end % r, years
    end % w, SiteGroups
  
    % Convert back to datetime
    MaxEvents = array2table(MaxEvents,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'BrStInd'});
    % Below for Details VBWIM
    %MaxEvents = array2table(MaxEvents,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'L1Veh', 'L2Veh', 'L1Load', 'L2Load', 'L1Ax', 'L2Ax', 'L1Sp', 'L2Sp'});
    MaxEvents.DTS = datetime(MaxEvents.Datenum,'ConvertFrom','datenum');
    MaxEvents.Datenum = [];
    
    if BaseData.StopSim(g)
        MaxEventsStop = array2table(MaxEventsStop,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'BrStInd'});
        MaxEventsStop.DTS = datetime(MaxEventsStop.Datenum,'ConvertFrom','datenum');
        MaxEventsStop.Datenum = [];
    end
    
    % m = 1 is ClassT 'All', m = 2 is 'ClassOW', and m = 3 is 'Class'
    % qInvestInitial Inputs
    BM = {'Daily', 'Weekly', 'Yearly'};             % j
    ClassType = {'All', 'ClassOW', 'Class'};        % i
    DistTypes = {'Lognormal'};
    [Max,~,~,~] = qInvestInitial(BM,ClassType,DistTypes,MaxEvents,ILData);
    
    TName = datestr(now,'mmmdd-yy HHMMSS');
    % Need to go back to original BaseData... no SITE switch
    BaseData = VBReadInputFile(FName);
    OutInfo.Name = TName; OutInfo.BaseData = BaseData(g,:);
    OutInfo.ESIA = ESIA;
    OutInfo.ILData = ILData;
    OutInfo.SimStop = false;
    OutInfo.Max = Max;
    
    % Plot BlockMax, find Design Values, Ed, using Beta, rather than 99th percentile
    for r = 1:Num.InfCases
        for i = 1:length(ClassType)
            Class = ClassType{i};
            BlockM = BM{2};
            [~,OutInfo.x_values.(Class).(BlockM)(:,r),OutInfo.y_valuespdf.(Class).(BlockM)(:,r),~] = GetBlockMaxFit(Max(r).(Class).(BlockM).Max,'Lognormal',BaseData.Plots(g));
            %[ECDF,ECDFRank,PPx,PPy,Fity,OutInfo.LNFitR2] = GetLogNormPPP(Max(r).(Class).(BlockM).Max,false);
            [OutInfo.EdLN.(Class).(BlockM)(r), OutInfo.AQ.(Class).(BlockM)(r), ~] = GetBlockMaxEd(Max(r).(Class).(BlockM).Max,BlockM,'Lognormal',ESIA.Total(r),ESIA.EQ(:,r),ESIA.Eq(:,r),0.6,0.5);
        end
    end
    
    % Create folders where there are none
    CreateFolders(BaseData.Folder{g},BaseData.VWIM(g),BaseData.Apercu(g),BaseData.Save(g))
    
    if BaseData.Save(g) == 1
        save(['Output' BaseData.Folder{g} '/' OutInfo.Name], 'OutInfo')
    end
    
    if BaseData.StopSim(g)
        MaxEventsStop(MaxEventsStop.MaxLE <= 0,:) = [];
        [Max,~,~,~] = qInvestInitial(BM,ClassType,DistTypes,MaxEventsStop,ILData);
        
        TName = datestr(now+1/86400,'mmmdd-yy HHMMSS');
        OutInfo.Name = TName;
        OutInfo.SimStop = true;
        OutInfo.Max = Max;
        
        % Plot BlockMax, find Design Values, Ed, using Beta, rather than 99th percentile
        for r = 1:Num.InfCases
            for i = 1:length(ClassType)
                Class = ClassType{i};
                BlockM = BM{2};
                %BlockM = BM{1};
                [~,OutInfo.x_values.(Class).(BlockM)(:,r),OutInfo.y_valuespdf.(Class).(BlockM)(:,r),~] = GetBlockMaxFit(Max(r).(Class).(BlockM).Max,'Lognormal',BaseData.Plots(g));
                %[ECDF,ECDFRank,PPx,PPy,Fity,OutInfo.LNFitR2] = GetLogNormPPP(Max(r).(Class).(BlockM).Max,false);
                [OutInfo.EdLN.(Class).(BlockM)(r), OutInfo.AQ.(Class).(BlockM)(r), ~] = GetBlockMaxEd(Max(r).(Class).(BlockM).Max,BlockM,'Lognormal',ESIA.Total(r),ESIA.EQ(:,r),ESIA.Eq(:,r),0.6,0.5);
            end
        end
        
        if BaseData.Save(g) == 1
            save(['Output' BaseData.Folder{g} '/' OutInfo.Name], 'OutInfo')
        end
    end
    
end % g, BaseData

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
FName = 'Input/VBWIMqInput60t-test.xlsx';
BaseData = VBReadInputFile(FName);

% Let's try to delete all the WIM records not around the 60t vehicles...
% That way we specifically find just the 60t. Try this with 60t... see 115

% Initialize parpool if necessary and initialize progress bar
if BaseData.Parallel(1) > 0, gcp; clc; end

% Each row of BaseData represents one analysis OR analysis 'set'/'group'
for g = 1:height(BaseData)
    
    % Initialize variables and start row counter
    MaxEvents = []; RamUsed = []; LenPrint = []; MaxEventsStop = []; load('SiteGroups')
    LostTrucks(g) = 0;
    
    % Recognize if BaseData.SITE(g) is actually a 'set'
    Sites = VBGetSiteSet(BaseData.SITE(g),BaseData.StopSim(g));
        
    for w = 1:length(Sites)
        
        % Modify BaseData.SITE(g) based on SiteSet
        BaseData.SITE(g) = Sites(w);
        % Update analysis data for current row of BaseData
        [Num,Lane,ILData,~,~,E] = VBUpdateData(BaseData(g,:));
        ESIA = E.ESPTR;
        
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
        
        if sum(BaseData.SITE(g)== SiteGroups.Bi4L) == 0 && sum(BaseData.SITE(g)== SiteGroups.Uni3L) == 0
            % Get Duplicates
            PDs = FindDup2(PDs,0,0);
            % Delete Duplicates - from L1
            PDs(PDs.Dup & PDs.LANE == 2,:) = [];
        end
        
        % Get Only the ClassType Specified
        try
            if strcmp(BaseData.ClassType(g),'Class')
                PDs = PDs(PDs.CLASS > 0 & (PDs.CLASS > 90 | PDs.CLASS < 40),:);
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
            
            % For 23 analyse
            SampleSize = 1; % How many 23 compare to 41 and 48 we want to keep for analysis. 1 ==> same sample size as 41+48 *samplesize
            Numb23 = sum(sum(PDsy.CLASS == [23] & PDsy.LANE == [1],2)); % Number of 23 actually inside the traffic sample
            NumbAna23 = sum(sum(PDsy.CLASS == [41,48],2)).*SampleSize; % Number of 23 that should be kept for the analysis.
            if Numb23 < NumbAna23                
            LostTrucks(g) = LostTrucks(g) + NumbAna23-Numb23;
            NumbAna23 = Numb23;
            end
            Rand23 = [ones(Numb23-NumbAna23,1); zeros(NumbAna23,1)];
            Rand23 = Rand23(randperm(length(Rand23)));
            Rand23 = find(PDsy.CLASS==23 & PDsy.LANE == [1]).*Rand23;
            Rand23(Rand23 == 0)= [];
            PDsy.CLASS(Rand23) = 0; % Set non keeped 23 trucks as 0 trucks type
            
            % Modify to only include type 41 and surrounding vehicles
            %PDsy = PDsy(logical(conv(PDsy.CLASS == 41,[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1],'same')),:);                      
            %PDsy = PDsy(logical(conv(sum(PDsy.CLASS == [23],2),[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1],'same')),:);
            % Add 20 vehicles + and -   
            PDsy = PDsy(logical(conv(sum(PDsy.CLASS == [23] & PDsy.LANE == [1],2),[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1],'same')),:);
            % Add 6 vehicles + and -   
            %PDsy(((PDsy.LANE == 2 | PDsy.LANE == 3)&sum(PDsy.CLASS == [41,48],2)),:) = [];
            % Removing special transport from 2/3 lanes
            
            % For 23 analyse
            % We need to amplify the axle loads of the 23 trucks to be sure
            % that they will be always determinant. To do this we multiply the loads
            % with "Multiplicater"
            Multiplicater = 100;
            PDsy.AWT01(PDsy.CLASS == [23] & PDsy.LANE == [1]) = PDsy.AWT01(PDsy.CLASS == [23] & PDsy.LANE == [1])*Multiplicater;
            PDsy.AWT02(PDsy.CLASS == [23] & PDsy.LANE == [1]) = PDsy.AWT02(PDsy.CLASS == [23] & PDsy.LANE == [1])*Multiplicater;
            PDsy.AWT03(PDsy.CLASS == [23] & PDsy.LANE == [1]) = PDsy.AWT03(PDsy.CLASS == [23] & PDsy.LANE == [1])*Multiplicater;
            PDsy.AWT04(PDsy.CLASS == [23] & PDsy.LANE == [1]) = PDsy.AWT04(PDsy.CLASS == [23] & PDsy.LANE == [1])*Multiplicater;
            PDsy.AWT05(PDsy.CLASS == [23] & PDsy.LANE == [1]) = PDsy.AWT05(PDsy.CLASS == [23] & PDsy.LANE == [1])*Multiplicater;
            PDsy.AWT06(PDsy.CLASS == [23] & PDsy.LANE == [1]) = PDsy.AWT06(PDsy.CLASS == [23] & PDsy.LANE == [1])*Multiplicater;
            PDsy.AWT07(PDsy.CLASS == [23] & PDsy.LANE == [1]) = PDsy.AWT07(PDsy.CLASS == [23] & PDsy.LANE == [1])*Multiplicater;
            PDsy.AWT08(PDsy.CLASS == [23] & PDsy.LANE == [1]) = PDsy.AWT08(PDsy.CLASS == [23] & PDsy.LANE == [1])*Multiplicater;
            PDsy.AWT09(PDsy.CLASS == [23] & PDsy.LANE == [1]) = PDsy.AWT09(PDsy.CLASS == [23] & PDsy.LANE == [1])*Multiplicater;
            PDsy.GW_TOT(PDsy.CLASS == [23] & PDsy.LANE == [1]) = PDsy.GW_TOT(PDsy.CLASS == [23] & PDsy.LANE == [1])*Multiplicater;
            
            if isempty(PDsy)
                user = memory;
                RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
                LenPrint = VBUpProgBar(st,g,length(UYears),RamUsed(end),r,LenPrint);
                continue
            end
            
            % Get TrLineUp, AllTrAx, Starti and Endi in sliced form
            [TrLineUpGr,PDsy] = GetSlicedPDs2AllTrAx(PDsy,MaxLength,Lane,BaseData.ILRes(g));

            % Perform search for maximums for each day
            %parfor (z = 1:max(PDsy.Group), BaseData.Parallel(g)*100)
            for z = 1:max(PDsy.Group)
                
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
                if length(AllTrAxGr) < 450/BaseData.ILRes(g)
                    continue
                end % modified temp by lucas true value 2000
                
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
                        
                        % Now the issue is that I have to figure out how
                        % much of MaxLE is from the ST... and remove
                        
                        % Subject Influence Line to Truck Axle Stream
                        [MaxLE,DLF,BrStInd,R] = VBGetMaxLE(AllTrAxSub,ILData(t).v(:,1:length(Lanes)),BaseData.RunDyn(g));
                        
                        if BaseData.RunFat(g) == 1
                            [c,hist,edges,rmm,idx] = rainflow(R);
                            T = array2table(c,'VariableNames',{'Count','Range','Mean','Start','End'});
                        end
                        % Next find a way to save T for every week
                        
                        if MaxLE == 0, k = k+1; continue, end
                        k = k+1; % Add to k
                        
                        % Adjust BrStInd for Starti [now TrLineUpSub(1,1)]
                        BrStIndx = BrStInd + TrLineUpSub(1,1) -1;
                        % Get BrEndIndx
                        BrEndIndx = BrStIndx + BrLengthInd - 1;
                        % Get Bridge Indices
                        BrIndsx = [BrStIndx:BrEndIndx]';
                        BrInds = [BrStInd:BrStInd+BrLengthInd - 1]';
                        %AxOnBr = sum(AllTrAxSub(BrInds,:),2);
                        %sum(sum(AllTrAxSub(BrInds,:).*flip(ILData(t).v(:,1:length(Lanes)))))
                        
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
                            
                            % Call VBWIMtoAllT.rAx w/ mods... must give stationary point or truck
                            [PDe, AllTrAxStop, TrLineUpStop] = VBWIMtoAllTrAxStop(PDe,MaxLength,Lane,BaseData.ILRes(g),find(TrNumsUE == TrIdMax));
                            
                            % Round TrLineUp first row, move unrounded to fifth row
                            TrLineUpStop(:,5) = TrLineUpStop(:,1); TrLineUpStop(:,1) = round(TrLineUpStop(:,1)/BaseData.ILRes(g));
                            % TrLineUpStop [ 1: AllTrAxIndex  2: AxleValue  3: Truck#  4: LaneID  5: Station(m) ]
                            
                            [MaxLEe,DLFe,BrStInde,Re] = VBGetMaxLE(AllTrAxStop,ILData(t).v,BaseData.RunDyn(g));
                        end
                        
                        % We use TrNums because they don't depend on Starti shift
                        MaxLETime = PDsy.DTS(TrNums(1));
                        Vehs = PDsy.CLASS(TrNumsU);
                        
                        %if sum(Vehs == 41) == 0
                        if sum(sum(Vehs == [23])) == 0
                            if sum(PDsy.Group == z & PDsy.CLASS == 23) >= 1
                            if sum(PDsy.Group == z & PDsy.LANE == 2) >= sum(PDsy.Group == z & PDsy.LANE == 1)
                            lucas = 1;
                            end
                            end
                            AllTrAxSub(TrLineUpOnBr(TrLineUpOnBr(:,3) == TrNumsU(1,1),1)-(TrLineUpSub(1,1)-1),mean(TrLineUpOnBr(TrLineUpOnBr(:,3) == TrNumsU(1,1),4))) = 0;
                            k = k + 1;
                            continue
                        end
                        % Get Modified MaxLE. first set alxes from 60t to 0
                        % Find indices of AllTrAxSub... must used
                        % TrLineUpOnBr cols 1 and 4 (but w/ x vs no x
                        % offset? hmm
                        B4 = sum(sum(AllTrAxSub(BrInds,:).*flip(ILData(t).v(:,1:length(Lanes)))));
                        % Set Vehs == 41 axles to 0
                        MaxLEContr = [];
                                                                                                                   
                        for v = 1:length(TrNumsU)
                            AllTrAxSubTemp = AllTrAxSub;
                            AllTrAxSubTemp(TrLineUpOnBr(logical(sum(TrLineUpOnBr(:,3) == (TrNumsU(v))',2)),1)-(TrLineUpSub(1,1)-1),sum((Lanes == (TrLineUpOnBr(logical(sum((TrLineUpOnBr(:,3) == TrNumsU(v)'),2)),4))').*[1:length(Lanes)]',1)') = 0;
                            %Calculate MaxLE again
                            MaxLEContr(v) = sum(sum(AllTrAxSubTemp(BrInds,:).*flip(ILData(t).v(:,1:length(Lanes)))));
                        end
                                                                     
                        %v = find(Vehs == 41);
                        v = find(sum(Vehs == [23],2));
                        %if max(B4-MaxLEContr(v)) == max(B4-MaxLEContr)
                        if max(B4-MaxLEContr(v)) == max(B4-MaxLEContr)
                           [~,vposi] = max(B4-MaxLEContr(v)); % NEW for more than 41
                           v = v(vposi); % NEW for more than 41
                           AllTrAxSub(TrLineUpOnBr(logical(sum(TrLineUpOnBr(:,3) == (TrNumsU(v))',2)),1)-(TrLineUpSub(1,1)-1),sum((Lanes == (TrLineUpOnBr(logical(sum((TrLineUpOnBr(:,3) == TrNumsU(v)'),2)),4))').*[1:length(Lanes)]',1)') = 0;
                           MaxLE = sum(sum(AllTrAxSub(BrInds,:).*flip(ILData(t).v(:,1:length(Lanes)))));
                        else
                            continue
                        end
                        
                        %AllTrAxSub(TrLineUpOnBr(logical(sum(TrLineUpOnBr(:,3) == (TrNumsU(Vehs == 41,1))',2)),1)-(TrLineUpSub(1,1)-1),sum((Lanes == (TrLineUpOnBr(logical(sum((TrLineUpOnBr(:,3) == TrNumsU(Vehs == 41,1)'),2)),4))').*[1:height(Lanes)]',1)') = 0;
                        % Calculate MaxLE again
                        %MaxLE = sum(sum(AllTrAxSub(BrInds,:).*flip(ILData(t).v(:,1:length(Lanes)))));
                        % First step is to find out how much of the time we
                        % actually have an accompanying vehicle...
%                         if length(TrNumsU) > 1
%                             k = 100;
%                         end
                        if MaxLE < 0
                        MaxLE = 0;
                        end
                                              
                        % Get ClassT (in m form for now)
                        if min(Vehs) == 0
                            m = 1;
                        elseif sum(sum(Vehs == [23],2)) > 0
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
                            TrLineUpSubTemp(find(TrLineUpSub(:,3) == PosiRemTruck),:) = [];
                            T = VBApercu(PDsy,'',ILData(t),BrStIndx,TrLineUpSubTemp,MaxLE/ESIA.Total(t),1,Lane,BaseData.ILRes(g));
                            VBApercu(PDsy,'',ILData(t),BrStIndx,TrLineUpSub,B4/ESIA.Total(t),1,Lane,BaseData.ILRes(g));
                            end
                            %exportgraphics(gcf,"Max"  + ".jpg",'Resolution',600)
                            if BaseData.StopSim(g)
                                TStop = VBApercu(PDe,'',ILData(t),BrStInde,TrLineUpStop,MaxLEe/ESIA.Total(t),1,Lane,BaseData.ILRes(g));
                            end
                        end
                        

                        % Save MaxEvents... save Times and Datenums and then convert
                       
                            MaxEvents1 = [MaxEvents1; datenum(MaxLETime), BaseData.SITE(g), MaxLE, t, m, k, BrStIndx];
                        
                        
                        % Rewrite line if DetailedVBWIM Desired
                        if BaseData.StopSim(g)
                            MaxEvents1Stop = [MaxEvents1Stop; datenum(MaxLETime), BaseData.SITE(g), MaxLEe, t, m, k, BrStInde];
                        end
                        
%                         if m == 2
%                             k = 100; % Bump k up so that analysis doesn't continue!
%                         end
                        
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
            LenPrint = VBUpProgBar(st,g,length(UYears),RamUsed(end),r,LenPrint);
                        
        end % r, years
    end % w, SiteGroups
  
    % If there was no determinent vehicules, initialise
    if height(MaxEvents) == 0
    MaxEvents = [datenum(MaxLETime), BaseData.SITE(g), MaxLE, t, m, k, BrStIndx];
    end
    
    % Convert back to datetime
    MaxEvents = array2table(MaxEvents,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'BrStInd'});
    
    MaxEvents.Datenum = datetime(MaxEvents.Datenum,'ConvertFrom','datenum');
    MaxEvents = renamevars(MaxEvents,"Datenum","DTS");
    
    if BaseData.StopSim(g)
        MaxEventsStop = array2table(MaxEventsStop,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'BrStInd'});
        MaxEventsStop.Datenum = datetime(MaxEventsStop.Datenum,'ConvertFrom','datenum');
        MaxEventsStop = renamevars(MaxEventsStop,"Datenum","DTS");
    end
    
    % m = 1 is ClassT 'All', m = 2 is 'ClassOW', and m = 3 is 'Class'
    % qInvestInitial Inputs
    BM = {'Monthly', 'Yearly'};             % j
    %BM = {'Monthly'};
    ClassType = {'ClassOW'};        % i
    DistTypes = {'Lognormal'};
    [Max,~,~,~] = qInvestInitial_60t(BM,ClassType,DistTypes,MaxEvents,ILData); % TROUBLESHOOT
    
    TName = datestr(now,'mmmdd-yy HHMMSS');
    % Need to go back to original BaseData... no SITE switch
    BaseData = VBReadInputFile(FName);
    OutInfo.Name = TName; OutInfo.BaseData = BaseData(g,:);
    % ATTENTION!
    OutInfo.ESIA = ESIA; %OutInfo.E41 = E41;
    OutInfo.ILData = ILData;
    OutInfo.SimStop = false;
    % OutInfo.Max = Max; Need to be done after
    
    % Plot BlockMax, find Design Values, Ed, using Beta, rather than 99th percentile
    for r = 1:Num.InfCases
        for i = 1:length(ClassType)
            for j = 1:length(BM)
            Class = ClassType{i};
            BlockM = BM{j};
            % Proportion of time when there is another truck on the bridge
            PropTrucks = 1-height(Max(r).(Class).(BlockM)(Max(r).(Class).(BlockM).Max == 0,:))/height(Max(r).(Class).(BlockM));
            OutInfo.PropTrucks.(Class).(BlockM)(r) = PropTrucks;
            Max(r).(Class).(BlockM) = Max(r).(Class).(BlockM)(Max(r).(Class).(BlockM).Max > 0,:);
            if height(Max(r).(Class).(BlockM))<= 30
                if r==1
                    OutInfo.x_values.(Class).(BlockM)(:,r) = zeros(100,1);
                    OutInfo.y_valuespdf.(Class).(BlockM)(:,r) = zeros(100,1);
                else
            OutInfo.x_values.(Class).(BlockM)(:,r) = 0;
            OutInfo.y_valuespdf.(Class).(BlockM)(:,r) = 0;
                end
                %{
            if r==1
                OutInfo.x_values.(Class).(BlockM)(:,r) = zeros(100,1);
                OutInfo.y_valuespdf.(Class).(BlockM)(:,r) = zeros(100,1);
            else
               OutInfo.x_values.(Class).(BlockM)(:,r) = 0;
               OutInfo.y_valuespdf.(Class).(BlockM)(:,r) = 0; 
            end
                %}
            OutInfo.EdLN.(Class).(BlockM)(r) = 0;
            OutInfo.AQ.(Class).(BlockM)(r) = 0;
            else
            [~,OutInfo.x_values.(Class).(BlockM)(:,r),OutInfo.y_valuespdf.(Class).(BlockM)(:,r),~] = GetBlockMaxFit_60t(Max(r).(Class).(BlockM).Max,'Lognormal',BaseData.Plots(g));
            %[~,OutInfo.x_values.(Class).(BlockM)(:,r),OutInfo.y_valuespdf.(Class).(BlockM)(:,r),~] = GetBlockMaxFit(Max(r).(Class).(BlockM).Max,'Lognormal',BaseData.Plots(g));
            %[ECDF,ECDFRank,PPx,PPy,Fity,OutInfo.LNFitR2] = GetLogNormPPP(Max(r).(Class).(BlockM).Max,false);
            [OutInfo.EdLN.(Class).(BlockM)(r), OutInfo.AQ.(Class).(BlockM)(r), ~] = GetBlockMaxEd_60t(Max(r).(Class).(BlockM).Max,BlockM,'Lognormal',ESIA.Total(r),ESIA.EQ(:,r),ESIA.Eq(:,r),0.6,0.5,PropTrucks);
            end
            end
         end
    end
    
    OutInfo.Max = Max;
    
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
                [OutInfo.EdLN.(Class).(BlockM)(r), OutInfo.AQ.(Class).(BlockM)(r), ~] = GetBlockMaxEd(Max(r).(Class).(BlockM).Max,BlockM,'Lognormal',ESIA.Total(r),ESIA.EQ(:,r),ESIA.Eq(:,r),0.6,0.5,PropTrucks);
            end
        end
        
        if BaseData.Save(g) == 1
            save(['Output' BaseData.Folder{g} '/' OutInfo.Name], 'OutInfo')
        end
    end
    
end % g, BaseData

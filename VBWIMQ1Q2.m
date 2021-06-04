% ------------------------------------------------------------------------
%                            VBWIMQ1Q2
% ------------------------------------------------------------------------
% Generate a summary of maximum effects on a strip length
% Purpose is to study Q1+Q2, in that way this script is a sister to
% AxleStatsBasic, which focusses on Q1 on its own. Both have live scripts
% which load vars and perform analyses.

% AxleStatsBasic >> AxTandem      >> Q1   Investigation
% VBWIMQ1Q2      >> Q1Q2MaxEvents >> Q1Q2 Investigation
% VBWIMq         >> qMaxEvents    >> q    Investigation

% Initializing commands
clear, clc, tic, format long g, rng('shuffle'), close all;

% Read Input File
BaseData = VBReadInputFile('Input/VBWIMQ1Q2Input.xlsx');

% Initialize variable
MaxEvents = nan(500000,13);

% Each row of BaseData represents one analysis
for g = 1:height(BaseData)

    MaxEvents1 = nan(500000,13);
    j = 1;
    
    % Update analysis data for current row of BaseData
    [Num,Lane,ILData,~,~,ESIA] = VBUpdateData(BaseData(g,:));
    
    % Fix ILData rounding issue for Axle Influence Line
    ILData.v < 1 == 0;
    
    % Load File
    load(['WIM/',num2str(BaseData.SITE(g)),'.mat']);
    
    % Delete any lanes > 2 from the WIM...this way we can do all except 456
    PDs = PDs(PDs.LANE < 3,:);
    
    %Stage2Prune
    UnderW = PDs.GW_TOT<6000;
    TotUnderW = sum(UnderW);
    PDs(UnderW,:) = [];
    
    % 2) Disqualification by Swiss10 Class (exclude 2,3,4,6)
    
    WrongC = PDs.CS == 2 | PDs.CS == 3 | PDs.CS == 4 | PDs.CS == 6;
    TotWrongC = sum(WrongC);
    % Only do the disqualification if we actually have SW10 Classification
    if sum(WrongC) < 0.7*height(PDs)
        PDs(WrongC,:) = [];
    end

    % Separate for each year...
    UYears = unique(year(PDs.DTS));
    
    % Start Progress Bar
    u = StartProgBar(length(UYears), 1, g, height(BaseData)); tic; st = now;
    
    parfor r = 1:length(UYears)
    %for r = 1:length(UYears)
        
        % Take just current year
        PDsy = PDs(year(PDs.DTS) == UYears(r),:);
    
        % Convert PDC to AllTrAx - Spacesave at MaxLength
        [PDsy, AllTrAx, TrLineUp] = VBWIMtoAllTrAx(PDsy,4,Lane,BaseData.ILRes(g));
        
        % Make groups out of each unique day
        PDsy.Group = findgroups(dateshift(PDsy.DTS,'start','day'));
        
        % Round TrLineUp first row, move unrounded to fifth row
        TrLineUp(:,5) = TrLineUp(:,1); TrLineUp(:,1) = round(TrLineUp(:,1)/BaseData.ILRes(g));
        % Expand TrLineUp to include groups
        TrLineUp(:,6) = PDsy.Group(TrLineUp(:,3));
        
        % TrLineUp [     1            2         3         4          5     ]
        %           AllTrAxIndex  AxleValue   Truck#    LaneID   Station(m)
        
        % Perform search for maximums for each day
        for z = 1:max(PDsy.Group)
            
            % Store starting and end indices
            Starti = max(1,min(TrLineUp(TrLineUp(:,6) == z,1))-30);
            Endi = min(max(TrLineUp(TrLineUp(:,6) == z,1)+30),length(AllTrAx));
            
            % Subdivide AllTrAx
            AllTrAxSub = AllTrAx(Starti:Endi,:);
            
            % Don't bother running if the segment is too small
            if length(AllTrAxSub) < 20000
                continue
            end
            k = 0;
            
            % Get length of bridge in number of indices
            BrLengthInd = size(ILData.v,1);
            
            % Eliminate the need for padding or BrStInd index issues
            AllTrAxSub(1:BrLengthInd,:) = 0;
            AllTrAxSub(end-BrLengthInd:end,:) = 0;
            
            % For each analysis
            while k < BaseData.NumAnalyses(g) & sum(AllTrAxSub,'all') > 0
                
                % Subject Influence Line to Truck Axle Stream
                [MaxLE,DLF,BrStInd,R] = VBGetMaxLE(AllTrAxSub,ILData.v,BaseData.RunDyn(g));
                
                % Now add to k
                k = k+1;
                
                % Adjust BrStInd for Starti
                BrStIndx = BrStInd + Starti -1;
                
                % Get BrEndIndx
                BrEndIndx = BrStIndx + BrLengthInd - 1;
                
                % Get Strip Inds
                StripInds = [BrStIndx:BrEndIndx]';
                % Really only take the strip
                StripInds = StripInds(flip(ILData.v(:,1)) == 1);
                %AxOnBr = sum(AllTrAxt(StripInds,:),2);
                
                % Get Key Info to Save
                TrNums = TrLineUp(TrLineUp(:,1) >= min(StripInds) & TrLineUp(:,1) <= max(StripInds),3);
                TrNumsU = unique(TrNums);
                
                % We use TrNums because they don't depend on Starti shift
                MaxLETime = PDsy.DTS(TrNums(1));
                Vehs = PDsy.CLASS(TrNumsU);
                Spds = PDsy.SPEED(TrNumsU);
                Lnes = PDsy.LANE(TrNumsU);
                L1Veh = Vehs(Lnes == 1);
                L2Veh = Vehs(Lnes == 2);
                L1Spd = Spds(Lnes == 1);
                L2Spd = Spds(Lnes == 2);
                % 99 is coded as empty (0 is taken by unclassified)
                if isempty(L1Veh); L1Veh = 99; end
                if isempty(L2Veh); L2Veh = 99; end
                if length(L1Veh) > 1; L1Veh(2) = []; end
                if length(L2Veh) > 1; L2Veh(2) = []; end
                if isempty(L1Spd); L1Spd = -1; end
                if isempty(L2Spd); L2Spd = -1; end
                if length(L1Spd) > 1; L1Spd(2) = []; end
                if length(L2Spd) > 1; L2Spd(2) = []; end
                
                % When retreiving from AllTrAxt we must shift back
                L1Load = sum(AllTrAxSub(StripInds-(Starti-1),1));
                L2Load = sum(AllTrAxSub(StripInds-(Starti-1),2));
                L1Ax = sum(AllTrAxSub(StripInds-(Starti-1),1)>0);
                L2Ax = sum(AllTrAxSub(StripInds-(Starti-1),2)>0);
                
                % Get ClassT (in m form for now)
                if min(Vehs) == 0
                    m = 1;
                elseif sum(Vehs > 39 & Vehs < 50) > 0
                    m = 2;
                else
                    m = 3;
                end
                
                % Optional Apercu
                if BaseData.Apercu(g) == 1
                    ApercuTitle = Lane.Sites.SName + " " + num2str(BaseData.SITE(g)) + " Max " + num2str(k);
                    T = VBApercu(PDsy,ApercuTitle,ILData,BrStIndx,TrLineUp,MaxLE/ESIA.Total,DLF,Lane,BaseData.ILRes(g));
                    %exportgraphics(gcf,"Apercu" + BaseData.Folder(g) + "/" + ApercuTitle + ".jpg",'Resolution',600)
                end
                
                % Save MaxEvents...
                % Save Times and Datenums and then convert
                MaxEvents1(j,:) = [datenum(MaxLETime), BaseData.SITE(g), MaxLE, m, k, L1Veh, L2Veh, L1Load, L2Load, L1Ax, L2Ax, L1Spd, L2Spd];
                j = j+1;
                
                % Prepare for next run - Set Axles to zero in AllTrAx (can't delete because indices are locations)
                AllTrAxSub(StripInds-(Starti-1),:) = 0;
                
            end % k, analyses
        end % z, groups

        % Update progress bar
        UpProgBar(u, st, g, 1, length(UYears), 1)
        
    end % r, years
    MaxEvents(:,:,g) = MaxEvents1;
end % g, BaseData

% Trim back
MaxEvents = reshape(permute(MaxEvents,[1 3 2]),[],13);
MaxEvents(isnan(MaxEvents(:,1)),:) = [];

% Convert back to datetime !!
% Delete empty rows and convert to table
MaxEvents = array2table(MaxEvents,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'm', 'DayRank', 'L1Veh', 'L2Veh', 'L1Load', 'L2Load', 'L1Ax', 'L2Ax', 'L1Sp', 'L2Sp'});
MaxEvents.DTS = datetime(MaxEvents.Datenum,'ConvertFrom','datenum');
MaxEvents.Datenum = [];

% Add in the description of MaxEvents that it is for 1.4 m (justified by previous memos)
MaxEvents.Properties.Description = '1.4 m strip length, 0.1 m ILRes';

% Add Column for All, Class, ClassOW and delete former m
MaxEvents.ClassT(MaxEvents.m == 1) = "All";
MaxEvents.ClassT(MaxEvents.m == 2) = "ClassOW";
MaxEvents.ClassT(MaxEvents.m == 3) = "Class";
MaxEvents.m = [];

% Manual Save

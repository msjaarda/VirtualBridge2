% ------------------------------------------------------------------------
%                            VBWIM
% ------------------------------------------------------------------------
% Find maximum load effects over a bridge from real WIM traffic

% Initializing commands
clear, clc, close all, format long g, rng('shuffle');

% Read Input File
BaseData = VBReadInputFile('Input/VBWIMInput.xlsx');

% Initialize parpool if necessary and initialize progress bar
if BaseData.Parallel(1) > 0, gcp; clc; end
u = StartProgBar(height(BaseData), 1, 1, 1); tic; st = now;

% After that we need to add back specific time period... daily,
% weekly, yearly, with logging of block maxima.
%For Q1/Q2 and q Investigations (but they will have their own names).

% Initialize variables and start row counter
MaxEvents = nan(500000,13); j = 1;

% Each row of BaseData represents one analysis
%parfor g = 1:height(BaseData)
for g = 1:height(BaseData)
    
    % Update analysis data for current row of BaseData
    [Num,Lane,ILData,~,~,ESIA] = VBUpdateData(BaseData(g,:));
    
    % Load File
    %PDx = 
    load(['WIM/',num2str(BaseData.SITE(g)),'.mat']);
    %PDs = PDx.PDs;
    
    % Get Only the ClassType Specified
    try
        if strcmp(BaseData.ClassType(g),'Class')
            PDs = PDs(PDs.CLASS > 0 & (PDs.CLASS > 50 | PDs.CLASS < 40),:);
        elseif strcmp(BaseData.ClassType(g),'ClassOW')
            PDs = PDs(PDs.CLASS > 0),:);
        end
    catch end
    
    % Let's choose a certain amount for a demo
    try
        PDs = PDs(1:800000,:);
    catch
        PDs = PDs(1:end,:);
    end
    
    % Convert PDC to AllTrAx - Spacesave at MaxLength
    MaxLength = (max(arrayfun(@(x) size(x.v,1),ILData))-1)/BaseData.ILRes(g);
    [PDs, AllTrAx, TrLineUp] = VBWIMtoAllTrAx(PDs,MaxLength,Lane,BaseData.ILRes(g));

    
    % Round TrLineUp first row, move unrounded to fifth row
    TrLineUp(:,5) = TrLineUp(:,1); TrLineUp(:,1) = round(TrLineUp(:,1)/BaseData.ILRes(g));
    % TrLineUp [     1            2         3         4          5     ]
    %           AllTrAxIndex  AxleValue   Truck#    LaneID   Station(m)
    
    % For each influence case
    for t = 1:Num.InfCases
        
        % Reset for each t
        AllTrAxt = AllTrAx;
        TrLineUpt = TrLineUp;
        k = 0;
        
        % For each analysis
        while k < BaseData.NumAnalyses(g) & sum(AllTrAxt,'all') > 0
            
            % Subject Influence Line to Truck Axle Stream
            [MaxLE,DLF,BrStInd,R] = VBGetMaxLE(AllTrAxt,ILData(t).v,BaseData.RunDyn(g));
            
            % Get length of bridge in number of indices
            BrLengthInds = size(ILData(t).v,1);
            
            % Add Padding if necessary
            if BrStInd < 1 | BrStInd + BrLengthInds - 1 > height(AllTrAxt)
                % Add Padding
                PadLen = BrLengthInds -1;
                AllTrAxt = [zeros(PadLen,size(AllTrAxt,2)); AllTrAxt; zeros(PadLen,size(AllTrAxt,2))];
                BrStInd = BrStInd + PadLen;
                % Also need to modify TrLineUp
                TrLineUpt(:,1) = TrLineUpt(:,1) + PadLen; TrLineUpt(:,5) = TrLineUpt(:,5) + BaseData.ILRes(g)*PadLen;
            end
            
            BrEndInds = BrStInd + BrLengthInds-1;
            BrInds = [BrStInd:BrEndInds]';
            AxOnBr = sum(AllTrAxt(BrInds,:),2);
            
            % Now add to k since continue has passed
            k = k+1;

            % Optional Apercu
            if BaseData.Apercu(g) == 1
                ApercuTitle = Lane.Sites.SName + " " + num2str(BaseData.SITE(g)) + " Max " + num2str(k);
                T = VBApercu(PDs,ApercuTitle,ILData(t),BrStInd,TrLineUpt,MaxLE/ESIA.Total(t),DLF,Lane,BaseData.ILRes(g));
                exportgraphics(gcf,"Apercu" + "/" + ApercuTitle + ".jpg",'Resolution',600)
            end
            
            % Delete vehicle entries from TrLineUp for re-analysis
            TrLineUpt(TrLineUpt(:,1) > BrStInd & TrLineUpt(:,1) < BrStInd + BrLengthInds,:) = [];
            % Prepare for next run - Set Axles to zero in AllTrAx (can't delete because indices are locations)
            AllTrAxt(BrInds,:) = 0;
            
        end
    end
    
    % Update progress bar
    UpProgBar(u, st, g, 1, height(BaseData), 1)
    
end

function [PDCx, AllTrAx, TrLineUp] = VBWIMtoAllTrAxStop(PDCx,SpaceSaver,Lane,ILRes,TrIndMax)
% WIMTOALLTRAX Translates WIM or VWIM data into AllTrAx and TrLineUp
% Also returns PDC (in the form of PDCx) with some mods

% Detect type of PDCx... basically, does it have time (DTS) or spacing (SpCu)
isWIM =  ismember('DTS', PDCx.Properties.VariableNames);

% Get Lanes
Lanes = unique(PDCx.LANE);

% If WIM, Convert time to distance
if isWIM
    
    OldSpCu = PDCx.SpCu;
    OldTrIndMaxSpCu = OldSpCu(TrIndMax);

    % Just to initialize both vars
    PDCx.Dist = PDCx.LENTH/100+3;
    PDCx.SpCu = cumsum(PDCx.Dist);
    
    % Calculate SpCu of TrIndMax
    TrIndMaxSpCu = PDCx.SpCu(TrIndMax);

    for i = 1:length(Lanes)
                
        % Find indices of the lane we are working in
        LaneInds = PDCx.LANE == Lanes(i);
        
        % Find closest vehicle to TrIndMax (can be TrIndMax itself) use SpCu
        [MinV, MinI] = min(abs(OldSpCu(LaneInds) - OldTrIndMaxSpCu));
        % Get actual mins (could be neg)
        %Act = PDCx.SpCu2(LaneInds) - TrIndMaxSpCu2;
        WOShift = cumsum(PDCx.Dist(LaneInds));
        % Get Lane Shift
        LShift = 1000-WOShift(MinI);
        
        % Cummulative distance in axle stream
        PDCx.SpCu(LaneInds) = cumsum(PDCx.Dist(LaneInds));
        %Temp = PDCx.SpCu(LaneInds);
        % Shift using LShift
        PDCx.SpCu(LaneInds) = PDCx.SpCu(LaneInds) + LShift;

    end
    
    % Cummulative distance in axle stream
    %PDCx.SpCu = cumsum(PDCx.Dist);
    
end

% Spacing is front of veh to front of veh
PDCx.LnTrSpacing = zeros(height(PDCx),1);
% Btw, or between, is rear of veh to front of next
PDCx.LnTrBtw = zeros(height(PDCx),1);

% Some kind of filter for making sure trucks don't encroach on one another
if isWIM
    for i = 1:length(Lanes)
        
        % Find indices of the lane we are working in
        LaneInds = PDCx.LANE == Lanes(i);
        
        % Find all locations where truck i and i - 1 arrived at the same time
        % 30 chosen as greater than max tr length (26) to avoid deleting r1
        AA = [30; diff(PDCx.SpCu(LaneInds))];
        
        PDCx.LnTrSpacing(LaneInds) = AA;
        % The following only makes sense in direction 1. We don't circshift
        % for the 2 direction... why not?
        if Lane.Details.Dir(Lanes(i)) == 1
            PDCx.LnTrBtw(LaneInds) = AA - PDCx.LENTH(circshift(find(LaneInds == 1),1))/100;
        else
            PDCx.LnTrBtw(LaneInds) = AA - PDCx.LENTH(LaneInds)/100;
        end

    end
    
    % If LnTrBtw is negative we delete entry
    PDCx(PDCx.LnTrBtw < -2,:) = [];
end

% Create wheelbase and axle load vectors
WBL = PDCx{:,strncmp(PDCx.Properties.VariableNames,'W',1)}/100;
AX = PDCx{:,strncmp(PDCx.Properties.VariableNames,'AW',2)}/102;

% Make wheelbase length cummulative
WBL = cumsum(WBL,2);

% Switch WBL for direction 2
for i = 1:length(Lanes)
    
    % Find indices of the lane we are working in
    LaneInds = PDCx.LANE == Lanes(i);
    
    % Change the sign of the WBL for those in direction 2
    if Lane.Details.Dir(Lanes(i)) == 2
        WBL(LaneInds,:) = -WBL(LaneInds,:);
    end 
end

WB = [PDCx.SpCu PDCx.SpCu + WBL];

% Must eliminate useless WB values
WB(AX == 0) = 0;
T = ones(size(AX)).*(AX > 0);
TrNum = [1:size(WB,1)]';
Q = repmat(TrNum,1,size(T,2));
TrNum = Q.*T;

% What is going on here with LaneNum?!
LaneNum = PDCx.LANE;
Q = repmat(LaneNum,1,size(T,2));
LaneNum = Q.*T;

x = WB'; WBv = x(:);
x = AX'; AXv = x(:);
x = TrNum'; TrNum = x(:);
x = LaneNum'; LaneNum = x(:);

% v stands for vector (not matrix)
WBv = WBv(WBv > 0);
AXv = AXv(AXv > 0);
TrNum = TrNum(TrNum > 0);
LaneNum = LaneNum(LaneNum > 0);

% Update the below
%AllLaneLineUp = [SpCu(1) AllAxLoads(2) AllVehNum(3) AllLaneNum(4)...
TrLineUp = [WBv AXv TrNum LaneNum];

% The way that the indexing and accumarray is working, we have wasted stuff
% at the start of the AllTrAx... and it is much too long (when using VWIM)
TrLineUp(:,1) = round(TrLineUp(:,1)/ILRes);

% Make a separate axle stream vector for each lane, and last one for all
% Put max() function in incase one lane has no representation in TrLineUp
AllTrAx = zeros(max(TrLineUp(:,1)),max(length(Lane.Details.Dir),length(Lanes)));

for i = 1:length(Lanes)
    A = accumarray(TrLineUp(TrLineUp(:,4)==Lanes(i),1),TrLineUp(TrLineUp(:,4)==Lanes(i),2));
    % Fixed issues 26/8/21 when PDCx only had lane 2 included it would place stuff
    % in lane 1 (for exmaple)
    %AllTrAx(1:length(A(1:1:end)),i) = A(1:1:end); 
    AllTrAx(1:length(A(1:1:end)),Lanes(i)) = A(1:1:end); 
end

% Return TrLineUp first row unrounded
TrLineUp(:,1) = WBv;

end





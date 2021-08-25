function [PDCx, AllTrAx, TrLineUp] = VBWIMtoAllTrAxStop(PDCx,SpaceSaver,Lane,ILRes,TrIndMax)
% WIMTOALLTRAX Translates WIM or VWIM data into AllTrAx and TrLineUp
% Also returns PDC (in the form of PDCx) with some mods

% Goal is to shrink distance between vehs... then add new time stamps.
% Right now we start with deterministic values
SpBtwnAll = 3;
CarLen = 5;
% TrNumCen is the truck which doesn't change time stamp... need one in each
% lane..
TimeDiff = PDCx.DTS-PDCx.DTS(11);
[MTD, MTDL] = min(abs(TimeDiff(PDCx.LANE ~= PDCx.LANE(11))));
% MTDL is the Minimum Time Difference Location... thus TrNumCen for off
% lane
% Summarize which truck doesn't move in each lane...
%PDCx.LANE(11) goes with index 11... note 11 is from 10 less and 10 more hard coded right now
%PDCx.LANE(MTDL) goes with index MTDL...
% Now we need to know the dominant lane, so we can put most cars in it...

% Detect type of PDCx... basically, does it have time (DTS) or spacing (SpCu)
isWIM =  ismember('DTS', PDCx.Properties.VariableNames);
% Else, assume it is VWIM, Apercu, or Det

% Change lane directions so that N and E go to the left (2), and SW right (1)

% Get Lanes
Lanes = unique(PDCx.LANE);

% If WIM, Convert time to distance
if isWIM
    
    % Sort? This adds time... let's presort.
    PDCx = sortrows(PDCx,2);
    
    % Convert time to distance
    PDCx.Dist = [1; seconds(diff(PDCx.DTS))].*((PDCx.SPEED)*0.2777777777778);
    
    % MAIN ISSUE HERE IS THAT IT IS NOT LANE SPECIFIC!! NEED TO USE A LOOP
    % TO ENTER THE LANES AND DO IT FOR EACH LANE...
    
    % Good to do this first as it anchors the lane 2 close...
    % maybe better is to do them separate, then all in the fast lane to be
    % centered around heaviest from fast. Could get main offender as input
    % to the function... actually don't think you can get around this. Must
    % get main offender somehow...
    
    % Just to initialize both vars
     PDCx.Dist = PDCx.LENTH/100+3;
     PDCx.SpCu = cumsum(PDCx.Dist);
     
     % Calculate SpCu of TrIndMax
     TrIndMaxSpCu = PDCx.SpCu(TrIndMax);
    
    for i = 1:length(Lanes)
                
        % Find indices of the lane we are working in
        LaneInds = PDCx.LANE == Lanes(i);
        
        % Do PDCx.Dist lane specifically...
        %PD.Dist(LaneInds)
        
        % Find closest vehicle to TrIndMax (can be TrIndMax itself) use SpCu
        [MinV, MinI] = min(abs(PDCx.SpCu(LaneInds) - TrIndMaxSpCu));
        % Get actual mins (could be neg)
        Act = PDCx.SpCu(LaneInds) - TrIndMaxSpCu;
        % Get Lane Shift
        LShift = 1000-Act(MinI);
        
        % Cummulative distance in axle stream
        PDCx.SpCu(LaneInds) = cumsum(PDCx.Dist(LaneInds));
        %Temp = PDCx.SpCu(LaneInds);
        % Shift using LShift
        PDCx.SpCu(LaneInds) = PDCx.SpCu(LaneInds) + LShift;
        
        % Shift things in lane so that vehicle is beside TrIdMax
        
        
        % Find all locations where truck i and i - 1 arrived at the same time
        % 30 chosen as greater than max tr length (26) to avoid deleting r1
        
       % AA = [30; diff(PDCx.Dist(LaneInds))];
        
        %PDCx.LnTrSpacing(LaneInds) = AA;

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
        if Lane.Dir(i) == 1
            PDCx.LnTrBtw(LaneInds) = AA - PDCx.LENTH(circshift(find(LaneInds == 1),1))/100;
        else
            PDCx.LnTrBtw(LaneInds) = AA - PDCx.LENTH(LaneInds)/100;
        end

    end
    
    % If LnTrBtw is negative we delete entry
    PDCx(PDCx.LnTrBtw < 1.5,:) = [];
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
    if Lane.Dir(i) == 2
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
AllTrAx = zeros(max(TrLineUp(:,1)),max(length(Lane.Dir),length(Lanes)));

for i = 1:length(Lanes)
    A = accumarray(TrLineUp(TrLineUp(:,4)==Lanes(i),1),TrLineUp(TrLineUp(:,4)==Lanes(i),2));
    AllTrAx(1:length(A(1:1:end)),i) = A(1:1:end); 
end

% Return TrLineUp first row unrounded
TrLineUp(:,1) = WBv;

end





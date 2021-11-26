function [PDCx, TrLineUp] = VBWIMtoAllTrAx_Exp(PDCx,SpaceSaver,Lane,ILRes)
% WIMTOALLTRAX Translates WIM or VWIM data into AllTrAx and TrLineUp
% Also returns PDC (in the form of PDCx) with some mods

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
    
    % Delete excess space according to IL... add max veh length
    if SpaceSaver > 0
        PDCx.Dist(PDCx.Dist > SpaceSaver + 26) = SpaceSaver + 26;
    end
    
    % Cummulative distance in axle stream
    PDCx.SpCu = cumsum(PDCx.Dist);
    
end

% ATTENTION FOR PLATOONING WIM INTEGRATION...
% WE NEED TO MAINTAIN THE LANE-2-LANE CORRELATION... NEED TO MODIFY
% INTRAPLATOON DIST WHILE NOT TOUCHING THIS...

% CAN WE DO THE CUMSUM LIKE
% NORMAL... then turn around and alter the PDCx.DIST for PDCx.PLAT == true

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
        if Lane.Dir(Lanes(i)) == 1
            PDCx.LnTrBtw(LaneInds) = AA - PDCx.LENTH(circshift(find(LaneInds == 1),1))/100;
        else
            PDCx.LnTrBtw(LaneInds) = AA - PDCx.LENTH(LaneInds)/100;
        end

    end
    
    % If LnTrBtw is negative we delete entry
    PDCx(PDCx.LnTrBtw < 1.5,:) = [];
end

% Create wheelbase and axle load vectors
WB = PDCx{:,strncmp(PDCx.Properties.VariableNames,'W',1)}/100;
AX = PDCx{:,strncmp(PDCx.Properties.VariableNames,'AW',2)}/102;

% Make wheelbase length cummulative
WB = cumsum(WB,2);

% Switch WBL for direction 2
for i = 1:length(Lanes)
    
    % Find indices of the lane we are working in
    LaneInds = PDCx.LANE == Lanes(i);
    
    % Change the sign of the WBL for those in direction 2
    if Lane.Dir(Lanes(i)) == 2
        WB(LaneInds,:) = -WB(LaneInds,:);
    end 
end

WB = [PDCx.SpCu PDCx.SpCu + WB];

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
% AllTrAx = zeros(max(TrLineUp(:,1)),max(length(Lane.Dir),length(Lanes)));
% 
% for i = 1:length(Lanes)
%     A = accumarray(TrLineUp(TrLineUp(:,4)==Lanes(i),1),TrLineUp(TrLineUp(:,4)==Lanes(i),2));
%     AllTrAx(1:length(A),i) = A;
%     %AllTrAx(1:length(A(1:1:end)),i) = A(1:1:end); 
% end

% Return TrLineUp first row unrounded
TrLineUp(:,1) = WBv;

end





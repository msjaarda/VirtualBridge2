function [PDs] = FindDup(PDs,Print,Save)
%FINDDUP Go through and identify duplicate vehicles in side-by-side
%lanes... created when a HV changes lanes and registers in both

% Look for instances of same-time vehicle arrival
TDiff = [1; seconds(diff(PDs.DTS))];

% FullLog means Full Logical, for the whole PDs
FullLogSame = TDiff == 0;
% Ind means the actual indexes... Log will correspond to the same length
IndSame = find(FullLogSame);

% Print Results
if Print
    fprintf('\nNumber of Same-Time Arrivals: \t%i/%i (%.2f %%)\n',sum(FullLogSame),height(PDs),100*sum(FullLogSame)/height(PDs));
end

% It is between IndSame and IndSame - 1
% datestr(PDs.DTS(IndSame(1)),'HH:MM:SS:FFF')
% datestr(PDs.DTS(IndSame(1)-1),'HH:MM:SS:FFF')
% PDs(IndSame(1)-2:IndSame(1)+1,:)

% Create two tables. 
% The first table is the first vehicle arrival PDs, the second is the second.
PD1 = PDs(IndSame-1,:); PD2 = PDs(IndSame,:);

% Now tell how many of these are from the same lane
LogSameLANE = PD1.LANE == PD2.LANE;

% Report
if Print
    fprintf('Number from the same LANE: \t\t%i/%i (%.2f %%)\n\n',sum(LogSameLANE),sum(FullLogSame),100*sum(LogSameLANE)/sum(FullLogSame));
end
% Delete
PD1(LogSameLANE,:) = []; PD2(LogSameLANE,:) = [];
% Must wait until the end to delete this one... changes PDs Indexes
DeleteAtEnd = IndSame(LogSameLANE); IndSame(LogSameLANE) = []; 

% Recreate with PD1 being those in Lane 1 and PD2, 2
PDs1 = sortrows([PD1(PD1.LANE == 1,:); PD2(PD2.LANE == 1,:)]);
PDs2 = sortrows([PD1(PD1.LANE == 2,:); PD2(PD2.LANE == 2,:)]);
clear PD1, clear PD2

% See how many are the same SPEED and AX
LogSameSPEED = PDs1.SPEED == PDs2.SPEED;
% Same AX and SPEED
LogSameSPEEDnAX = PDs1.AX == PDs2.AX & LogSameSPEED;
% Same AX different SPEED
LogSameAXDiffSPEED = PDs1.AX == PDs2.AX & ~LogSameSPEED;

% Report
if Print
    fprintf('Number with the same SPEED: \t%i/%i (%.2f %%)\n',sum(LogSameSPEED),height(PDs1),100*sum(LogSameSPEED)/height(PDs1));
    fprintf('Same AX of same SPEED: \t\t\t%i/%i (%.2f %%)\n',sum(LogSameSPEEDnAX),sum(LogSameSPEED),100*sum(LogSameSPEEDnAX)/sum(LogSameSPEED));
    fprintf('Same AX of not same SPEED: \t\t%i/%i (%.2f %%)\n\n',sum(LogSameAXDiffSPEED),sum(~LogSameSPEED),100*sum(LogSameAXDiffSPEED)/sum(~LogSameSPEED));
end

% Find those with GW_TOT within 7.5%
LogSameGW = abs(PDs1.GW_TOT-PDs2.GW_TOT)./PDs1.GW_TOT < 0.1;
LogSameGWnSPEEDnAX = LogSameGW & LogSameSPEEDnAX;
% Report
if Print
    fprintf('Same GW of same SPEED/AX: \t\t\t%i/%i (%.2f %%)\n',sum(LogSameGWnSPEEDnAX),sum(LogSameSPEEDnAX),100*sum(LogSameGWnSPEEDnAX)/sum(LogSameSPEEDnAX));
end

% Find those with AWT within 10%
% See if all axles within 10%
% Get column names starting with AWT
CIndAxs = contains(PDs1.Properties.VariableNames, 'AWT');
AXCheck = abs((PDs1{:,CIndAxs}-PDs2{:,CIndAxs})./PDs1{:,CIndAxs});
% Replace all NaNs
AXCheck(isnan(AXCheck)) = 0;
LogAXCheck = AXCheck < 0.125;
LogSameAWT = all(LogAXCheck,2);
LogSameAWTnGWnSPEEDnAX = LogSameAWT & LogSameGWnSPEEDnAX;
% Report
if Print
    fprintf('Same AWT and GWT of same SPEED/AX: \t%i/%i (%.2f %%)\n',sum(LogSameAWTnGWnSPEEDnAX),sum(LogSameSPEEDnAX),100*sum(LogSameAWTnGWnSPEEDnAX)/sum(LogSameSPEEDnAX));
end

% Let's do the same for wheelbase!
% See if all wheelbases within 10%
% Get column names for WB
CIndWbs = logical(contains(PDs.Properties.VariableNames, 'W').*contains(PDs.Properties.VariableNames, '_'));
CIndWbs(find(string(PDs.Properties.VariableNames) == 'GW_TOT')) = 0; %#ok<FNDSB>
WBCheck = abs((PDs1{:,CIndWbs}-PDs2{:,CIndWbs})./PDs1{:,CIndWbs});
% Replace all NaNs
WBCheck(isnan(WBCheck)) = 0;
LogWBCheck = WBCheck < 0.1;
LogSameWB = all(LogWBCheck,2);
LogSameWBnAWTnGWnSPEEDnAX = LogSameWB & LogSameAWTnGWnSPEEDnAX;
% Report
if Print
    fprintf('Same WB/AWT/GWT of same SPEED/AX: \t%i/%i (%.2f %%)\n\n',sum(LogSameWBnAWTnGWnSPEEDnAX),sum(LogSameSPEEDnAX),100*sum(LogSameWBnAWTnGWnSPEEDnAX)/sum(LogSameSPEEDnAX));
end

PDs.Dup = ~logical(1:height(PDs))';
PDs.Dup(IndSame(LogSameWBnAWTnGWnSPEEDnAX)) = true;
PDs.Dup(IndSame(LogSameWBnAWTnGWnSPEEDnAX)-1) = true;
% as well as -1!

% Get histogram of weights for offenders
%histogram(PDs.GW_TOT(IndSame(LogSameWBnAWTnGWnSPEEDnAX))/1000);
% Show
[~, b] = max(PDs.GW_TOT(IndSame(LogSameWBnAWTnGWnSPEEDnAX))/1000);
IndB = IndSame(LogSameWBnAWTnGWnSPEEDnAX);
%PDs(IndB(b)-1:IndB(b),:)

PDs(DeleteAtEnd,:) = []; 

if Save
    % Save
    save(strcat('WIMx\',num2str(Name)),'PDs','-v7.3')
end


end


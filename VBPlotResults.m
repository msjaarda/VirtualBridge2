clear, clc, close all

% Get VBResults using VBOutput2Struct
[VBResults] = VBOutput2Struct('NewILLib');

% Now we simply say what we want to plot!!
% Format is:
% VBResults.AQ.Type.SubType.Width.Layout.Support.Trans.AE.Traffic
% VBResults.x.Type.SubType.Width.Layout.Support.Trans.AE.Traffic

% General Info
Type = 'Slab';         SubType = 'Standard';
Width = 'Stand9';     Layout = 'Unidirectional';
Support = 'Fixed'; Trans = 'p1';
% Subplot Info (m 1:3 is the subplots)
AE(1) = "Mn"; AE(2) = "Mp"; AE(3) = "V";
% Series Info (i 1:4 are the series')
Traffic(1) = "Denges2019";   Traffic(2) = "Oberburen2019";
Traffic(3) = "Gotthard2019"; Traffic(4) = "Ceneri2019";

% Other deets
FName = 'Slab Bridges'; FigNum = 0;
Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
D = linspecer(24);      % Set colours

% Prepare plot data for VBTriPlot
for m = 1:3
    % All ClassTypes
    for i = 1:4
        % Get ydata
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        % xdata should always be equal... give a warning if one isn't?
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        % DisplayName
        PDets.DN{m,i} = Traffic(i);
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDets.MEC{m}(i,:) = [0 0 0]; PDets.MFC{m}(i,:) = D(i*6-3,:); 
    end
end

% Plot
FigNum = VBTriPlot(xdata,ydata,PDets,Title,'SIM',FigNum,FName);


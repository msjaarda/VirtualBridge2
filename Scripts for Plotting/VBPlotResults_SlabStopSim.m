%clear, clc, close all

% Get VBResults using VBOutput2Struct
[VBResults] = VBOutput2Struct('Uni2LSlab');
[VBResultsx] = VBOutput2Struct('Uni2LSlabStop');

% Now we simply say what we want to plot!!
% Format is:
% VBResults.AQ.Type.SubType.Width.Layout.Support.Trans.AE.Traffic
% VBResults.x.Type.SubType.Width.Layout.Support.Trans.AE.Traffic

% General Info
Type = 'Slab';         SubType = 'Long';
Width = 'Wid9';     Layout = 'Uni';
Support = 'Semi'; Trans = 'p1';
% Subplot Info (m 1:3 is the subplots)
AE(1) = "Mn"; AE(2) = "Mp"; AE(3) = "V";
% Series Info (i 1:4 are the series')
Traffic(1) = "All";   Traffic(2) = "ClassOW";
Traffic(3) = "Class"; %Traffic(4) = "Ceneri2019";

% Other deets
FName = 'Long Slab Bridge'; FigNum = 0;
Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
D = linspecer(24);      % Set colours

% Prepare plot data for VBTriPlot
for m = 1:3
    % All ClassTypes
    for i = 1:length(Traffic)
        % Get ydata
        ydata{m,i+3} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        ydata{m,i} = VBResultsx.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        % xdata should always be equal... give a warning if one isn't?
        xdata{m,i+3} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        % DisplayName
        PDets.DN{m,i+3} = Traffic(i); PDets.DN{m,i} = 'off';
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i+3) = 5; PDets.MEC{m}(i,:) = [0 0 0]; PDets.MFC{m}(i+3,:) = D(i*6-3,:);
        PDets.MS{m}(i) = 5; PDets.MEC{m}(i+3,:) = [0 0 0]; PDets.MFC{m}(i,:) = [1 1 1]; % [.98 .98 .98] is 'none'
    end
end

% Plot
FigNum = VBTriPlot(xdata,ydata,PDets,Title,'WIM',FigNum,FName);
clear, clc, close all

% Get VBResults using VBOutput2Struct
%[VBResults] = VBOutput2Struct('FindProb02');
[VBResults] = VBOutput2Struct('Slab10mProb');

% Now we simply say what we want to plot!!
% Format is:
% VBResults.AQ.Type.SubType.Width.Layout.Support.Trans.AE.Traffic
% VBResults.x.Type.SubType.Width.Layout.Support.Trans.AE.Traffic

% General Info
Type = 'Slab';         SubType = 'Short';
Width = 'Wid9';     Layout = 'Uni';
Support = 'Semi'; Trans = 'p1';
% Subplot Info (m 1:3 is the subplots)
AE(1) = "Mn"; AE(2) = "Mp"; AE(3) = "V";
% Series Info (i 1:4 are the series')
Traffic(1) = "All";   Traffic(2) = "ClassOW";
Traffic(3) = "Class"; %Traffic(4) = "Ceneri2019";

% Other deets
FName = '0.2 Short Semi Slab Bridges'; FigNum = 0;
Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
D = linspecer(24);      % Set colours

% Prepare plot data for VBTriPlot
for m = 1:3
    % All ClassTypes
    for i = 1:length(Traffic)
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
FigNum = VBTriPlot(xdata,ydata,PDets,Title,'WIM',FigNum,FName);




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
FName = '0.2 Long Semi Slab Bridges'; FigNum = 0;
Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
D = linspecer(24);      % Set colours

% Prepare plot data for VBTriPlot
for m = 1:3
    % All ClassTypes
    for i = 1:length(Traffic)
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
FigNum = VBTriPlot(xdata,ydata,PDets,Title,'WIM',FigNum,FName);



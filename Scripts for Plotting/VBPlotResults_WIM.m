clear, clc, close all

% Global Definitions
FigNum = 0; D = linspecer(30);      % Set colours
Type = 'Box';         SubType = 'Stand';
Support = 'Simp';     Trans = 'p0';
% Subplot Info (m)
AE(1) = "Mn"; AE(2) = "Mp"; AE(3) = "V";
Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
%#ok<*SAGROW> ... to prevent warning messages
% Get VBResults using VBOutput2Struct
[VBResults] = VBOutput2Structx('WIMOct18pr');
%[VBResultsPP] = VBOutput2Struct('WIMClassPP');
%load('VBResultsWIMOct5Redo.mat')
% Class
Class = 'ClassOW'; ClassS = 'StopClassOW';

Export = true;



% BOX --------------------------------------------------------------------

% General Info
Width = 'Wid12';      Layout = 'Uni';
% Series Info (i)
Traffic = fieldnames(VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).Mn);
% Set overall at the end to be top layer
Temp = Traffic{1}; Traffic(1) = []; Traffic{end+1} = Temp;

FName = 'Box Girder Bridge Unidirectional 2L (Class+)';
FNameS = 'Box Girder Bridge Unidirectional 2L (Class+) StopSimR';
for m = 1:3
    for i = 1:length(Traffic)
        % Get ydata and xdata
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(Class);
        ydataS{m,i} = 0.8660*VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(ClassS);
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(Class);
        % DisplayName
        PDets.DN{m,i} = Traffic{i}; PDetsS.DN{m,i} = Traffic{i};
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDetsS.MS{m}(i) = 5; 
        PDets.MEC{m}(i,:) = [0 0 0]; PDetsS.MEC{m}(i,:) = [.5 .5 .5]; 
        PDets.MFC{m}(i,:) = D(i*2,:); PDetsS.MFC{m}(i,:) = D(i*2,:);
        if i == length(Traffic)
            PDets.MFC{m}(i,:) = [0 0 0]; PDetsS.MFC{m}(i,:) = [0 0 0];
        end
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIM',FigNum,FName); 
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end
FigNum = VBTriPlotPro(xdata,ydataS,PDetsS,Title,'WIM',FigNum,FNameS); clear ydata ydataS xdata PDets PDetsS
if Export
    exportgraphics(gcf,FNameS + ".jpg",'Resolution',600);
end


% General Info
Width = 'Wid12';      Layout = 'Bi';
clear Traffic
Traffic = fieldnames(VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).Mn);
% Set overall at the end to be top layer
Temp = Traffic{1}; Traffic(1) = []; Traffic{end+1} = Temp;

FName = 'Box Girder Bridge Bidirectional 2L (All)';
FNameS = 'Box Girder Bridge Bidirectional 2L (All) StopSimR';
for m = 1:3
    for i = 1:length(Traffic)
        % Get ydata and xdata
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(Class);
        ydataS{m,i} = 0.8660*VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(ClassS);
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(Class);
        % DisplayName
        PDets.DN{m,i} = Traffic{i}; PDetsS.DN{m,i} = Traffic{i};
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDetsS.MS{m}(i) = 5; 
        PDets.MEC{m}(i,:) = [0 0 0]; PDetsS.MEC{m}(i,:) = [.5 .5 .5]; 
        PDets.MFC{m}(i,:) = D(i*2,:); PDetsS.MFC{m}(i,:) = D(i*2,:);
        if i == length(Traffic)
            PDets.MFC{m}(i,:) = [0 0 0]; PDetsS.MFC{m}(i,:) = [0 0 0];
        end
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIM',FigNum,FName); 
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end
FigNum = VBTriPlotPro(xdata,ydataS,PDetsS,Title,'WIM',FigNum,FNameS); clear ydata ydataS xdata PDets PDetsS
if Export
    exportgraphics(gcf,FNameS + ".jpg",'Resolution',600);
end


% General Info
Width = 'Wid18';      Layout = 'Bi';
clear Traffic
Traffic = fieldnames(VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).Mn);
% Set overall at the end to be top layer
Temp = Traffic{1}; Traffic(1) = []; Traffic{end+1} = Temp;

FName = 'Box Girder Bridge Bidirectional 4L (All)';
FNameS = 'Box Girder Bridge Bidirectional 4L (All) StopSimR';
for m = 1:3
    for i = 1:length(Traffic)
        % Get ydata and xdata
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(Class);
        ydataS{m,i} = 0.8660*VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(ClassS);
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(Class);
        % DisplayName
        PDets.DN{m,i} = Traffic{i}; PDetsS.DN{m,i} = Traffic{i};
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDetsS.MS{m}(i) = 5; 
        PDets.MEC{m}(i,:) = [0 0 0]; PDetsS.MEC{m}(i,:) = [.5 .5 .5]; 
        PDets.MFC{m}(i,:) = D(i*2,:); PDetsS.MFC{m}(i,:) = D(i*2,:);
        if i == length(Traffic)
            PDets.MFC{m}(i,:) = [0 0 0]; PDetsS.MFC{m}(i,:) = [0 0 0];
        end
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIM',FigNum,FName);
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end
FigNum = VBTriPlotPro(xdata,ydataS,PDetsS,Title,'WIM',FigNum,FNameS); clear ydata ydataS xdata PDets PDetsS
if Export
    exportgraphics(gcf,FNameS + ".jpg",'Resolution',600);
end


% General Info
Width = 'Wid12';      Layout = 'Uni';
% Series Info
clear Traffic
Traffic = fieldnames(VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).Mn);
% Set overall at the end to be top layer
Temp = Traffic{1}; Traffic(1) = []; Traffic{end+1} = Temp;
clear Class
Class{1} = 'Class'; Class{2} = 'ClassOW'; Class{3} = 'All'; 
Class{4} = 'StopClass'; Class{5} = 'StopClassOW'; Class{6} = 'StopAll'; 

FName = 'Box Girder Bridge Unidirectional 2L Filter Comparison';
for m = 1:3
    for i = 1:length(Class)
        % Get ydata and xdata
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{end}).(Class{i});
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{end}).(Class{i});
        % DisplayName
        PDets.DN{m,i} = Class{i};
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDets.MEC{m}(i,:) = [0 0 0]; 
        PDets.MFC{m}(i,:) = D(i*4,:);
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIM',FigNum,FName); clear ydata xdata PDets
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end






% TWIN -------------------------------------------------------------------

% Global Definitions
FigNum = 0; D = linspecer(35);      % Set colours
Type = 'Twin';         SubType = 'Stand';
Support = 'Simp';     Trans = 'p0';
% Subplot Info (m)
AE(1) = "Mn"; AE(2) = "Mp"; AE(3) = "V";
Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
clear Class
% Class
Class = 'ClassOW'; ClassS = 'StopClassOW';

% General Info
Width = 'Wid9';      Layout = 'Uni';
% Series Info (i)
Traffic = fieldnames(VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).Mn);
% Set overall at the end to be top layer
Temp = Traffic{1}; Traffic(1) = []; Traffic{end+1} = Temp;

FName = 'Twin Girder Bridge Unidirectional 2L (Class+)';
FNameS = 'Twin Girder Bridge Unidirectional 2L (Class+) StopSimR';
for m = 1:3
    for i = 1:length(Traffic)
        % Get ydata and xdata
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(Class);
        ydataS{m,i} = 0.8660*VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(ClassS);
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(Class);
        % DisplayName
        PDets.DN{m,i} = Traffic{i}; PDetsS.DN{m,i} = Traffic{i};
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDetsS.MS{m}(i) = 5; 
        PDets.MEC{m}(i,:) = [0 0 0]; PDetsS.MEC{m}(i,:) = [.5 .5 .5]; 
        PDets.MFC{m}(i,:) = D(i*2,:); PDetsS.MFC{m}(i,:) = D(i*2,:);
        if i == length(Traffic)
            PDets.MFC{m}(i,:) = [0 0 0]; PDetsS.MFC{m}(i,:) = [0 0 0];
        end
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIM',FigNum,FName);
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end
FigNum = VBTriPlotPro(xdata,ydataS,PDetsS,Title,'WIM',FigNum,FNameS); clear ydata ydataS xdata PDets PDetsS
if Export
    exportgraphics(gcf,FNameS + ".jpg",'Resolution',600);
end


% General Info
SubType = 'Conc';
% Series Info (i)
Traffic = fieldnames(VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).Mn);
% Set overall at the end to be top layer
Temp = Traffic{1}; Traffic(1) = []; Traffic{end+1} = Temp;

FName = 'Concrete Twin Girder Bridge Unidirectional 2L (Class+)';
FNameS = 'Concrete Twin Girder Bridge Unidirectional 2L (Class+) StopSimR';
for m = 1:3
    for i = 1:length(Traffic)
        % Get ydata and xdata
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(Class);
        ydataS{m,i} = 0.8660*VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(ClassS);
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(Class);
        % DisplayName
        PDets.DN{m,i} = Traffic{i}; PDetsS.DN{m,i} = Traffic{i};
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDetsS.MS{m}(i) = 5; 
        PDets.MEC{m}(i,:) = [0 0 0]; PDetsS.MEC{m}(i,:) = [.5 .5 .5]; 
        PDets.MFC{m}(i,:) = D(i*2,:); PDetsS.MFC{m}(i,:) = D(i*2,:);
        if i == length(Traffic)
            PDets.MFC{m}(i,:) = [0 0 0]; PDetsS.MFC{m}(i,:) = [0 0 0];
        end
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIM',FigNum,FName);
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end
FigNum = VBTriPlotPro(xdata,ydataS,PDetsS,Title,'WIM',FigNum,FNameS); clear ydata ydataS xdata PDets PDetsS
if Export
    exportgraphics(gcf,FNameS + ".jpg",'Resolution',600);
end












% SLAB -------------------------------------------------------------------


% Global Definitions
FigNum = 0; D = linspecer(30);      % Set colours
%Support = 'Semi';     Trans = 'p3';
Support = 'Fixed';     Trans = 'p1';
% Subplot Info (m)
AE(1) = "Mn"; AE(2) = "Mp"; AE(3) = "V";
%AE(1) = "Mn"; AE(2) = "MxMid"; AE(3) = "MxEdg";
Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
%Title{1} = 'M-'; Title{2} = 'MxMid'; Title{3} = 'MxEdg';
clear Class
% Class
Class = 'ClassOW'; ClassS = 'StopClassOW';

Type = 'Slab';         SubType = 'Short';
% General Info
Width = 'Wid18';      Layout = 'Bi';

% Series Info (i)
Traffic = fieldnames(VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).Mn);
% Set overall at the end to be top layer
Temp = Traffic{1}; Traffic(1) = []; Traffic{end+1} = Temp;

FName = 'Slab Bridge Fixed p1 Bidirectional 4L (Class+)';
for m = 1:3
    for i = 1:length(Traffic)
        % Get ydata and xdata
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(Class); 
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(Class);
        % DisplayName
        PDets.DN{m,i} = Traffic{i};
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDets.MEC{m}(i,:) = [0 0 0]; 
        PDets.MFC{m}(i,:) = D(i*2,:);
        % Turn overall (end) black
        if i == length(Traffic)
            PDets.MFC{m}(i,:) = [0 0 0];
        end
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIM',FigNum,FName); clear ydata xdata PDets
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end

Type = 'Slab';         SubType = 'Long';
% General Info
Width = 'Wid18';      Layout = 'Bi';
% Series Info (i)
Traffic = fieldnames(VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).Mn);
% Set overall at the end to be top layer
Temp = Traffic{1}; Traffic(1) = []; Traffic{end+1} = Temp;

FName = 'Slab Bridge Fixed p1 Bidirectional 4L (Class+)';
FNameS = 'Slab Bridge Fixed p1 Bidirectional 4L (Class+) StopSimR';
for m = 1:3
    for i = 1:length(Traffic)
        % Get ydata and xdata
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(Class);
        ydataS{m,i} = 0.8660*VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(ClassS);
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{i}).(Class);
        % DisplayName
        PDets.DN{m,i} = Traffic{i}; PDetsS.DN{m,i} = Traffic{i};
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDetsS.MS{m}(i) = 5; 
        PDets.MEC{m}(i,:) = [0 0 0]; PDetsS.MEC{m}(i,:) = [.5 .5 .5]; 
        PDets.MFC{m}(i,:) = D(i*2,:); PDetsS.MFC{m}(i,:) = D(i*2,:);
        if i == length(Traffic)
            PDets.MFC{m}(i,:) = [0 0 0]; PDetsS.MFC{m}(i,:) = [0 0 0];
        end
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIM',FigNum,FName);
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end
FigNum = VBTriPlotPro(xdata,ydataS,PDetsS,Title,'WIM',FigNum,FNameS); clear ydata ydataS xdata PDets PDetsS
if Export
    exportgraphics(gcf,FNameS + ".jpg",'Resolution',600);
end



% GLOBAL COMPARISON -------------------------------------------------------


% Global Definitions
FigNum = 0; D = linspecer(30);      % Set colours

clear Type SubType Support Trans

Type{1} = 'Box';         SubType{1} = 'Stand';
Support{1} = 'Simp';     Trans{1} = 'p0';

Type{2} = 'Twin';         SubType{2} = 'Stand';
Support{2} = 'Simp';     Trans{2} = 'p0';

Type{3} = 'Twin';         SubType{3} = 'Conc';
Support{3} = 'Simp';     Trans{3} = 'p0';

Type{4} = 'Slab';         SubType{4} = 'Short';
Support{4} = 'Fixed';     Trans{4} = 'p1';

Type{5} = 'Slab';         SubType{5} = 'Long';
Support{5} = 'Fixed';     Trans{5} = 'p1';

% Subplot Info (m)
AE(1) = "Mn"; AE(2) = "Mp"; AE(3) = "V";
Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
%#ok<*SAGROW> ... to prevent warning messages
% Get VBResults using VBOutput2Struct
%[VBResults] = VBOutput2Struct('WIM');
% Class
Class = 'ClassOW'; ClassS = 'StopClassOW';

% General Info
clear Width
Width{1} = 'Wid12';  Width{2} = 'Wid9';  Width{3} = 'Wid9';   Width{4} = 'Wid9';  Width{5} = 'Wid9';    

Layout = 'Uni';
Traffic = 'Uni2L';

FName = 'Unidirectional 2L Comparison (Class+)';
for m = 1:3
    t
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end

FName = 'Unidirectional 2L Comparison (Class+) StopSimR';
for m = 1:3
    for i = 1:length(Type)
        % Get ydata and xdata
        try
            ydata{m,i} = 0.8660*VBResults.AQ.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(ClassS);
            xdata{m,i} = VBResults.x.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(ClassS);
        catch
            ydata{m,i} = VBResults.AQ.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(Class);
            xdata{m,i} = VBResults.x.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(Class);
        end
        % DisplayName
        PDets.DN{m,i} = [Type{i} + " " + SubType{i} + " " +  Support{i} + " " + Trans{i}];
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDets.MEC{m}(i,:) = [0 0 0]; 
        PDets.MFC{m}(i,:) = D(i*4,:);
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIMComp',FigNum,FName); clear ydata xdata PDets
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end




clear Type SubType Support Trans

Type{1} = 'Box';         SubType{1} = 'Stand';
Support{1} = 'Simp';     Trans{1} = 'p0';

Type{2} = 'Slab';         SubType{2} = 'Short';
Support{2} = 'Fixed';     Trans{2} = 'p1';

Type{3} = 'Slab';         SubType{3} = 'Long';
Support{3} = 'Fixed';     Trans{3} = 'p1';

% General Info
Width{1} = 'Wid18';  Width{2} = 'Wid18';  Width{3} = 'Wid18'; 

Layout = 'Bi';
Traffic = 'Bi4L';

FName = 'Bidirectional 4L Comparison (Class+)';
for m = 1:3
    for i = 1:length(Type)
        % Get ydata and xdata
        ydata{m,i} = VBResults.AQ.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(Class);
        xdata{m,i} = VBResults.x.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(Class);
        % DisplayName
        PDets.DN{m,i} = [Type{i} + " " + SubType{i} + " " +  Support{i} + " " + Trans{i}];
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDets.MEC{m}(i,:) = [0 0 0]; 
        PDets.MFC{m}(i,:) = D(i*4,:);
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIMComp',FigNum,FName); clear ydata xdata PDets
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end



FName = 'Bidirectional 4L Comparison (Class+) StopSimR';
for m = 1:3
    for i = 1:length(Type)
        % Get ydata and xdata
        try
            ydata{m,i} = 0.8660*VBResults.AQ.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(ClassS);
            xdata{m,i} = VBResults.x.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(ClassS);
        catch
            ydata{m,i} = VBResults.AQ.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(Class);
            xdata{m,i} = VBResults.x.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(Class);
        end
        % DisplayName
        PDets.DN{m,i} = [Type{i} + " " + SubType{i} + " " +  Support{i} + " " + Trans{i}];
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDets.MEC{m}(i,:) = [0 0 0]; 
        PDets.MFC{m}(i,:) = D(i*4,:);
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIMComp',FigNum,FName); clear ydata xdata PDets
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end






% Global Definitions
FigNum = 0; D = linspecer(30);      % Set colours

clear Type SubType Support Trans Class Width

Type{1} = 'Box';         SubType{1} = 'Stand';
Support{1} = 'Simp';     Trans{1} = 'p0';

Type{2} = 'Twin';         SubType{2} = 'Stand';
Support{2} = 'Simp';     Trans{2} = 'p0';

Type{3} = 'Twin';         SubType{3} = 'Conc';
Support{3} = 'Simp';     Trans{3} = 'p0';

Type{4} = 'Slab';         SubType{4} = 'Short';
Support{4} = 'Fixed';     Trans{4} = 'p1';

Type{5} = 'Slab';         SubType{5} = 'Long';
Support{5} = 'Fixed';     Trans{5} = 'p1';

% Subplot Info (m)
AE(1) = "Mn"; AE(2) = "Mp"; AE(3) = "V";
Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
%#ok<*SAGROW> ... to prevent warning messages
% Get VBResults using VBOutput2Struct
%[VBResults] = VBOutput2Struct('WIM');
% Class
Class = 'ClassOW'; ClassS = 'StopClassOW';

% General Info
Width{1} = 'Wid12';  Width{2} = 'Wid9';  Width{3} = 'Wid9';   Width{4} = 'Wid9';  Width{5} = 'Wid9';    

Layout = 'Bi';
Traffic = 'Bi2L';

FName = 'Bidirectional 2L Comparison (Class+)x';
for m = 1:3
    for i = 1:length(Type)
        % Get ydata and xdata
        ydata{m,i} = VBResults.AQ.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(Class);
        xdata{m,i} = VBResults.x.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(Class);
        % DisplayName
        PDets.DN{m,i} = [Type{i} + " " + SubType{i} + " " +  Support{i} + " " + Trans{i}];
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDets.MEC{m}(i,:) = [0 0 0]; 
        PDets.MFC{m}(i,:) = D(i*4,:);
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIMComp',FigNum,FName); clear ydata xdata PDets
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end

FName = 'Bidirectional 2L Comparison (Class+) StopSimRx';
for m = 1:3
    for i = 1:length(Type)
        % Get ydata and xdata
        try
            ydata{m,i} = 0.8660*VBResults.AQ.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(ClassS);
            xdata{m,i} = VBResults.x.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(ClassS);
        catch
            ydata{m,i} = VBResults.AQ.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(Class);
            xdata{m,i} = VBResults.x.(Type{i}).(SubType{i}).(Width{i}).(Layout).(Support{i}).(Trans{i}).(AE(m)).(Traffic).(Class);
        end
        % DisplayName
        PDets.DN{m,i} = [Type{i} + " " + SubType{i} + " " +  Support{i} + " " + Trans{i}];
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDets.MEC{m}(i,:) = [0 0 0]; 
        PDets.MFC{m}(i,:) = D(i*4,:);
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIMComp',FigNum,FName); clear ydata xdata PDets
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end
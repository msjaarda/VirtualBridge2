clear, clc%, close all

% Global Definitions
FigNum = 0; D = linspecer(30);      % Set colours
Type = 'Box';         SubType = 'Stand';
Support = 'Simp';     Trans = 'p0';
% Subplot Info (m)
AE(1) = "Mn"; AE(2) = "Mp"; AE(3) = "V";
Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
% Get VBResults using VBOutput2Struct
[VBResultsPx] = VBOutput2Struct('PlatooningFull');
load('VBResults.mat')
% Class
Class = 'ClassOW';

Export = true;

% BOX --------------------------------------------------------------------

% General Info
Width = 'Wid12';      Layout = 'Uni';
% Series Info (i)
Traffic = fieldnames(VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).Mn);
% Platooning Info
PSize = 'Size3'; PRate = 'Rate2'; PFolDist{1} = 'FolDist50'; PFolDist{2} = 'FolDist50'; PFolDist{3} = 'FolDist75'; PFolDist{4} = 'FolDist100';

FName = 'Box Girder Bridge Uni 2L w 3 Tr Platoons @ 20% (Class+)';
for m = 1:3
    for i = 1:4
        % Get ydata and xdata
        if i == 1
            ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{1}).(Class);
            PDets.MFC{m}(i,:) = [0 0 0];
            PDets.DN{m,i} = Traffic{1};
        elseif i == 2
            ydata{m,i} = VBResultsPx.P.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{1}).(PSize).(PRate).(PFolDist{i}).(Class);
            PDets.MFC{m}(i,:) = D(i*4,:);
            PDets.DN{m,i} = PFolDist{i};
        elseif i == 3
            ydata{m,i} = VBResultsPx.P.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{1}).(PSize).(PRate).(PFolDist{i}).(Class);
            PDets.MFC{m}(i,:) = D(i*4,:);
            PDets.DN{m,i} = PFolDist{i};
        else
            ydata{m,i} = VBResultsPx.P.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{1}).(PSize).(PRate).(PFolDist{i}).(Class);
            PDets.MFC{m}(i,:) = D(i*4,:);
            PDets.DN{m,i} = PFolDist{i};
        end
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{1}).(Class);
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5;
        PDets.MEC{m}(i,:) = [0 0 0];
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIM',FigNum,FName); 
if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end

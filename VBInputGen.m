% ------------------------------------------------------------------------
%                              VBInputGen
% ------------------------------------------------------------------------
% Summarize traffic characteristics from WIM or VirtualWIM
%       - Output is a MATSimInput spreadsheet to be used with MATSim
%       - Class+ Option Added 19/8/2022

tic, clear, clc, hold off, close all, format long g, rng('shuffle')

% ----- INPUT -----

% Year, #, Station Name, string, save and plot toggles
%Year = 2017; SName = 'Ger9625'; Save = 0; PlotFits = 1; Stage2P = 1;
Year = 2019; SName = 'Ceneri'; PlotFits = 1; Save = 0; Stage2P = 1; ClassPlus = 1;

% ----- ENDIN -----

load('Sites.mat')
Sitesx = Sites.SITE(strcmp(Sites.SName,SName) & datetime(Year,1,1) >= Sites.StartDate & datetime(Year,1,1) <= Sites.EndDate & Sites.Core);

if length(Sitesx) == 1
    load(['WIM\' num2str(Sitesx) '.mat']) 
elseif length(Sitesx) == 2
    load(['WIM\' num2str(Sitesx(1)) '.mat'])
    PDs = PDs(year(PDs.DTS) == Year,:);  PD1 = PDs;
    load(['WIM\' num2str(Sitesx(2)) '.mat'])
    PDs = PDs(year(PDs.DTS) == Year,:);  PDs = [PD1; PDs];
else; disp('Error: WIM > 3'); end

% Stage 2 Prune
if Stage2P; PDs = Stage2Prune(PDs); end

% Get Input file!

[TrTyps,TrName] = VBTypes2Names;
TrName = cellstr(TrName); TrName(end-1:end) = []; TrTyps(end-1:end) = [];

if ~ClassPlus
    TrTyps(end-3:end) = [];  TrName(end-3:end) = [];
    TrAxPerGr = [11; 12; 22; 23; 111; 1111; 112; 1211; 122; 1112; 112; 113; 123];
    TrTypPri = [21; 21; 21; 21; 321; 2341; 231; 2341; 231; 2341; 321; 321; 321];
else
    TrAxPerGr = [11; 12; 22; 23; 111; 1111; 112; 1211; 122; 1112; 112; 113; 123; 14; 15; 223; 332];
    TrTypPri = [21; 21; 21; 21; 321; 2341; 231; 2341; 231; 2341; 321; 321; 321; 21; 21; 321; 321];
    del = [];
    for j = 1:4
        if sum(PDs.CLASS == TrTyps(13+j)) < 10
            del = [del; 13+j];
        end
    end
    TrTyps(del) = []; TrName(del) = []; TrAxPerGr(del) = []; TrTypPri(del) = [];
end

% Get three tables
BetAx_Excel = BetAx(PDs,TrName,TrTyps,TrAxPerGr,TrTypPri,SName,Year,PlotFits);
LinFit_Excel = LinFit(PDs,TrName,TrTyps,TrAxPerGr,TrTypPri,SName,Year,PlotFits);
Distr_Excel = Distr(PDs,TrName,TrTyps,TrAxPerGr,TrTypPri,SName,Year,PlotFits);
% Summarize Weight by Axles
[STaTr,AllAx] = AxleStats(PDs,TrAxPerGr,TrTyps,SName,Year,PlotFits);

load('TrLib.mat')
TrLib.([SName,num2str(Year)]).TrDistr = Distr_Excel;
TrLib.([SName,num2str(Year)]).TrLinFit = LinFit_Excel;
TrLib.([SName,num2str(Year)]).TrAllo = TrLib.Ceneri2018.TrAllo;
TrLib.([SName,num2str(Year)]).TrBetAx = BetAx_Excel;
TrLib.([SName,num2str(Year)]).TrWitAx = TrLib.Ceneri2018.TrWitAx;
%TrLib.([SName,num2str(Year)]).TrDistr.PlatPct = [];
TrLib = orderfields(TrLib);

% If importing new TrLib to TrLib
if Save == 1
    save('Misc/TrLib.mat','TrLib')
end

% Optional
%LaneDistBr(PDC,TrTyps,TrAxPerGr,Station)
%ExportOpen

% ------------------------------------------------------------------------
%                              VBInputGen
% ------------------------------------------------------------------------
% Summarize traffic characteristics from WIM or VirtualWIM
%       - Output is a MATSimInput spreadsheet to be used with MATSim
%       - Class+ Option Added 19/8/2022

tic, clear, clc, hold off, close all, format long g, rng('shuffle'), load('TrLib.mat')

% ----- INPUT -----

% Year, #, Station Name, string, save and plot toggles
%Year = 2017; SName = 'Ger9625'; Save = 0; PlotFits = 1; Stage2P = 1;
Year = 2019; SName = 'Gotthard'; PlotFits = 0; Save = 1; Stage2P = 1; ClassPlus = 1;

% ----- ENDIN -----

load('Sites.mat')
Sitesx = Sites.SITE(strcmp(Sites.SName,SName) & datetime(Year,1,1) >= Sites.StartDate & datetime(Year,1,1) <= Sites.EndDate & Sites.Core);

if length(Sitesx) == 1
    load(['WIM/' num2str(Sitesx) '.mat']) 
elseif length(Sitesx) == 2
    load(['WIM/' num2str(Sitesx(1)) '.mat'])
    PDs = PDs(year(PDs.DTS) == Year,:);  PD1 = PDs;
    load(['WIM/' num2str(Sitesx(2)) '.mat'])
    PDs = PDs(year(PDs.DTS) == Year,:);  PDs = [PD1; PDs];
else; disp('Error: WIM > 3'); end

% Stage 2 Prune
if Stage2P; PDs = Stage2Prune(PDs); end

% Get Input file!

[TrTyps,TrName] = VBTypes2Names;
TrName = cellstr(TrName); TrName(end-1:end) = []; TrTyps(end-1:end) = [];
TrAllo = TrLib.Ceneri2018.TrAllo; TrWitAx = TrLib.Ceneri2018.TrWitAx;

if ~ClassPlus
    TrTyps(end-3:end) = [];  TrName(end-3:end) = [];
    TrAxPerGr = [11; 12; 22; 23; 111; 1111; 112; 1211; 122; 1112; 112; 113; 123];
    TrTypPri = [21; 21; 21; 21; 321; 2341; 231; 2341; 231; 2341; 321; 321; 321];
else
    load('Misc/Cranes.mat')
    TrAxPerGr = [11; 12; 22; 23; 111; 1111; 112; 1211; 122; 1112; 112; 113; 123; 14; 15; 223; 332];
    TrTypPri = [21; 21; 21; 21; 321; 2341; 231; 2341; 231; 2341; 321; 321; 321; 21; 21; 321; 321];
    TrAllo.j = repmat(nan,height(TrAllo),1);
    TrAllo.a(14:15) = 1; TrAllo.e(14) = .25; TrAllo.f(14) = .25; TrAllo.g(14) = .25; TrAllo.h(14) = .25;
    TrAllo.f(15) = .2; TrAllo.g(15) = .2; TrAllo.h(15) = .2; TrAllo.i(15) = .2; TrAllo.j(15) = .2;
    TrAllo.a(16) = 0.5; TrAllo.b(16) = 0.5; TrAllo.d(16) = 0.5; TrAllo.e(16) = 0.5; 
    TrAllo.g(16) = 0.33; TrAllo.h(16) = 0.34; TrAllo.i(16) = 0.33;
    TrAllo.a(17) = 0.3; TrAllo.b(17) = 0.34; TrAllo.c(17) = 0.33;
    TrAllo.d(17) = 0.3; TrAllo.e(17) = 0.34; TrAllo.f(17) = 0.33;
    TrAllo.g(17) = 0.5; TrAllo.h(17) = 0.5;
    TrWitAx = renamevars(TrWitAx,'EndD','d67');
    TrWitAx.d78 = repmat(nan,height(TrWitAx),1); TrWitAx.EndD = repmat(nan,height(TrWitAx),1);
    TrWitAx.d23(14) = mean(Crane.t60.W2); TrWitAx.d34(14) = mean(Crane.t60.W3); TrWitAx.d45(14) = mean(Crane.t60.W4);
    TrWitAx.d23(15) = mean(Crane.t72.W2); TrWitAx.d34(15) = mean(Crane.t72.W3); TrWitAx.d45(15) = mean(Crane.t72.W4); TrWitAx.d56(15) = mean(Crane.t72.W5);
    TrWitAx.d23(16) = mean(Crane.t84.W2); TrWitAx.d34(16) = mean(Crane.t84.W3); TrWitAx.d45(16) = mean(Crane.t84.W4); TrWitAx.d56(16) = mean(Crane.t84.W5); TrWitAx.d67(16) = mean(Crane.t84.W6);
    TrWitAx.d23(17) = mean(Crane.t96.W2); TrWitAx.d34(17) = mean(Crane.t96.W3); TrWitAx.d45(17) = mean(Crane.t96.W4); TrWitAx.d56(17) = mean(Crane.t96.W5); TrWitAx.d67(17) = mean(Crane.t96.W6); TrWitAx.EndD(17) = mean(Crane.t96.W7);
    del = [];
    for j = 1:4
        if sum(PDs.CLASS == TrTyps(13+j)) < 10
            del = [del; 13+j];
        end
    end
    TrTyps(del) = []; TrName(del) = []; TrAxPerGr(del) = []; TrTypPri(del) = []; TrAllo(del,:) = [];
    TrWitAx(del,:) = [];
end

% Get three tables
BetAx_Excel = BetAx(PDs,TrName,TrTyps,TrAxPerGr,TrTypPri,SName,Year,PlotFits);
LinFit_Excel = LinFit(PDs,TrName,TrTyps,TrAxPerGr,TrTypPri,SName,Year,PlotFits);
Distr_Excel = Distr(PDs,TrName,TrTyps,TrAxPerGr,TrTypPri,SName,Year,PlotFits);
% Summarize Weight by Axles
[STaTr,AllAx] = AxleStats(PDs,TrAxPerGr,TrTyps,SName,Year,PlotFits);

if ClassPlus
    NAME = [SName,num2str(Year),'CP'];
else
    NAME = [SName,num2str(Year)];
end

TrLib.(NAME).TrDistr = Distr_Excel;
TrLib.(NAME).TrLinFit = LinFit_Excel;
TrLib.(NAME).TrBetAx = BetAx_Excel;
TrLib.(NAME).TrAllo = TrAllo;
TrLib.(NAME).TrWitAx = TrWitAx;
TrLib.(NAME).TrDistr.PlatPct = [];
TrLib = orderfields(TrLib);

% If importing new TrLib to TrLib
if Save == 1
    save('Misc/TrLib.mat','TrLib')
end

% Optional
%ExportOpen

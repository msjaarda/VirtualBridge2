clear all, clc

%% Script to add to existing VBResults table, OutInfo deterministics results

%% INPUT
Folder_Names = 'Det'; % Name of the folder inside Output, containing deterministic results 
VBResults_Names = 'VBResults'; % Name of the VBResults matlab file where we will add the deterministic results
DetTraffic_Names = 'Deterministic'; % Name of the deterministic truck(s) inside table in Traffic column

%% CODE
%NameFileSave = append('VBResults.mat');
NameFileSave = append(VBResults_Names,'.mat');
load(append(VBResults_Names,'.mat'));
warning('off','MATLAB:table:RowsAddedExistingVars');

LenPrint = []; RamUsed = [];

% Open folders
Dir_List = dir(append('Output/',Folder_Names));
File_List = {Dir_List(:).name}';

%cleaning file list
File_List = File_List(~strcmp(File_List,'.')&~strcmp(File_List,'..')&contains(File_List,'.mat'));
File_List = erase(File_List,'.mat');

% Start Progress Bar
u = StartProgBar(height(File_List), 1, 1, 1); tic; st = now;

for i=1:height(File_List)
    
    load(append('Output/',Folder_Names,'/',File_List{i},'.mat'));
    
    NameCodes = fieldnames(VBResults);
    NameCodes(strcmp(NameCodes,'AQ')) = []; % Remove AQ struct containing old format results
    
    E = VBGetECode(OutInfo.ILData(:),OutInfo.BaseData.ILRes(1));
    ILNames = string({OutInfo.ILData.Name}');
    for j=1:length(ILNames)
        ILSplit(j,:) = strsplit(ILNames(j),'.');
    end
    %Delete first column (ILLib)
    ILSplit(:,1) = [];
    %Remove "S"s from spans
    ILSplit(:,8) = extractAfter(ILSplit(:,8),1);
    
        for k=1:height(NameCodes)
            BlockMax = fieldnames(VBResults.(NameCodes{k}));
            for l=1:height(BlockMax)
                VBResults.(NameCodes{k}).(BlockMax{l}).Type(end+1:end+length(ILNames)) = ILSplit(:,1);
                VBResults.(NameCodes{k}).(BlockMax{l}).SubType(end-length(ILNames)+1:end) = ILSplit(:,2);
                VBResults.(NameCodes{k}).(BlockMax{l}).Width(end-length(ILNames)+1:end) = ILSplit(:,3);
                VBResults.(NameCodes{k}).(BlockMax{l}).Layout(end-length(ILNames)+1:end) = ILSplit(:,4);
                VBResults.(NameCodes{k}).(BlockMax{l}).Support(end-length(ILNames)+1:end) = ILSplit(:,5);
                VBResults.(NameCodes{k}).(BlockMax{l}).Trans(end-length(ILNames)+1:end) = ILSplit(:,6);
                VBResults.(NameCodes{k}).(BlockMax{l}).AE(end-length(ILNames)+1:end) = ILSplit(:,7);
                VBResults.(NameCodes{k}).(BlockMax{l}).Traffic(end-length(ILNames)+1:end) = DetTraffic_Names;
                VBResults.(NameCodes{k}).(BlockMax{l}).Span(end-length(ILNames)+1:end) = str2double(ILSplit(:,8));
                VBResults.(NameCodes{k}).(BlockMax{l}).BestFitAll(end-length(ILNames)+1:end) = 'None';
                VBResults.(NameCodes{k}).(BlockMax{l}).BestFitClassOW(end-length(ILNames)+1:end) = 'None';
                VBResults.(NameCodes{k}).(BlockMax{l}).BestFitClass(end-length(ILNames)+1:end) = 'None';
                for j=1:length(ILNames)
                VBResults.(NameCodes{k}).(BlockMax{l}).Class(end-length(ILNames)+j) = OutInfo.OverMax(j,2)./E(j).(NameCodes{k}).Total;
                end
                VBResults.(NameCodes{k}).(BlockMax{l}).ClassOW(end-length(ILNames)+1:end) = VBResults.(NameCodes{k}).(BlockMax{l}).Class(end-length(ILNames)+1:end);
                VBResults.(NameCodes{k}).(BlockMax{l}).All(end-length(ILNames)+1:end) = VBResults.(NameCodes{k}).(BlockMax{l}).Class(end-length(ILNames)+1:end);
            end
        end
    
    clear ILSplit
    % Update progress bar
    user = memory;
    RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
    LenPrint = VBUpProgBar(st,RamUsed(end),i,LenPrint);
end

save(NameFileSave,'VBResults');


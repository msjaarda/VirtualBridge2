%  April 2019 | Matthew Sjaarda
%  Take daily traffic files *.V00 from OFROU, combine into a yearly *.trd
%  Keep in mind one must rename the header row in a text editor because
%  matlab cannot write hyphens, and the .trd file has hyphens
%  Say 'Yes' to pruning if the file hasn't already been compiled. If it
%  has, use the PruneYearlymatFile instead.

% April 2021 | Matthew Sjaarda
% Now we should do the same thing, but this time we make a few structural
% changes... we no longer need to maintain compatibility with previous
% versions, since we have ALL the raw data!
% Could be .V00 or .VA0 or .V01 or .V0#
% 1) We do not have to extract the files from their folders
% 2) We will name the Raw variables as we name the Processed (without Station ID)
% 3) We will fix 'Headx' and other issues
% 4) We will use Datetime rather than separate columns
% see Swiss WIM Data for column details

clear, clc, close all, format long g
fprintf('\n Time: %s \n\n',datetime('now')), tic

% -------- INPUT --------

% Tell Matlab where the directory is that holds the daily files
Directory = 'C:\Users\mjsja\switchdrive\WIMDrop';
%E:\Rohdaten
Folder_List = dir(Directory);

% Get List of Folders, including SNames (strings) and Stations (nums)

% Initialize count
count = 1;
for i = 1:length(Folder_List)
    try
        % Folder must be of M660 variety
        if strcmp(Folder_List(i).name(end-4:end-1),'M660')
            KeepFolder(i) = true;
            % Name is 7th to first space or comma
            FC = strfind(Folder_List(i).name,',');
            SNames(count) = string(Folder_List(i).name(7:FC(1)-1));
            % Station is first 3 digits
            Stations(count) = str2num(Folder_List(i).name(1:3));
            % Increase count
            count = count + 1;
        else
            KeepFolder(i) = false;
        end
    catch
        KeepFolder(i) = false;
    end
end
% Refine folder list to only those to be processed
Folder_List = Folder_List(KeepFolder); SName = SNames'; Stations = Stations';

% Loop through each folder
for i = 6%:length(Folder_List)
    
    % Identify subfolders
    SFolder_List = dir([Folder_List(i).folder '\' Folder_List(i).name]);
    SFolder_List(1:2) = [];
    
    % Loop through subfolders, recording the Year in question
    for j = 1:length(SFolder_List)
        % Now we go through each individual folder looking for .VA0, .V00, .V01
        Files = [dir([SFolder_List(j).folder '\' SFolder_List(j).name '\**\*.VA0']);...
            dir([SFolder_List(j).folder '\' SFolder_List(j).name '\**\*.V0*'])];
        FilesEx = [dir([SFolder_List(j).folder '\' SFolder_List(j).name '\**\*.csv'])];
        
        % Get total number of files
        NumFiles = length(Files);
        
        % Get Station, SName, and Year
        Station = Stations(i);
        SName = SNames(i);
        Year = str2double(SFolder_List(j).name(end-3:end));
        
        % Create final file name
        FileName = strcat(SName,'_',num2str(Year));
        
        RDI = [];
        
        % Loop through all files to read and concatenate
        parfor k = 1:NumFiles
            
            % Have to do try catch because some files may be courrupted
            % Example: Gotthard June 21, 2012
            try
                
                % Get individual filename
                FullFileName = strcat(Files(k).folder,'\',Files(k).name);
                
                % Open file
                fileID = fopen(FullFileName,'r');
                
                formatSpec = '%8s%[^\n\r]';
                startRow = 1;
                endRow = 22;
                dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                % Reset filepointer
                frewind(fileID)
                % Define format
                if strcmp(dataArray{1,1}(15),'* UNITS ')
                    formatSpec = '%5s%3s%8s%3s%2s%2s%3s%2s%3s%3s%7s%3s%2s%5s%5s%4s%6s%3s%3s%6s%5s%5s%5s%5s%5s%5s%5s%5s%6s%6s%6s%6s%6s%6s%6s%6s%6s%s%[^\n\r]';
                    HR = 21;
                else
                    formatSpec = '%5s%3s%8s%3s%2s%2s%3s%2s%3s%3s%7s%3s%2s%5s%5s%4s%6s%3s%3s%6s%5s%5s%5s%5s%5s%5s%5s%5s%6s%5s%5s%5s%5s%5s%5s%5s%5s%s%[^\n\r]';
                    HR = 20;
                end
                
                % Read data
                dataArrayBlock = textscan(fileID, formatSpec, Inf, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', HR, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                if size(dataArrayBlock{1},1) ~= size(dataArrayBlock{end},1)
                    frewind(fileID)
                    dataArrayBlock = textscan(fileID, formatSpec, size(dataArrayBlock{end},1), 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', HR, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                end
                % Define inline function
                fn = @(x) cellfun(@str2doubleq,x(1:end-1));
                % Create an array from the read
                RDi = cell2mat(cellfun(fn,dataArrayBlock(3:end),'UniformOutput',false));
                % Close file
                fclose(fileID);
                
                % Get size of RDi
                [RDisize, RDirows] = size(RDi);
                
                % If very small, skip
                if RDisize <= 1
                    continue
                end
                
                RDI = [RDI; RDi];
                
            catch
                continue
            end
            
        end
        
        if isempty(RDI)
            continue
        end
        
        if istable(RDI)
            
            RD = table();
            RD.SITE = Station*ones(height(RDI),1);

            D = mod(RDI.YYYYMMDD,100);
            Mo = (mod(RDI.YYYYMMDD,10000) - D)/100;
            Y = (RDI.YYYYMMDD - mod(RDI.YYYYMMDD,10000))/10000;
            S = mod(RDI.HHMMSS,100);
            M = (mod(RDI.HHMMSS,10000) - S)/100;
            H = (RDI.HHMMSS - mod(RDI.HHMMSS,10000))/10000;
            
            RD.DTS = datetime(Y,Mo,D,H,M,S);
            
            RDI.VarName9 = [];
            RDI.W9_10 = [];
            RDI.AWT10 = [];
            RDI{:,11:19} = RDI{:,11:19}*10;
            
            RD = [RD RDI(:,4:end)];
            RD.SPEED = RD.SPD/100;
            RD = sortrows(RD,2);
            
        else
            
            RDI = real(RDI);
            % RD will be the final variable to be saved (table)
            RD = table;
            % Create SITE and COUNTID
            RD.SITE = Station*ones(size(RDI,1),1);
            % Create datetime variable
            RD.DTS = datetime(Year*ones(size(RDI,1),1),RDI(:,3),RDI(:,2),RDI(:,5),RDI(:,6),RDI(:,7),10*RDI(:,8));
            RD.COUNTID = RDI(:,1);
            % Delete column 18 (total wheelbase) and 1-9 (RESCOD and others)
            RDI(:,18) = []; RDI(:,1:9) = [];
            % Reorder RDI
            RDI(:,18:26) = RDI(:,18:26)*10;
            RDI = RDI(:,[1 2 5 7 17 3 4 6 8 18 19 20 21 22 23 24 25 26 9 10 11 12 13 14 15 16]);
            % Create RD
            RD = [RD array2table(RDI,'VariableNames',{'LANE' 'DIR' 'SPEED' 'AX' 'GW_TOT' 'HEADT' 'GAPT' 'LENTH' 'CS' 'AWT01'...
                'AWT02' 'AWT03' 'AWT04' 'AWT05' 'AWT06' 'AWT07' 'AWT08' 'AWT09' 'W1_2' 'W2_3' 'W3_4' 'W4_5'...
                'W5_6' 'W6_7' 'W7_8' 'W8_9'})];
            
        end
        
        RD.GW_TOT = RD.GW_TOT*10;
        RD = sortrows(RD,2);
        RD(RD.GW_TOT < 1000,:) = [];
        RD(isnan(RD.COUNTID),:) = [];
        RD(isnan(RD.LANE),:) = [];
        RD(isnan(RD.LANE),:) = [];
        % Look for folder (SName) in Raw WIM
        if ~isfolder(strcat('Raw WIM\',SName))
            % Add if not there
            mkdir(strcat('Raw WIM\',SName))
        end
        
        % Look for pre-existing file with FileName
        if ~isfile(strcat('Raw WIM\',SName,'\',FileName,'.mat'))
            save(strcat('Raw WIM\',SName,'\',FileName),'RD')
        else
            RDx = RD;
            RDx.DIR = zeros(length(RDx.DIR),1);
            % If there, append
            load(strcat('Raw WIM\',SName,'\',FileName,'.mat'))
            RD = [RD; RDx];
            RD = sortrows(RD,2);
            save(strcat('Raw WIM\',SName,'\',FileName),'RD')
        end
        
        
        clear RD
        fprintf('%s %d %d %s mins %.0f secs\n',SName,Station,Year,num2str(floor(toc/60)),rem(toc,60))
        tic
    end
end
fprintf('\n Time: %s \n\n',datetime('now'))

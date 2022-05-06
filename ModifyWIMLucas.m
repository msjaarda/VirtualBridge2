clear, clc, close all

% Script to change the fitting method. You can open previous results and
% then save them with a new fitting method

% Folder from where to load previous results
OpenFolder = 'WIM60tv10MS';
% Folder to save the new results
SaveFolder = 'WIM10MSNormalWithZero';

Dir_List = dir(append('Output\',OpenFolder,'\'));
AllName = {Dir_List(:).name};

%cleaning file list
AllName = AllName(~strcmp(AllName,'.')&~strcmp(AllName,'..')&contains(AllName,'.mat'));
%File_List = erase(File_List,'.mat');

%AllName = {'Nov27-21 011254.mat'};
for z = 1:width(AllName)
    Name = AllName{z};
    load(append('Output\',OpenFolder,'\',Name))
    %load('Output\WIM\Nov01-21 130953.mat')
    
    for  i = 1:width(OutInfo.ILData) %i = [16:20,41:45,66:70,82:83,92:93,102:103,121:125,146:150,171:175] %i = [16:20,41:45,66:70,91:95,116:120,141:145]
        % Call it Out1 for now (just in case we load 2 OutInfos, for comparison)
        Out1 = OutInfo;
        
        % We set the BlockM we want, and the Period
        BlockM = {'ClassOW'};%{'All','ClassOW','Class'}; %{'ClassOW'};           % Could also be 'ClassOW' or 'All'
        Period = {'Yearly','Monthly'}; %{'Weekly'}; %Period = {'Yearly','Monthly'};          % Could also be 'Yearly' or 'Daily' or 'Monthly'
        % Set influence line you want (take a look at OutInfo.ILData to know which
        AE = i; %10;
        
        for j=1:width(Period)
            for k =1:width(BlockM)
                % Get a subset of the maximum results for the three selections made
                Max1 = OutInfo.Max(AE).(BlockM{k}).(Period{j});
                
                if height(Max1) == 0
                    
                else
                    
                    % Check the fit
                    %[~,~,~,~] = GetBlockMaxFit(Max1.Max,'Lognormal',1);
                    % Calculate Ed and AQ (can compare
                    %[Ed, AQ, ~] = GetBlockMaxEd(Max1,Period,'Lognormal',Out1.ESIA.Total(AE),Out1.ESIA.EQ(:,AE),Out1.ESIA.Eq(:,AE),.7,.6);
                    %Ed1 = Ed;
                    if height(Max1.Max) > 35
                        [Ed, AQ, ~] = GetBlockMaxEd(Max1.Max,Period{j},'NormalZero',Out1.ESIA.Total(AE),Out1.ESIA.EQ(:,AE),Out1.ESIA.Eq(:,AE),.7,.6,Out1.PropTrucks.(BlockM{k}).(Period{j})(i),1);
                        if ishandle(1)
                            title(Out1.ILData(i).Name);
                        end
                        %if Ed(4)>Ed1
                        %   go = 1;
                        %end
                        %OutInfo.EdLN.(BlockM{j}).(Period)(i) = Ed(4);
                        %Lucas = OutInfo.EdLN.(BlockM{j}).(Period)(i)-Ed;
                        OutInfo.EdLN.(BlockM{k}).(Period{j})(i) = Ed;
                    end
                end
                
            end
        end
    end
    
    Dir_List = dir('Output\');
    File_List = {Dir_List.name}';
    if sum(strcmp(File_List,SaveFolder))>=1
    else
        mkdir(append('Output\',SaveFolder));
    end
    
    save(append('Output\',SaveFolder,'\',Name),'OutInfo');
    disp(append(int2str(z),'/',int2str(width(AllName))));
end

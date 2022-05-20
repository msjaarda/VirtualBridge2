clear, clc, close all

% Script to change the fitting method. You can open previous results and then save them with a new fitting method

% Folder from where to load previous results
OpenFolder = 'WIM60tv10MS';
% Folder to save the new results
SaveFolder = 'WIM10MSNormalWithZero';

Dir_List = dir(append('Output\',OpenFolder,'\'));
AllName = {Dir_List(:).name};

%cleaning file list
AllName = AllName(~strcmp(AllName,'.')&~strcmp(AllName,'..')&contains(AllName,'.mat'));
%File_List = erase(File_List,'.mat');

for z = 1:width(AllName)
    Name = AllName{z};
    load(append('Output\',OpenFolder,'\',Name))
    %load('Output\WIM\Nov01-21 130953.mat')
    
    for  i = 1:width(OutInfo.ILData)
        % Call it Out1 for now (just in case we load 2 OutInfos, for comparison)
        Out1 = OutInfo;
        
        % We set the BlockM we want, and the Period
        BlockM = {'ClassOW'};
        Period = {'Yearly','Monthly'};
        % Set influence line you want (take a look at OutInfo.ILData to know which
        AE = i; %10;
        
        for j=1:width(Period)
            for k =1:width(BlockM)
                % Get a subset of the maximum results for the three selections made
                Max1 = OutInfo.Max(AE).(BlockM{k}).(Period{j});
                
                if height(Max1) == 0
                    
                else
                    if height(Max1.Max) > 35
                        [Ed, AQ, ~] = GetBlockMaxEd(Max1.Max,Period{j},'NormalZero',Out1.ESIA.Total(AE),Out1.ESIA.EQ(:,AE),Out1.ESIA.Eq(:,AE),.7,.6,Out1.PropTrucks.(BlockM{k}).(Period{j})(i),1);
                        if ishandle(1)
                            title(Out1.ILData(i).Name);
                        end
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

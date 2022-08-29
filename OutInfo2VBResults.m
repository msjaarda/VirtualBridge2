clear all, clc

%% Script to create VBResults table from OutInfo results

%% INPUT : Open folder who contains WIM results
% First place is for the WIM or SIM folder,
% Second is for the models folder (if needed)

Folder_Names{1} = 'WIMv17pr';
Folder_Names{2} = Folder_Names{1};
%Folder_Names{2} = 'WIMDynPO'; %second folder will import the ECodes of the 2sd file

% Select parameters for alpha analysis
AlphaAnalys = 4; % 1)Blended Alpha 2)AlphaQ1 3)AlphaQ2 4)Alphaq
AlphaQ1 = 0.6;
AlphaQ2 = 0.4;
Alphaq = 0.5;

%% Running script
%open folders
Dir_List = dir('Output'); % Folder containing the Output results
File_List = {Dir_List.name}';

% check if WIM folder exist
if sum(strcmp(File_List,Folder_Names{1}))>=1
else
disp("------------------NO FOLDER IN OUTPUT TRY AGAIN----------------------");
end

% check fast loop
if strcmp(Folder_Names{1},Folder_Names{2})
   FasterLoop = 1;
else
   FasterLoop = 0;
end

if AlphaAnalys == 1
   NameFileSave = append('VBResults.mat');
elseif AlphaAnalys == 2
   NameFileSave = append('VBResultsAlphaQ1(Q2=',num2str(AlphaQ2),'etq=',num2str(Alphaq),').mat');
elseif AlphaAnalys == 3
   NameFileSave = append('VBResultsAlphaQ2(Q1=',num2str(AlphaQ1),'etq=',num2str(Alphaq),').mat');
elseif AlphaAnalys == 4
   NameFileSave = append('VBResultsAlphaq(Q1=',num2str(AlphaQ1),'etQ2=',num2str(AlphaQ2),').mat');
end

%set factors
LambdaQ = 1.5; %SIA
%AlphaQ1 = 0.6; %Fixed at 0.7/0.6
%AlphaQ2 = 0.4; %Fixed at 0.5/0.4
%Alphaq = 0.5; %editable 0.5

LenPrint = []; RamUsed = [];

%% Read inside each WIM reasults and save new format for analyse, same for BoxSimNuma

Dir_List = dir(append('Output/',Folder_Names{1}));
File_List = {Dir_List(:).name}';

%cleaning file list
File_List = File_List(~strcmp(File_List,'.')&~strcmp(File_List,'..')&contains(File_List,'.mat'));
File_List = erase(File_List,'.mat');

% Start Progress Bar
u = StartProgBar(height(File_List), 1, 1, 3-FasterLoop); tic; st = now;

for i=1:height(File_List)
    
    load(append('Output/',Folder_Names{1},'/',File_List{i},'.mat'));
    
    if OutInfo.BaseData.StopSim == 0
    Fields1 = fieldnames(OutInfo.pd);
    
    ILDataNames = {OutInfo.ILData.Name}';
    
    for j=1:height(Fields1)
    
    templ = OutInfo.pd.(Fields1{j});
    Fields2 = fieldnames(templ);
    
        for k=1:height(Fields2)
            
            for l=1:height(ILDataNames)
        
        % Create and save inside CombInfo struct
        BestFit = OutInfo.pd(l).(Fields1{j}).(Fields2{k}).Best;
        MatName = append(ILDataNames{l},'.EdLN.',Fields1{j},'.',Fields2{k},' = ','OutInfo.pd(',int2str(l),').',Fields1{j},'.',Fields2{k},'.','(BestFit)','.Ed;');
        MatName = erase(MatName,'ILLib');
        MatName = append('CombInfo',MatName);
        eval(MatName);
        % Save best fit
        MatName = append(ILDataNames{l},'.EdLN.',Fields1{j},'.',Fields2{k},'BestFit',' = BestFit;');
        MatName = erase(MatName,'ILLib');
        MatName = append('CombInfo',MatName);
        eval(MatName);
        % Saving prop trucks
        MatName = append(ILDataNames{l},'.PropTrucks.',Fields1{j},'.',Fields2{k},' = ','OutInfo.PropTrucks.',Fields1{j},'.',Fields2{k},'(',int2str(l),');');
        MatName = erase(MatName,'ILLib');
        MatName = append('CombInfo',MatName);
        try
        eval(MatName);        
        catch
        end
            end
            
        end
        
    end
    
    if FasterLoop
        Comp_List = fieldnames(OutInfo.E);
        ILDataNames = {OutInfo.ILData.Name}';
        
        for j=1:height(Comp_List)
            
            temp2 = OutInfo.E.(Comp_List{j});
            Fields1 = fieldnames(temp2);
            
            for k=1:height(Fields1)
                
                for l=1:height(ILDataNames)
                    
                    % Create and save inside CombInfo struct
                    MatName = append(ILDataNames{l},'.',Comp_List{j},'.',Fields1{k},' = ','OutInfo.E','(',int2str(l),')','.',Comp_List{j},'.',Fields1{k},';');
                    MatName = erase(MatName,'ILLib');
                    MatName = append('CombInfo',MatName);
                    eval(MatName);
                    
                    % Check if EdLN exist with EModels, if yes do the ratio
                    if isfield(eval(append('CombInfo',erase(ILDataNames{l},'ILLib'))),'EdLN') && strcmp(Fields1{k},'Total')
                        
                        Fields2 = fieldnames(eval(append('CombInfo',erase(ILDataNames{l},'ILLib'),'.EdLN')));
                        
                        for m=1:height(Fields2)
                            
                            Fields3 = fieldnames(eval(append('CombInfo',erase(ILDataNames{l},'ILLib'),'.EdLN.',Fields2{m})));
                            
                            for n=1:height(Fields3)
                                
                                Actual = append('CombInfo',erase(ILDataNames{l},'ILLib'),'.EdLN',Fields2{m},Fields3{n},'0',Comp_List{j},' = ','CombInfo',erase(ILDataNames{l},'ILLib'),'.EdLN.',Fields2{m},'.',Fields3{n},'/','CombInfo',erase(ILDataNames{l},'ILLib'),'.',Comp_List{j},'.Total;');
                                eval(Actual);
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
    end
    
    end
    % Update progress bar
    user = memory;
    RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
    LenPrint = VBUpProgBar(st,RamUsed(end),i,LenPrint);
end

if FasterLoop == 0
% Go inside Model folder
Dir_List = dir(append('Output/',Folder_Names{2}));
File_List = {Dir_List(:).name}';

%cleaning file list
File_List = File_List(~strcmp(File_List,'.')&~strcmp(File_List,'..')&contains(File_List,'.mat'));
File_List = erase(File_List,'.mat');

% Start Progress Bar
u = StartProgBar(height(File_List), 1, 2, 3); tic; st = now;

for i=1:height(File_List)
    
    load(append('Output/',Folder_Names{2},'/',File_List{i},'.mat'));
    if OutInfo.BaseData.StopSim == 0
    Comp_List = fieldnames(OutInfo.E); 
    ILDataNames = {OutInfo.ILData.Name}';
            
    for j=1:height(Comp_List)
    
    temp2 = OutInfo.E.(Comp_List{j});
    Fields1 = fieldnames(temp2);
    
        for k=1:height(Fields1)
                      
            for l=1:height(ILDataNames)
        
        % Create and save inside CombInfo struct
        MatName = append(ILDataNames{l},'.',Comp_List{j},'.',Fields1{k},' = ','OutInfo.E','(',int2str(l),')','.',Comp_List{j},'.',Fields1{k},';');
        MatName = erase(MatName,'ILLib');
        MatName = append('CombInfo',MatName);
        eval(MatName);
        
        % Check if EdLN exist with EModels, if yes do the ratio
        if isfield(eval(append('CombInfo',erase(ILDataNames{l},'ILLib'))),'EdLN') && strcmp(Fields1{k},'Total')
            
            Fields2 = fieldnames(eval(append('CombInfo',erase(ILDataNames{l},'ILLib'),'.EdLN')));
            
            for m=1:height(Fields2)
                
                Fields3 = fieldnames(eval(append('CombInfo',erase(ILDataNames{l},'ILLib'),'.EdLN.',Fields2{m})));
                
                for n=1:height(Fields3)
                    
                    Actual = append('CombInfo',erase(ILDataNames{l},'ILLib'),'.EdLN',Fields2{m},Fields3{n},'0',Comp_List{j},' = ','CombInfo',erase(ILDataNames{l},'ILLib'),'.EdLN.',Fields2{m},'.',Fields3{n},'/','CombInfo',erase(ILDataNames{l},'ILLib'),'.',Comp_List{j},'.Total;');
                    eval(Actual);
                    
                end
                
            end
                        
        end
        
            end
            
        end
        
    end
    end
    % Update progress bar
    user = memory;
    RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
    LenPrint = VBUpProgBar(st,RamUsed(end),i,LenPrint);
end
end
%NameAnala = {'Monthly','Yearly'};
%NameAnala = {'Daily','Weekly','Monthly','Yearly'};
NameAnala = (Fields3(~contains(Fields3,'Fit')))';
%NameAnala = {'Monthly'};
warning('off','MATLAB:table:RowsAddedExistingVars');

% Start Progress Bar
u = StartProgBar(width(NameAnala), 1, 3-FasterLoop, 3-FasterLoop); tic; st = now;

[NumInfCases, ILData] = findILNamesStr(CombInfo);
CodesName = fieldnames(eval(append(ILData(1).Name,';'))); CodesName = CodesName(contains(CodesName,'0'));CodesName = cellfun(@(x) strsplit(x, '0'), CodesName, 'UniformOutput', false);CodesName = vertcat(CodesName{:});CodesName = CodesName(:,2);CodesName = unique(CodesName);
ClassName = fieldnames(eval(append(ILData(1).Name,'.EdLN',';')));

for b = 1:width(NameAnala)
    
    for i = 1:height(CodesName)
        % Save new format for plotting all results with AlphaSummaryPlot
        VBResults.(CodesName{i}).(NameAnala{b}) = table('Size',[1,15],'VariableTypes',["string","string","string","string","string","string","string","string","double","double","double","double","string","string","string"]);%,"double"]);
        VBResults.(CodesName{i}).(NameAnala{b}).Properties.VariableNames = {'Type','SubType','Width','Layout','Support','Trans','AE','Traffic','Span','All','ClassOW','Class','BestFitAll','BestFitClassOW','BestFitClass'};%,'PropTrucks'};
        
        for j = 1:NumInfCases
            
            try
                if AlphaAnalys ~= 1
                    Q1 = eval(append(ILData(j).Name,'.',CodesName{i},'.EQ(1)'));
                    Q2 = eval(append(ILData(j).Name,'.',CodesName{i},'.EQ(2)'));
                    Qq = sum(eval(append(ILData(j).Name,'.',CodesName{i},'.Eq')));
                end
                
                for k = 1:height(ClassName)
                    
                    Ed = eval(append(ILData(j).Name,'.EdLN.',ClassName{k},'.',NameAnala{b}));
                    if AlphaAnalys == 1
                        VBResults.(CodesName{i}).(NameAnala{b}).(ClassName{k})(j) = eval(append(ILData(j).Name,'.EdLN',ClassName{k},NameAnala{b},'0',CodesName{i}));
                        VBResults.(CodesName{i}).(NameAnala{b}).(append('BestFit',ClassName{k}))(j) = eval(append(ILData(j).Name,'.EdLN.',ClassName{k},'.',NameAnala{b},'BestFit'));
                        %VBResults.(CodesName{i}).(NameAnala{b}).(fields10{r})(a) = Ed./(LambdaQ.*(AlphaQ1.*Q1+Alphaq.*Qq+AlphaQ2.*Q2));
                    elseif AlphaAnalys == 2
                        VBResults.(CodesName{i}).(NameAnala{b}).(ClassName{k})(j) = (Ed./LambdaQ-(AlphaQ2.*Q2+Alphaq.*Qq))./Q1;
                        VBResults.(CodesName{i}).(NameAnala{b}).(append('BestFit',ClassName{k}))(j) = eval(append(ILData(j).Name,'.EdLN.',ClassName{k},'.',NameAnala{b},'BestFit'));
                    elseif AlphaAnalys == 3
                        VBResults.(CodesName{i}).(NameAnala{b}).(ClassName{k})(j) = (Ed./LambdaQ-(AlphaQ1.*Q1+Alphaq.*Qq))./Q2;
                        VBResults.(CodesName{i}).(NameAnala{b}).(append('BestFit',ClassName{k}))(j) = eval(append(ILData(j).Name,'.EdLN.',ClassName{k},'.',NameAnala{b},'BestFit'));
                    elseif AlphaAnalys == 4
                        VBResults.(CodesName{i}).(NameAnala{b}).(ClassName{k})(j) = (Ed./LambdaQ-(AlphaQ1.*Q1+AlphaQ2.*Q2))./Qq;
                        VBResults.(CodesName{i}).(NameAnala{b}).(append('BestFit',ClassName{k}))(j) = eval(append(ILData(j).Name,'.EdLN.',ClassName{k},'.',NameAnala{b},'BestFit'));
                    else
                    end
                end
            catch
                for k=1:height(ClassName)
                    VBResults.(CodesName{i}).(NameAnala{b}).(ClassName{k})(j) = 0;
                end
            end
            InflName = strsplit(ILData(j).Name,'.'); InflName = InflName(2:end);
            VBResults.(CodesName{i}).(NameAnala{b}).Type(j) = InflName{1};
            VBResults.(CodesName{i}).(NameAnala{b}).SubType(j) = InflName{2};
            VBResults.(CodesName{i}).(NameAnala{b}).Width(j) = InflName{3};
            VBResults.(CodesName{i}).(NameAnala{b}).Layout(j) = InflName{4};
            VBResults.(CodesName{i}).(NameAnala{b}).Support(j) = InflName{5};
            VBResults.(CodesName{i}).(NameAnala{b}).Trans(j) = InflName{6};
            VBResults.(CodesName{i}).(NameAnala{b}).AE(j) = InflName{7};
            if (contains(InflName{4},'Uni')||contains(InflName{4},'Bi'))&&(contains(InflName{3},'12')||contains(InflName{3},'9'))
                VBResults.(CodesName{i}).(NameAnala{b}).Traffic(j) = append(InflName{4},'2L');
            elseif ((contains(InflName{4},'Uni')||contains(InflName{4},'Bi'))&&contains(InflName{3},'15'))||(contains(InflName{4},'PUN')&&(contains(InflName{3},'12')||contains(InflName{3},'9')))
                VBResults.(CodesName{i}).(NameAnala{b}).Traffic(j) = append(InflName{4},'3L');
            elseif ((contains(InflName{4},'Uni')||contains(InflName{4},'Bi'))&&contains(InflName{3},'18'))||(contains(InflName{4},'PUN')&&contains(InflName{3},'15'))||(contains(InflName{4},'Chan')&&contains(InflName{3},'12'))
                VBResults.(CodesName{i}).(NameAnala{b}).Traffic(j) = append(InflName{4},'4L');
            elseif (contains(InflName{4},'PUN')&&contains(InflName{3},'18'))||(contains(InflName{4},'Chan')&&contains(InflName{3},'15'))
                VBResults.(CodesName{i}).(NameAnala{b}).Traffic(j) = append(InflName{4},'5L');
            elseif contains(InflName{4},'Chan')&&contains(InflName{3},'18')
                VBResults.(CodesName{i}).(NameAnala{b}).Traffic(j) = append(InflName{4},'6L');
            elseif contains(InflName{1},'Multi')
                VBResults.(CodesName{i}).(NameAnala{b}).Traffic(j) = append(InflName{4},'XL');
            end
            VBResults.(CodesName{i}).(NameAnala{b}).Span(j) = str2double(erase(InflName{8},'S'));
        end
    end
    
    VBResults.AQ = CombInfo;
    
    % Update progress bar
    user = memory;
    RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
    LenPrint = VBUpProgBar(st,RamUsed(end),b,LenPrint);
end

save(NameFileSave,'VBResults');


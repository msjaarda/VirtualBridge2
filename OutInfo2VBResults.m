clear all, clc

%% Script to create VBResults table from OutInfo results

disp("---------------------------------------------------------------------");
disp("Bonjour bienvenue dans le script OutInfo2VBResults");
disp("Vous utilisez la version 1.0 du 16/06/2022, Lucas");
disp("---------------------------------------------------------------------");

%% Open folder who contains WIM results
% First place is for the WIM or SIM folder,
% Second is for the models folder (if needed)
Dir_List = dir('Output'); % Folder containing the Output results
File_List = {Dir_List.name}';

boucle = 0; LenghtClear = 70;
while boucle == 0
disp("---------------------------------------------------------------------");
Prompt = 'Veuillez entrer le nom du fichier contenant les résultats OutInfo : \n';
Folder_Names{1} = input(Prompt,"s");
Folder_Names{2} = Folder_Names{1};
fprintf(repmat('\b',1,LenghtClear+length(Prompt)+length(Folder_Names{1})));
% check if WIM folder exist
if sum(strcmp(File_List,Folder_Names{1}))>=1
boucle = 1;
else
disp("------------------NO FOLDER IN OUTPUT TRY AGAIN----------------------");
LenghtClear = 140;
end
end

boucle = 0; LenghtClear = 70; FasterLoop = 0;
while boucle == 0
disp("---------------------------------------------------------------------");
Prompt = 'Voulez-vous importer les ECodes d''un autre fichier? 1)Oui 2)Non \n';
Entry = input(Prompt);
fprintf(repmat('\b',1,LenghtClear+length(Prompt)+numel(num2str(Entry))));
if Entry == 1
while boucle == 0
    disp("---------------------------------------------------------------------");
    Prompt = 'Veuillez entrer le nom du fichier contenant les ECodes : \n';
    Folder_Names{2} = input(Prompt,"s");
    fprintf(repmat('\b',1,LenghtClear+length(Prompt)+length(Folder_Names{2})));
        if sum(strcmp(File_List,Folder_Names{2}))>=1
        boucle = 1;
        else
        disp("------------------NO FOLDER IN OUTPUT TRY AGAIN----------------------");
        LenghtClear = 140;
        end
end
elseif Entry == 2  
FasterLoop = 1;
boucle = 1;
else
end
end

boucle = 0; LenghtClear = 70; AlphaAnalys = 0;
while boucle == 0
    disp("---------------------------------------------------------------------");
    Prompt = 'Quel alpha étudier? 1)Blended Alpha 2)AlphaQ1 3)AlphaQ2 4)Alphaq \n';
    Entry = input(Prompt);
    fprintf(repmat('\b',1,LenghtClear+length(Prompt)+numel(num2str(Entry))));
    if Entry == 1
        AlphaAnalys = 1;
        boucle = 1;
        NameFileSave = append('VBResults.mat');
    elseif Entry == 2 %AlphaQ1
        AlphaAnalys = 2;
        while boucle == 0
            disp("---------------------------------------------------------------------");
            Prompt = 'Entrez la valeur pour AlphaQ2 : \n';
            AlphaQ2 = input(Prompt);
            fprintf(repmat('\b',1,LenghtClear+length(Prompt)+numel(num2str(AlphaQ2))));
            if (AlphaQ2>=0)
                boucle = 1;
            else
                disp("----------------------------WRONG ALPHA------------------------------");
                LenghtClear = 140;
            end
        end
        boucle = 0;
        while boucle == 0
            disp("---------------------------------------------------------------------");
            Prompt = 'Entrez la valeur pour Alphaq : \n';
            Alphaq = input(Prompt);
            fprintf(repmat('\b',1,LenghtClear+length(Prompt)+numel(num2str(Alphaq))));
            if (Alphaq>=0)
                boucle = 1;
            else
                disp("----------------------------WRONG ALPHA------------------------------");
                LenghtClear = 140;
            end
        end
        NameFileSave = append('VBResultsAlphaQ1(Q2=',num2str(AlphaQ2),'etq=',num2str(Alphaq),').mat');
    elseif Entry == 3 %AlphaQ2
        AlphaAnalys = 3;
        while boucle == 0
            disp("---------------------------------------------------------------------");
            Prompt = 'Entrez la valeur pour AlphaQ1 : \n';
            AlphaQ1 = input(Prompt);
            fprintf(repmat('\b',1,LenghtClear+length(Prompt)+numel(num2str(AlphaQ1))));
            if (AlphaQ1>=0)
                boucle = 1;
            else
                disp("----------------------------WRONG ALPHA------------------------------");
                LenghtClear = 140;
            end
        end
        boucle = 0;
        while boucle == 0
            disp("---------------------------------------------------------------------");
            Prompt = 'Entrez la valeur pour Alphaq : \n';
            Alphaq = input(Prompt);
            fprintf(repmat('\b',1,LenghtClear+length(Prompt)+numel(num2str(Alphaq))));
            if (Alphaq>=0)
                boucle = 1;
            else
                disp("----------------------------WRONG ALPHA------------------------------");
                LenghtClear = 140;
            end
        end
        NameFileSave = append('VBResultsAlphaQ2(Q1=',num2str(AlphaQ1),'etq=',num2str(Alphaq),').mat');
    elseif Entry == 4 %Alphaq
        AlphaAnalys = 4;
        while boucle == 0
            disp("---------------------------------------------------------------------");
            Prompt = 'Entrez la valeur pour AlphaQ1 : \n';
            AlphaQ1 = input(Prompt);
            fprintf(repmat('\b',1,LenghtClear+length(Prompt)+numel(num2str(AlphaQ1))));
            if (AlphaQ1>=0)
                boucle = 1;
            else
                disp("----------------------------WRONG ALPHA------------------------------");
                LenghtClear = 140;
            end
        end
        boucle = 0;
        while boucle == 0
            disp("---------------------------------------------------------------------");
            Prompt = 'Entrez la valeur pour AlphaQ2 : \n';
            AlphaQ2 = input(Prompt);
            fprintf(repmat('\b',1,LenghtClear+length(Prompt)+numel(num2str(AlphaQ2))));
            if (AlphaQ2>=0)
                boucle = 1;
            else
                disp("----------------------------WRONG ALPHA------------------------------");
                LenghtClear = 140;
            end
        end
        NameFileSave = append('VBResultsAlphaq(Q1=',num2str(AlphaQ1),'etQ2=',num2str(AlphaQ2),').mat');
    else
        disp("---------------------------WRONG NUMBER------------------------------");
        LenghtClear = 140;
    end
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

for b = 1:width(NameAnala)

% Save new format for plotting all results with AlphaSummaryPlot

fields1 = fieldnames(CombInfo);
a = 1; %initialise a
for i=1:height(fields1)
    fields2 = fieldnames(CombInfo.(fields1{i}));
    for j=1:height(fields2)
        fields3 = fieldnames(CombInfo.(fields1{i}).(fields2{j}));
        for k=1:height(fields3)
            fields4 = fieldnames(CombInfo.(fields1{i}).(fields2{j}).(fields3{k}));
            for l=1:height(fields4)
                fields5 = fieldnames(CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}));
                for m=1:height(fields5)
                    fields6 = fieldnames(CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}).(fields5{m}));
                    for n=1:height(fields6)
                        fields7 = fieldnames(CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}).(fields5{m}).(fields6{n}));
                        for o=1:height(fields7)
                            fields8 = fieldnames(CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}).(fields5{m}).(fields6{n}).(fields7{o}));
                            for p=1:height(fields8)
                                fields9 = fieldnames(CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}).(fields5{m}).(fields6{n}).(fields7{o}).(fields8{p})); fields9 = fields9(contains(fields9,'0'));fields9 = cellfun(@(x) strsplit(x, '0'), fields9, 'UniformOutput', false);fields9 = vertcat(fields9{:});fields9 = fields9(:,2);fields9 = unique(fields9);
                                fields10 = fieldnames(CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}).(fields5{m}).(fields6{n}).(fields7{o}).(fields8{p}).EdLN);
                                for q=1:height(fields9)
                                    if a==1 % create table at the begining
                                        VBResults.(fields9{q}).(NameAnala{b}) = table('Size',[1,15],'VariableTypes',["string","string","string","string","string","string","string","string","double","double","double","double","string","string","string"]);%,"double"]);
                                        VBResults.(fields9{q}).(NameAnala{b}).Properties.VariableNames = {'Type','SubType','Width','Layout','Support','Trans','AE','Traffic','Span','All','ClassOW','Class','BestFitAll','BestFitClassOW','BestFitClass'};%,'PropTrucks'};
                                    end
                                    try
                                        if AlphaAnalys ~= 1
                                    Q1 = CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}).(fields5{m}).(fields6{n}).(fields7{o}).(fields8{p}).(fields9{q}).EQ(1);
                                    Q2 = CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}).(fields5{m}).(fields6{n}).(fields7{o}).(fields8{p}).(fields9{q}).EQ(2);
                                    Qq = sum(CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}).(fields5{m}).(fields6{n}).(fields7{o}).(fields8{p}).(fields9{q}).Eq);
                                        end
                                    for r=1:height(fields10)
                                    Ed = CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}).(fields5{m}).(fields6{n}).(fields7{o}).(fields8{p}).EdLN.(fields10{r}).(NameAnala{b});
                                    if AlphaAnalys == 1
                                        VBResults.(fields9{q}).(NameAnala{b}).(fields10{r})(a) = CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}).(fields5{m}).(fields6{n}).(fields7{o}).(fields8{p}).(append('EdLN',fields10{r},NameAnala{b},'0',fields9{q}));
                                        VBResults.(fields9{q}).(NameAnala{b}).(append('BestFit',fields10{r}))(a) = CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}).(fields5{m}).(fields6{n}).(fields7{o}).(fields8{p}).('EdLN').(fields10{r}).(append(NameAnala{b},'BestFit'));
                                        %VBResults.(fields9{q}).(NameAnala{b}).(fields10{r})(a) = Ed./(LambdaQ.*(AlphaQ1.*Q1+Alphaq.*Qq+AlphaQ2.*Q2));
                                    elseif AlphaAnalys == 2
                                        VBResults.(fields9{q}).(NameAnala{b}).(fields10{r})(a) = (Ed./LambdaQ-(AlphaQ2.*Q2+Alphaq.*Qq))./Q1;
                                        VBResults.(fields9{q}).(NameAnala{b}).(append('BestFit',fields10{r}))(a) = CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}).(fields5{m}).(fields6{n}).(fields7{o}).(fields8{p}).('EdLN').(fields10{r}).(append(NameAnala{b},'BestFit'));
                                    elseif AlphaAnalys == 3
                                        VBResults.(fields9{q}).(NameAnala{b}).(fields10{r})(a) = (Ed./LambdaQ-(AlphaQ1.*Q1+Alphaq.*Qq))./Q2;
                                        VBResults.(fields9{q}).(NameAnala{b}).(append('BestFit',fields10{r}))(a) = CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}).(fields5{m}).(fields6{n}).(fields7{o}).(fields8{p}).('EdLN').(fields10{r}).(append(NameAnala{b},'BestFit'));
                                    elseif AlphaAnalys == 4
                                        VBResults.(fields9{q}).(NameAnala{b}).(fields10{r})(a) = (Ed./LambdaQ-(AlphaQ1.*Q1+AlphaQ2.*Q2))./Qq;
                                        VBResults.(fields9{q}).(NameAnala{b}).(append('BestFit',fields10{r}))(a) = CombInfo.(fields1{i}).(fields2{j}).(fields3{k}).(fields4{l}).(fields5{m}).(fields6{n}).(fields7{o}).(fields8{p}).('EdLN').(fields10{r}).(append(NameAnala{b},'BestFit'));
                                    else
                                    end                           
                                    end
                                    catch
                                       for r=1:height(fields10) 
                                       VBResults.(fields9{q}).(NameAnala{b}).(fields10{r})(a) = 0;
                                       end
                                    end
                                    VBResults.(fields9{q}).(NameAnala{b}).Type(a) = fields1{i};
                                    VBResults.(fields9{q}).(NameAnala{b}).SubType(a) = fields2{j};
                                    VBResults.(fields9{q}).(NameAnala{b}).Width(a) = fields3{k};
                                    VBResults.(fields9{q}).(NameAnala{b}).Layout(a) = fields4{l};
                                    VBResults.(fields9{q}).(NameAnala{b}).Support(a) = fields5{m};
                                    VBResults.(fields9{q}).(NameAnala{b}).Trans(a) = fields6{n};
                                    VBResults.(fields9{q}).(NameAnala{b}).AE(a) = fields7{o};
                                    if (contains(fields4{l},'Uni')||contains(fields4{l},'Bi'))&&(contains(fields3{k},'12')||contains(fields3{k},'9'))
                                        VBResults.(fields9{q}).(NameAnala{b}).Traffic(a) = append(fields4{l},'2L');
                                    elseif ((contains(fields4{l},'Uni')||contains(fields4{l},'Bi'))&&contains(fields3{k},'15'))||(contains(fields4{l},'PUN')&&(contains(fields3{k},'12')||contains(fields3{k},'9')))
                                        VBResults.(fields9{q}).(NameAnala{b}).Traffic(a) = append(fields4{l},'3L');
                                    elseif ((contains(fields4{l},'Uni')||contains(fields4{l},'Bi'))&&contains(fields3{k},'18'))||(contains(fields4{l},'PUN')&&contains(fields3{k},'15'))||(contains(fields4{l},'Chan')&&contains(fields3{k},'12'))
                                        VBResults.(fields9{q}).(NameAnala{b}).Traffic(a) = append(fields4{l},'4L');
                                    elseif (contains(fields4{l},'PUN')&&contains(fields3{k},'18'))||(contains(fields4{l},'Chan')&&contains(fields3{k},'15'))
                                        VBResults.(fields9{q}).(NameAnala{b}).Traffic(a) = append(fields4{l},'5L');
                                    elseif contains(fields4{l},'Chan')&&contains(fields3{k},'18')
                                        VBResults.(fields9{q}).(NameAnala{b}).Traffic(a) = append(fields4{l},'6L');
                                    elseif contains(fields1{i},'Multi')
                                        VBResults.(fields9{q}).(NameAnala{b}).Traffic(a) = append(fields4{l},'XL');
                                    end
                                    VBResults.(fields9{q}).(NameAnala{b}).Span(a) = str2double(erase(fields8{p},'S'));
                                    %clc
                                end
                                a = a+1;
                            end
                        end
                    end
                end
            end
        end
    end      
end

VBResults.AQ = CombInfo;

% Update progress bar
    user = memory;
    RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
    LenPrint = VBUpProgBar(st,RamUsed(end),b,LenPrint);
end

save(NameFileSave,'VBResults');


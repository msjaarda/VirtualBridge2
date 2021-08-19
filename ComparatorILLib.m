clc
clear all
close all

load("ILLib.mat");

disp("-----------------------------------------------------");
disp("Bonjour bienvenue dans le script COMPARATOR");
disp("Vous utilisez la version 2.0 du 18/08/2021, Lucas");
disp("-----------------------------------------------------");
a = 1;

%% Inputs

while a ~= 0

% Type de pont
disp('-----------------------------');
Prompt = 'Type de pont : 1)Box 2)Twin 3)Multi 4)Slab \n';
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+31));
Donnee.BridgeType = {'Box';'Twin';'Multi';'Slab'};
Donnee.BridgeType = char(Donnee.BridgeType(Entry));
Fields = fieldnames(ILLib);
Fields = Fields(contains(Fields,Donnee.BridgeType));

% Si besoin catch le type de donnée
if height(Fields) > 1
    Prompt = 'Type de données : ';
    Donnee.TypDonn = erase(Fields,Donnee.BridgeType);
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Donnee.TypDonn(i)),' ');
    end
    Prompt = append(Prompt,'\n');
    Entry = input(Prompt);
    fprintf(repmat('\b',1,length(Prompt)+1));
    Donnee.TypDonn = char(Donnee.TypDonn(Entry));
    Fields = Fields(Entry);
end

Path = append('ILLib.',char(Fields));

% Sous type du pont
disp(append('Disposition ',int2str(a),' : ',Path));
disp('-----------------------------');
Prompt = 'Sous type du pont : ';
Fields = fieldnames(eval(Path));
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Fields(i)),' ');
    end
Prompt = append(Prompt,'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+1+length(Path)+47));
Fields = Fields(Entry);
Donnee.SSType = char(Fields);
    
Path = append(Path,'.',char(Fields));

% Largeur du pont
disp(append('Disposition ',int2str(a),' : ',Path));
disp('-----------------------------');
Prompt = 'Largeur du Pont : ';
Fields = fieldnames(eval(Path));
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Fields(i)),' ');
    end
Prompt = append(Prompt,'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+1+length(Path)+47));
Fields = Fields(Entry);
Donnee.Largeur = char(Fields);
    
Path = append(Path,'.',char(Fields));

% Disposition du trafic
disp(append('Disposition ',int2str(a),' : ',Path));
disp('-----------------------------');
Prompt = 'Disposition du trafic : ';
Fields = fieldnames(eval(Path));
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Fields(i)),' ');
    end
Prompt = append(Prompt,'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+1+length(Path)+47));
Fields = Fields(Entry);
Donnee.DispoTrafic = char(Fields);
    
Path = append(Path,'.',char(Fields));

% Conditions d'appuis
disp(append('Disposition ',int2str(a),' : ',Path));
disp('-----------------------------');
Prompt = 'Conditions d''appuis : ';
Fields = fieldnames(eval(Path));
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Fields(i)),' ');
    end
Prompt = append(Prompt,'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+1+length(Path)+47));
Fields = Fields(Entry);
Donnee.Support = char(Fields);
    
Path = append(Path,'.',char(Fields));

% Positition du point d'étude
disp(append('Disposition ',int2str(a),' : ',Path));
disp('-----------------------------');
Prompt = 'Position d''étude : ';
Fields = fieldnames(eval(Path));
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Fields(i)),' ');
    end
Prompt = append(Prompt,'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+1+length(Path)+47));
Fields = Fields(Entry);
Donnee.TransPos = char(Fields);
    
Path = append(Path,'.',char(Fields));

% Effort étudié
disp(append('Disposition ',int2str(a),' : ',Path));
disp('-----------------------------');
Prompt = 'Effort étudié : ';
Fields = fieldnames(eval(Path));
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Fields(i)),' ');
    end
Prompt = append(Prompt,'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+1+length(Path)+47));
Fields = Fields(Entry);
Donnee.AEstress = char(Fields);
    
Path = append(Path,'.',char(Fields));

% Portée
disp(append('Disposition ',int2str(a),' : ',Path));
disp('-----------------------------');
Prompt = 'Portée : ';
Fields = fieldnames(eval(Path));
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Fields(i)),' ');
    end
Prompt = append(Prompt,'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+1+length(Path)+47));
Fields = Fields(Entry);
Donnee.Span = char(Fields);
    
Path = append(Path,'.',char(Fields));

Graph = eval(Path);

X = Graph(:,1);
Y = Graph(:,2:end);

%% Plot des lignes d'influences

figure(1);
xlabel('Longueur du pont [m]');
ylabel('Poids d influence');
TracGraph = {'-';'--';':';'-.';'*';'s';'d'};

for i=1:width(Y)
plot(X,Y(:,i),char(TracGraph(a)));
hold on;
Legende{i} = append('Voie',int2str(i));
end

xlim([0 X(end)]);
grid on;

position = 'northeast';
if strcmp(Donnee.AEstress,'Mn')
position = 'southeast';
end

legend(Legende,'Location',position);

title(append(Donnee.BridgeType,' : portée ',num2str(X(end)),' mètres, Effort : ',Donnee.AEstress,', ',Donnee.TransPos,' ',Donnee.Support));

a = a+1;

disp(append('Disposition ',int2str(a-1),' : ',Path,' ','''',char(TracGraph(a-1)),''''));
% Continue ?
Prompt = append('Voulez vous ajouter un graphique : 1)Oui 2)Non ',', Nombre de graphique(s) actuellement : ',int2str(a-1),'\n');
Entry = input(Prompt);

if Entry == 2
   a = 0;
   fprintf(repmat('\b',1,91));
else
fprintf(repmat('\b',1,91));
end

end

clc
clear all
close all

load("ILLib.mat");

disp("-----------------------------------------------------");
disp("Bonjour bienvenue dans le script COMPARATOR");
disp("Vous utilisez la version 2.5 du 06/03/2023, Lucas");
disp("-----------------------------------------------------");
a = 1;

%% Inputs

while a ~= 0

% Type de pont
b = 0;
while b == 0
disp("-----------------------------------------------------");
Prompt = 'Type de pont : 1)Box 2)Twin 3)Multi 4)Slab 5)Dalle de roulement\n';
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+55));
if Entry > 0 && Entry < 6
b = 1;
end
end
Donnee.BridgeType = {'Box';'Twin';'Multi';'Slab';'DalleRoulem'};
Donnee.BridgeType = char(Donnee.BridgeType(Entry));
Fields = fieldnames(ILLib);
Fields = Fields(contains(Fields,Donnee.BridgeType));
Fields(contains(Fields,'Test')) = [];
Fields = flip(Fields);

b = 0;
while b == 0
disp("-----------------------------------------------------");
Prompt = 'Type de données : 1)New 2)AGB \n';
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+55));
if Entry > 0 && Entry < 3
b = 1;
end
end
Donnee.TypDonn = {'';'AGB'};
Donnee.TypDonn = char(Donnee.TypDonn(Entry));
Fields = Fields(Entry);


Path = append('ILLib.',char(Fields));

% Sous type du pont
if isempty(Donnee.TypDonn) || (strcmp(Donnee.BridgeType,'Twin') && strcmp(Donnee.TypDonn,'AGB'))
b = 0;
while b == 0
disp(append('Disposition ',int2str(a),' : ',Path));
disp("-----------------------------------------------------");
Prompt = 'Sous type du pont : ';
Fields = fieldnames(eval(Path));
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Fields(i)),' ');
    end
Prompt = append(Prompt,'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+1+length(Path)+71));
if Entry > 0 && Entry < height(Fields)+1
b = 1;
end
end
Fields = Fields(Entry);
Donnee.SSType = char(Fields);
    
Path = append(Path,'.',char(Fields));
end

% Largeur du pont
if isempty(Donnee.TypDonn)
b = 0;
while b == 0
disp(append('Disposition ',int2str(a),' : ',Path));
disp("-----------------------------------------------------");
Prompt = 'Largeur du Pont : ';
Fields = fieldnames(eval(Path));
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Fields(i)),' ');
    end
Prompt = append(Prompt,'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+1+length(Path)+71));
if Entry > 0 && Entry < height(Fields)+1
b = 1;
end
end
Fields = Fields(Entry);
Donnee.Largeur = char(Fields);
    
Path = append(Path,'.',char(Fields));
end

% Disposition du trafic
if isempty(Donnee.TypDonn)
b = 0;
while b == 0
disp(append('Disposition ',int2str(a),' : ',Path));
disp("-----------------------------------------------------");
Prompt = 'Disposition du trafic : ';
Fields = fieldnames(eval(Path));
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Fields(i)),' ');
    end
Prompt = append(Prompt,'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+1+length(Path)+71));
if Entry > 0 && Entry < height(Fields)+1
b = 1;
end
end
Fields = Fields(Entry);
Donnee.DispoTrafic = char(Fields);
    
Path = append(Path,'.',char(Fields));
end

% Conditions d'appuis
if isempty(Donnee.TypDonn) || (strcmp(Donnee.BridgeType,'Slab') && strcmp(Donnee.TypDonn,'AGB'))
b = 0;
while b == 0
disp(append('Disposition ',int2str(a),' : ',Path));
disp("-----------------------------------------------------");
Prompt = 'Conditions d''appuis : ';
Fields = fieldnames(eval(Path));
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Fields(i)),' ');
    end
Prompt = append(Prompt,'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+1+length(Path)+71));
if Entry > 0 && Entry < height(Fields)+1
b = 1;
end
end
Fields = Fields(Entry);
Donnee.Support = char(Fields);
    
Path = append(Path,'.',char(Fields));
end

% Positition du point d'étude
if isempty(Donnee.TypDonn) || (strcmp(Donnee.BridgeType,'Slab') && strcmp(Donnee.TypDonn,'AGB')) || (strcmp(Donnee.BridgeType,'Multi') && strcmp(Donnee.TypDonn,'AGB'))
b = 0;
while b == 0
disp(append('Disposition ',int2str(a),' : ',Path));
disp("-----------------------------------------------------");
Prompt = 'Position d''étude : ';
Fields = fieldnames(eval(Path));
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Fields(i)),' ');
    end
Prompt = append(Prompt,'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+1+length(Path)+71));
if Entry > 0 && Entry < height(Fields)+1
b = 1;
end
end
Fields = Fields(Entry);
Donnee.TransPos = char(Fields);
    
Path = append(Path,'.',char(Fields));
end


% Effort étudié
b = 0;
while b == 0
disp(append('Disposition ',int2str(a),' : ',Path));
disp("-----------------------------------------------------");
Prompt = 'Effort étudié : ';
Fields = fieldnames(eval(Path));
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Fields(i)),' ');
    end
Prompt = append(Prompt,'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+1+length(Path)+71));
if Entry > 0 && Entry < height(Fields)+1
b = 1;
end
end
Fields = Fields(Entry);
Donnee.AEstress = char(Fields);
    
Path = append(Path,'.',char(Fields));

% Portée
b = 0;
while b == 0
disp(append('Disposition ',int2str(a),' : ',Path));
disp("-----------------------------------------------------");
Prompt = 'Portée : ';
Fields = fieldnames(eval(Path));
    for i = 1:height(Fields)
    Prompt = append(Prompt,int2str(i),')',char(Fields(i)),' ');
    end
Prompt = append(Prompt,'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,length(Prompt)+1+length(Path)+71));
if Entry > 0 && Entry < height(Fields)+1
b = 1;
end
end
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
TracGraph = {'-';'--';':';'-.';'-';'-';'-'};
MarcGraph = {'none';'none';'none';'none';'o';'+';'.'};
CoulGraph = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840];

for i=1:width(Y)
plot(X,Y(:,i),'Color',CoulGraph(i,:),'LineStyle',char(TracGraph(a)),'Marker',char(MarcGraph(a)),'LineWidth',1.5);
hold on;
Legende{i} = append('Voie',int2str(i));
end

xlim([0 X(end)]);
grid on;

if strcmp(Donnee.BridgeType,'DalleRoulem')
    xlim([18 32]);
end

position = 'northeast';
if strcmp(Donnee.AEstress,'Mn')
position = 'southeast';
end

legend(Legende,'Location',position);

% title(append(Donnee.BridgeType,' : portée ',num2str(X(end)),' mètres, Effort : ',Donnee.AEstress,', ',Donnee.TransPos,' ',Donnee.Support));

a = a+1;

disp(append('Disposition ',int2str(a-1),' : ',Path,' ','''',char(TracGraph(a-1)),''''));
% Continue ?
disp("-----------------------------------------------------");
Prompt = append('Voulez vous ajouter un graphique : 1)Oui 2)Non 3)Supprimer ',', Nombre de graphique(s) actuellement : ',int2str(a-1),'\n');
Entry = input(Prompt);
fprintf(repmat('\b',1,103));

if Entry == 2
   a = 0;
elseif Entry == 3
   fprintf(repmat('\b',1,54)); 
   fprintf(repmat('\b',1,20+length(char(TracGraph(a-1)))+length(Path)));
   h = findobj('Type','line');
   sizeplotrem = height(h)-sizeplot;
   h = h(1:sizeplotrem);
   delete(h);
   a = a-1;
   disp("-----------------------------------------------------");
   Prompt = append('Supprimé! Continuer ? 1)Oui 2)Non ','\n');
   Entry = input(Prompt);
   fprintf(repmat('\b',1,91));
   if Entry == 2
   a = 0;
   disp("-----------------------------------------------------");
   end
else
fprintf(repmat('\b',1,54)); 
end
sizeplot = height(findobj('Type','line'));


end

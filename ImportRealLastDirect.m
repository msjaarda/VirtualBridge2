% Etude des données de trafic simulées par le RealLast
clear, clc, close all

% Nom du fichier
PDsFileNameRealLast = 'C:\Users\lmatheve\Desktop\Travail_PostPDM_RealLast\DonneesRealLast\Fall2\mGSV\Fliessverkehr\';


% Chemin du fichier
FilePathReal = 'CH_A1_Fall2_mGSV_R01_FS'; %0_f.txt';
FilePathReal2 = 'CH_A1_Fall2_mGSV_R02_FS';
Bidirectio = true;

% Nom du site
SiteName = 2001201;

for i=1:6

    if i <=3
FilePathRealTemp = append(FilePathReal,int2str(i-1),'_f.txt');
    else
FilePathRealTemp = append(FilePathReal2,int2str(6-i),'_f.txt');
    end

opts = detectImportOptions(append(PDsFileNameRealLast,FilePathRealTemp),'Delimiter', ';','DecimalSeparator', ',');
T = readtable(append(PDsFileNameRealLast,FilePathRealTemp), opts);
T.SITE(:) = SiteName;
T.DTS = datetime(strrep(string(T.Datum) + " " + string(T.Uhrzeit), ',', '.'),'InputFormat','d/M/yyyy HH:mm:ss.SS');
T.DTS.Format = 'dd/MM/yyyy HH:mm:ss';
T.Datum = [];
T.Uhrzeit = [];
T.COUNTID(:) = 1:height(T);
T.Properties.VariableNames{'Spur_Nummer'} = 'LANE';
T.LANE(:) = i;
T.DIR = ismember(T.LANE, [1 2 3]);
T.Properties.VariableNames{'Geschwindigkeit'} = 'SPEED';
T.Properties.VariableNames{'Achsanzahl'} = 'AX';
T.Properties.VariableNames{'Gesamtgewicht'} = 'GW_TOT';
T.GW_TOT = T.GW_TOT/9.81*1000; %kg
%T.Properties.VariableNames{'Fahrzeugabstand'} = 'HEADT'; %Fahrzeugabstand c'est la distance entre les véhicules en mètres alors qu'on veut HEADT en secondes
%T.HEADT2(:) = round([0;(seconds(T.DTS(2:end)-T.DTS(1:end-1)).*T.SPEED(2:end)/3.6)-T.Fzl(1:end-1)],1); méthode pour calculer les HEADT en mètres comme déjà donné
T.HEADT(:) = [0;min(seconds(T.DTS(2:end)-T.DTS(1:end-1)),99.9)];
T.GAPT(:) = [0;T.HEADT(2:end)-T.Fzl(1:end-1)./T.SPEED(2:end)*3.6];
T.GAPT(:) = max(T.GAPT(:),0);
T.GAPT(:) = min(T.GAPT(:),99.9);
T.LENTH(:) = T.Fzl*100; % [cm] a voir si on doit rajouter la longeur de la cabine et de l'arriere du camion
T.CS(:) = NaN; %Classification Swiss10. Pas encore de méthode pour passer de type à Swiss10.
T.Achslast1(isnan(T.Achslast1(:)))=0; T.AWT01(:) = T.Achslast1(:)/9.81*1000;
T.Achslast2(isnan(T.Achslast2(:)))=0; T.AWT02(:) = T.Achslast2(:)/9.81*1000;
T.Achslast3 = str2double(string(T.Achslast3)); T.Achslast3(isnan(T.Achslast3(:)))=0; T.AWT03(:) = T.Achslast3(:)/9.81*1000;
T.Achslast4 = str2double(string(T.Achslast4)); T.Achslast4(isnan(T.Achslast4(:)))=0; T.AWT04(:) = T.Achslast4(:)/9.81*1000;
T.Achslast5 = str2double(string(T.Achslast5)); T.Achslast5(isnan(T.Achslast5(:)))=0; T.AWT05(:) = T.Achslast5(:)/9.81*1000;
T.Achslast6 = str2double(string(T.Achslast6)); T.Achslast6(isnan(T.Achslast6(:)))=0; T.AWT06(:) = T.Achslast6(:)/9.81*1000;
T.Achslast7 = str2double(string(T.Achslast7)); T.Achslast7(isnan(T.Achslast7(:)))=0; T.AWT07(:) = T.Achslast7(:)/9.81*1000;
T.Achslast8 = str2double(string(T.Achslast8)); T.Achslast8(isnan(T.Achslast8(:)))=0; T.AWT08(:) = T.Achslast8(:)/9.81*1000;
T.AWT09(:) = 0;
T.Achsabst2(isnan(T.Achsabst2(:)))=0; T.W1_2(:) = T.Achsabst2(:)*100;
T.Achsabst3 = str2double(string(T.Achsabst3)); T.Achsabst3(isnan(T.Achsabst3(:)))=0; T.W2_3(:) = T.Achsabst3(:)*100;
T.Achsabst4 = str2double(string(T.Achsabst4)); T.Achsabst4(isnan(T.Achsabst4(:)))=0; T.W3_4(:) = T.Achsabst4(:)*100;
T.Achsabst5 = str2double(string(T.Achsabst5)); T.Achsabst5(isnan(T.Achsabst5(:)))=0; T.W4_5(:) = T.Achsabst5(:)*100;
T.Achsabst6 = str2double(string(T.Achsabst6)); T.Achsabst6(isnan(T.Achsabst6(:)))=0; T.W5_6(:) = T.Achsabst6(:)*100;
T.Achsabst7 = str2double(string(T.Achsabst7)); T.Achsabst7(isnan(T.Achsabst7(:)))=0; T.W6_7(:) = T.Achsabst7(:)*100;
T.Achsabst8 = str2double(string(T.Achsabst8)); T.Achsabst8(isnan(T.Achsabst8(:)))=0; T.W7_8(:) = T.Achsabst8(:)*100;
T.W8_9(:) = 0;
T = VBClassify(T);
%T.Properties.VariableNames{'Fz_typ'} = 'TYPE'; %pas tout à fait les mêmes résultats que CLASS mais OK. A étudier peut-etre un jour...
T = removevars(T, "Fz_Typ");
T = removevars(T, "Fzl");
T = removevars(T, "Fahrzeugabstand");
T = removevars(T, "Achslast1"); T = removevars(T, "Achslast2"); T = removevars(T, "Achslast3"); T = removevars(T, "Achslast4"); T = removevars(T, "Achslast5"); T = removevars(T, "Achslast6"); T = removevars(T, "Achslast7"); T = removevars(T, "Achslast8");
T = removevars(T, "Achsabst1"); T = removevars(T, "Achsabst2"); T = removevars(T, "Achsabst3"); T = removevars(T, "Achsabst4"); T = removevars(T, "Achsabst5"); T = removevars(T, "Achsabst6"); T = removevars(T, "Achsabst7"); T = removevars(T, "Achsabst8");
T = T(:, {'SITE', 'DTS', 'COUNTID','LANE','DIR','SPEED','AX','GW_TOT','HEADT','GAPT','LENTH','CS','AWT01','AWT02','AWT03','AWT04','AWT05','AWT06','AWT07','AWT08','AWT09','W1_2','W2_3','W3_4','W4_5','W5_6','W6_7','W7_8','W8_9','TYPE','CLASS'});

if i==1
PDs = T;
else
PDs = [PDs;T];
end

end

PDs = sortrows(PDs,"DTS","ascend");
save(append('WIM/',string(SiteName)), 'PDs', 'PDsFileNameRealLast', '-v7.3');

load("Sites.mat");
load("SiteLanes.mat");

if sum(PDs.SITE(1) == Sites.SITE)==0
positemp = height(Sites)+1;
else
positemp = find(PDs.SITE(1) == Sites.SITE);
end

Sites.SName(positemp) = 'REALLAST-DIRECT';
Sites.SITE(positemp) = PDs.SITE(1);
Sites.NumLanes(positemp) = i;
Sites.HWY(positemp) = 'A';
Sites.STATE(positemp) = append(FilePathReal(1:2),'-REAL');
Sites.StartDate(positemp) = min(PDs.DTS);
Sites.EndDate(positemp) = max(PDs.DTS);
Sites.TrRateEst(positemp) = sum(PDs.GW_TOT>=6000)/(height(PDs)); % Estimation avec les véhicules plus lourds que 6 tonnes
Sites.ADTT(positemp) = sum(PDs.GW_TOT>=6000)/days(max(PDs.DTS)-min(PDs.DTS)); % On retient que les véhicules lourds...
Sites.Light(positemp) = 1; %A voir mais on pourrait optimiser la taille en supprimant les véhicules légers avant de sauvegarder
Sites.Core(positemp) = true;
Sites.COUNTRY(positemp) = append(FilePathReal(1:2),'-REAL');
if Bidirectio
    Sites.Layout(positemp) = str2double([repmat('1',1,ceil(i/2)), ...
                                repmat('2',1,floor(i/2))]);
else
    Sites.Layout(positemp) = str2double(repmat('1',1,i));
end

if sum(PDs.SITE(1) == SiteLanes.SITE)==0
positemp = height(SiteLanes)+1;
else
positemp = find(PDs.SITE(1) == SiteLanes.SITE);
positemp = positemp(1);
end

for j=1:i
SiteLanes.SName(positemp+j-1) = 'REALLAST-DIRECT';
SiteLanes.SITE(positemp+j-1) = PDs.SITE(1);
SiteLanes.LANE(positemp+j-1) = j;
SiteLanes.DIR(positemp+j-1) = 'REALLAST-DIRECT';
    if j>ceil(i/2) && Bidirectio
        SiteLanes.NSEW(positemp+j-1) = 2;
        SiteLanes.LANETYPE(positemp+j-1) = i+1-j;
    else
        SiteLanes.NSEW(positemp+j-1) = 1;
        SiteLanes.LANETYPE(positemp+j-1) = j;
    end
SiteLanes.FROM(positemp+j-1) = 'REALLAST-DIRECT';
SiteLanes.WIM(positemp+j-1) = 1;
end

save('Misc/Sites.mat', 'Sites')
save('Misc/SiteLanes.mat', 'SiteLanes')

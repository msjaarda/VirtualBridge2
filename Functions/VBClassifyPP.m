function [PD] = VBClassifyPP(PD)
%VB Classify PP (plus plus)

load('ClassPPInfo')       % CPPInfo
load('ClassPPWheelBases') % UWBr_CT

% Use KUBA-ST database to increase the number of special transports
% classified (call this new one Class++)

% Get wheelbases and put into m
PDw = PD{:,22:29}; PDw = PDw/100;
% Round to nearest 0.2
PDwr = 2*round(PDw/2,1);

% Grab corresponding matrix from KubaST
PDkr = UWBr_CT(:,1:8);

% Rename for convenience
A = PDwr;
B = PDkr;

[C,iA,iB] = intersect(A,B,'rows');
% If the 'rows' option is specified, then C = A(ia,:) and C = B(ib,:), C = A(iA,:), C = B(iB,:)

ABID = zeros(length(A),1);
for k = 1:size(C,1)
    ABID(ismember(A,C(k,:),'rows')) = iB(k);
end

% We can only add new CLASS to those that don't already have a class
% AND, I propse we make sure that the newly classified are of at least 50
% tonnes...

% Total number assigned class by procedure
ZZ = PD(ABID > 0 & PD.GW_TOT > 50000,:);
% Total number of these that didn't already have a class
YY = ZZ(ZZ.CLASS == 0,:);

% Before procedure
ClassPBe = 100*sum(PD.CLASS > 0)/height(PD);
% After
ClassPAf = 100*(sum(PD.CLASS > 0)+height(YY))/height(PD);

% Original Number of Class+
ClassPlusBe = sum(PD.CLASS > 40 & PD.CLASS < 80);
fprintf('\nTotal Before: %i\n',ClassPlusBe)
% After
ClassPlusAf = ClassPlusBe + height(YY);
fprintf('Total After:  %i (%% %0.1f)\n\n ',ClassPlusAf,100*(ClassPlusAf-ClassPlusBe)/ClassPlusBe)

% Change Class (add an offset of 50)
PD.CLASS(ABID > 0 & PD.CLASS == 0) = ABID(ABID > 0 & PD.CLASS == 0)+ 50;

% New ++ Classes Histogram
figure; histogram(PD.CLASS(PD.CLASS > 40 & PD.CLASS < 80)); ylabel('#')
xlabel('Vehicle Class')
% Weights of only classified
figure; histogram(YY.GW_TOT/1000,40,'normalization','pdf','DisplayName','Newly Classified')
set(gca,'ytick',[],'yticklabel',[],'ycolor','k'); ylabel('Normalized PDFs')
xlabel('Weight (tonnes)')
% Weights of all
hold on; histogram(PD.GW_TOT(PD.GW_TOT > 6000)/1000,30,'normalization','pdf','DisplayName','All Vehicles')

% Save new var? Not yet...


end


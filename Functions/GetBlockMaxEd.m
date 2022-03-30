function [Ed, AQ, Aq] = GetBlockMaxEd(Data,BlockM,Dist,ESIAT,ESIAEQ,ESIAEq,AQ1,AQ2,PropTruck,FitType)
%GETBLOCKMAXEd Fits, and optionally plots, BlockMaximumData
%   Data    - simply the block maximum data (max moment, shear, etc. during period)
%   BlockM  - string, 'Daily', 'Weekly', 'Monthly', 'Yearly', or 'Lifetime'
%   Dist    - string, 'Nomral, 'Lognormal'
%   PropTruck - double, proportion of special transport yearly maxima with
%   accompaniment, = 1 if dont needed
%   FitType - double, 1 : original fit method, 2 : tail fitting method, 3 :
%   Alain's fitting method (fitting each class individually)


% Data = Max.Max for Matt [old code] or Data = Max for Lucas

% --- Calculate Real Design Value Ed ---

if strcmp(BlockM,'Yearly')
    n = 1;
elseif strcmp(BlockM,'Weekly')
    n = 1/50;
elseif strcmp(BlockM,'Daily')
    n = 1/(5*50);
elseif strcmp(BlockM,'Monthly')
    n = 1/12;
elseif strcmp(BlockM,'Lifetime')
    n = 50;
else
    n = BlockM;
end
    
% Beta.Yearly = 4.700; % Here we use the Beta annual - so we should use annual max effects
% Beta.Weekly = 5.444; % See Tail Fitting > Beta Conversion
% Beta.Daily = 5.724;
% Beta.Lifetime = 3.830;


Beta = norminv(1-n*0.0000013/PropTruck);
Alpha = 0.7;

if contains(Dist,'Zero') && PropTruck ~= 0 % Consider the zero if needed 
Data(end+1:round(length(Data)/PropTruck)) = 0;
Beta = norminv(1-n*0.0000013);
Dist = erase(Dist,'Zero');
elseif contains(Dist,'Zero')
Dist = erase(Dist,'Zero');    
end

if FitType == 1 || FitType == 2 
Prop = 0.95;
if FitType == 2
if strcmp(Dist,'Normal')
    %mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear');
    mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),sort(Data),'linear','Weights',[0.1*ones(round(length(Data)*Prop),1);1*ones(length(Data)-round(length(Data)*(Prop)),1)]);
    pd = makedist('normal',mdlx.Coefficients.Estimate(1),mdlx.Coefficients.Estimate(2));
    Em = mean(Data);
    Stdev = std(Data);
    COV = Stdev/Em;
    Delta2 = log(COV^2+1);
elseif strcmp(Dist,'Lognormal')
    mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear','Weights',[0.1*ones(round(length(Data)*Prop),1);1*ones(length(Data)-round(length(Data)*(Prop)),1)]);
    %mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear');
    muu = mdlx.Coefficients.Estimate(1);
    sig = mdlx.Coefficients.Estimate(2);
    pd = makedist('lognormal',mdlx.Coefficients.Estimate(1),mdlx.Coefficients.Estimate(2));  
    Em = exp(muu+sig^2/2);
    %Stdev = sqrt(exp(2*muu+sig^2)*(exp(sig^2)-1));
    %COV = Stdev/Em;
    Delta2 = sig^2;
end
elseif FitType == 1
    Em = mean(Data);
    Stdev = std(Data);
    COV = Stdev/Em;
    Delta2 = log(COV^2+1);
end



if strcmp(Dist,'Normal')
    Ed = Em*(1+Alpha*Beta*COV);
    % FYI ESIAT has the 1.5 in it already...
    AQ = Ed/(ESIAT);
    Aq = ((Ed/1.5)-AQ1*ESIAEQ(1)-AQ2*ESIAEQ(2))/(sum(ESIAEq));
elseif strcmp(Dist,'Lognormal')
    Ed = Em*exp(Alpha*Beta*sqrt(Delta2)-0.5*Delta2);
    AQ = Ed/(ESIAT);
    Aq = ((Ed/1.5)-AQ1*ESIAEQ(1)-AQ2*ESIAEQ(2))/(sum(ESIAEq));
elseif strcmp(Dist,'Extreme Value')
    Ed = Em*(1 + COV*(0.45 + 0.78*log(-log(normpdf(Alpha*Beta)))));    %            exp(Alpha*Beta*sqrt(Delta2)-0.5*Delta2);
    AQ = Ed/(ESIAT);
    Aq = 1;   
end

elseif FitType == 3
Types = Data.m; Data = Data.Max; m1 = Types == 1; m2 = Types == 2; m3 = Types == 3;

Prop = 0.95;
for i=1:3
   DataTemp = eval(append('Data(m',int2str(i),')'));
if strcmp(Dist,'Normal')
    mdlx{i} = fitlm(norminv((1:length(DataTemp))/(length(DataTemp) + 1)),sort(DataTemp),'linear','Weights',[0.1*ones(round(length(DataTemp)*Prop),1);1*ones(length(DataTemp)-round(length(DataTemp)*(Prop)),1)]);
    pd{i} = makedist('normal',mdlx{i}.Coefficients.Estimate(1),mdlx{i}.Coefficients.Estimate(2));
    Em(i) = mean(DataTemp);
    Stdev(i) = std(DataTemp);
elseif strcmp(Dist,'Lognormal')
    mdlx{i} = fitlm(norminv((1:length(DataTemp))/(length(DataTemp) + 1)),log(sort(DataTemp)),'linear','Weights',[0.1*ones(round(length(DataTemp)*Prop),1);1*ones(length(DataTemp)-round(length(DataTemp)*(Prop)),1)]);
    muu(i) = mdlx{i}.Coefficients.Estimate(1);
    sig(i) = mdlx{i}.Coefficients.Estimate(2);
    pd{i} = makedist('lognormal',mdlx{i}.Coefficients.Estimate(1),mdlx{i}.Coefficients.Estimate(2));  
    Em(i) = exp(muu(i)+sig(i)^2/2);
    Stdev(i) = sqrt(exp(2*muu(i)+sig(i)^2)*(exp(sig(i)^2)-1));
     if mdlx{i}.Rsquared.Ordinary*100<85 %check if the Rsquared is not to bad. If yes, check if its not due to the near zero values (remove them if needed).
         if Em(i)>3
             DataTempNZer = DataTemp(DataTemp>=0.8);
             mdlxNZer = fitlm(norminv((1:length(DataTempNZer))/(length(DataTempNZer) + 1)),log(sort(DataTempNZer)),'linear','Weights',[0.1*ones(round(length(DataTempNZer)*Prop),1);1*ones(length(DataTempNZer)-round(length(DataTempNZer)*(Prop)),1)]);
             if mdlxNZer.Rsquared.Ordinary*100>mdlx{i}.Rsquared.Ordinary*100
                 DataTemp(DataTemp<=0.8) = [];
                 mdlx{i} = fitlm(norminv((1:length(DataTemp))/(length(DataTemp) + 1)),log(sort(DataTemp)),'linear','Weights',[0.1*ones(round(length(DataTemp)*Prop),1);1*ones(length(DataTemp)-round(length(DataTemp)*(Prop)),1)]);
                 muu(i) = mdlx{i}.Coefficients.Estimate(1);
                 sig(i) = mdlx{i}.Coefficients.Estimate(2);
                 pd{i} = makedist('lognormal',mdlx{i}.Coefficients.Estimate(1),mdlx{i}.Coefficients.Estimate(2));
                 Em(i) = exp(muu(i)+sig(i)^2/2);
                 Stdev(i) = sqrt(exp(2*muu(i)+sig(i)^2)*(exp(sig(i)^2)-1));
             end
         end
     end
end

COV(i) = Stdev(i)/Em(i);
Delta2(i) = log(COV(i)^2+1);

if strcmp(Dist,'Normal')
    Ed(i) = Em(i)*(1+Alpha*Beta*COV(i));
    % FYI ESIAT has the 1.5 in it already...
    AQ(i) = Ed(i)/(ESIAT);
    Aq(i) = ((Ed(i)/1.5)-AQ1*ESIAEQ(1)-AQ2*ESIAEQ(2))/(sum(ESIAEq));
elseif strcmp(Dist,'Lognormal')
    Ed(i) = Em(i)*exp(Alpha*Beta*sqrt(Delta2(i))-0.5*Delta2(i));
    AQ(i) = Ed(i)/(ESIAT);
    Aq(i) = ((Ed(i)/1.5)-AQ1*ESIAEQ(1)-AQ2*ESIAEQ(2))/(sum(ESIAEq));
elseif strcmp(Dist,'Extreme Value')
    Ed(i) = Em(i)*(1 + COV(i)*(0.45 + 0.78*log(-log(normpdf(Alpha*Beta)))));    %            exp(Alpha*Beta*sqrt(Delta2)-0.5*Delta2);
    AQ(i) = Ed(i)/(ESIAT);
    Aq(i) = 1;   
end
end
Ed(4) = Ed(1)*sum(m1)/height(m1)+Ed(2)*sum(m2)/height(m2)+Ed(3)*sum(m3)/height(m3);
    
end

end
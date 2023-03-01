function [model2,targetRID,extype] = modelSetting(model,targetMet)

target=findMetIDs(model,targetMet);
m=size(model.mets,1);
n=size(model.rxns,1);
if isempty(find(strcmp(model.rxns,strcat('EX_',targetMet))))==0
    targetRID=find(strcmp(model.rxns,strcat('EX_',targetMet)));
    model2=model;
    extype=1;
elseif isempty(find(strcmp(model.rxns,strcat('DM_',targetMet))))==0
    targetRID=find(strcmp(model.rxns,strcat('DM_',targetMet)));
    model2=model;
    extype=2;
else
    [model2,rxnIDexists]=addReaction(model,'Transport',{targetMet},[-1]);
    m=size(model2.mets,1);
    n=size(model2.rxns,1);
    model2.S(target,n)=-1;
    model2.ub(n)=999999;
    model2.lb(n)=0;
    model2.rev(n)=0;
    targetRID=n;
    extype=3;
end
%save('modelSetting.mat');
end


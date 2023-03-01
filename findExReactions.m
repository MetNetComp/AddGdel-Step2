function [ex] = findExReactions(model)

m=size(model.mets,1);
n=size(model.rxns,1);

k=1;
for i=1:n
    if size(find(model.S(:,i)),1)==1
       ex.R(k,1)=i;
       if strcmp(model.rxns{i},'Transport')~=1
           ex.R2{k,1}=model.rxns{i};
       else
           ex.R2{k,1}='Auxiliary production reaction';
       end
       ex.met(k,1)=find(model.S(:,i));
       ex.met2{k,1}=model.mets{ex.met(k,1)};
       k=k+1;
    end
end

end


function [replaceList] = findReplacableGenes(model)

g=size(model.genes,1);

[term,ng,nt,nr,nko,reactionKO,reactionKO2term] = readGeneRules(model);
t=size(term,2);

for i=1:g
    occ{i,1}=find(contains(model.grRules,model.genes{i}));
end
list0=[];
for i=1:g
    i
    for j=1:g
        if i~=j
            if isequal(occ{i,1},occ{j,1})==1
                list0(i,j)=1;
            else
                list0(i,j)=0;
            end
        end
    end
end
same=zeros(g,g);
for i=1:g
    i
    for j=1:g
        if i==j
            same(i,j)=0;
        else i~=j
            if list0(i,j)==1
                same(i,j)=1;
                for k=1:t
                    s=size(term(k).input,1);
                    flag1=0;flag2=0;
                    for p=1:s
                        tmp1=isempty(find(strcmp(term(k).input{p},model.genes{i})));
                        if tmp1==0
                            flag1=1;
                        end
                        tmp2=isempty(find(strcmp(term(k).input{p},model.genes{j})));
                        if tmp2==0
                            flag2=1;
                        end
                    end
                    if flag1~=flag2
                        same(i,j)=0;
                    end
                    
                    if flag1==1 && flag2==1 && strcmp(term(k).function,'and')~=1
                        same(i,j)=0;
                    end
                    
                end
            end
        end
    end
end
replaceList=zeros(g,1);
replaceList2=model.genes;
replaceList3=[];
for j=1:g
    j
    flag=0;
    for i=1:g
        if same(i,j)==1 && flag==0 && i<j
            replaceList(j,1)=i;
            replaceList2{j,2}=model.genes{i};
            flag=1;
        end
    end
end
for i=1:g
    tmp=find(replaceList==i);
    tmp2=(model.genes(tmp))'
    tmp3=horzcat({model.genes{i}},tmp2);
    tmp4=size(tmp3,2);
    if tmp4>0
        for j=1:tmp4
            replaceList3{i,j}=tmp3(1,j);
        end
    end
end

%save('findReplacableGenes.mat');
end


function [outputArg1,outputArg2] = example()
initCobraToolbox;
load('gDel-minRN-strategy.mat');
load('e_coli_core.mat');
model=e_coli_core;

[replaceList] = findReplacableGenes(model)
[gvalue] = addGdelStep2(model,'succ_e',gvalue,replaceList)

save('example.mat');
end


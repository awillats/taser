%% ask for input order

fprintf('total file length = %d\n',sum(sum(Ls)))
allSuperDirNames
asdn = allSuperDirNames

i=0;
fprintf('\n')
fprintf('\n')
fprintf('\n')
for fn = allSuperDirNames
    fprintf('%i: %s  ',i,fn{1})
    i = i+1;
    if (mod(i,6) == 0)
        fprintf('\n')
    end
end
fprintf('\n')
nI = 0;
order = [];

while (nI < sdN)
    o = input('next stim: ');
    order = [order,o];
    nI = length(order);
    
end
%%
fprintf('\n')
revOrder = fliplr(order)
for o = revOrder
    fprintf('%s,',allSuperDirNames{o+1})
    
end
fprintf('\n')
for o = revOrder
    fprintf('%i,',o)
end
fprintf('\n')

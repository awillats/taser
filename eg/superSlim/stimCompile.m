%% TO DO
%{
-improve ordering (display how many left?
-display compiled results
-turn these into functions??
%}

clear

D = dir;
allNames = {D.name};
isDir = [D.isdir];
allSuperDirNames = allNames(find(isDir));
allSuperDirNames = allSuperDirNames(3:end)

sdN = length(allSuperDirNames);

CWO_labels = {'C_','W_','O_','ST_', 'CT_'};
C=1; W=2; O=3; ST=4; CT=5;
nLab = length(CWO_labels);


%file lengths
Ls = zeros(sdN,nLab)-1;

%start indices
StrtIs = zeros(sdN,nLab)-1;

%Stim names
StimNames = cell(sdN,1);

%Compiled Stim Files
Compiled = cell(1,nLab);
%CompiledFileNames = cell(sdN,3);


% CA_i = 0;
% WA_i = 0;
% OA_i = 0;
% TA_i = 0;

Lprev = 0;
for sdI = 1:sdN
    sdName = allSuperDirNames{sdI};
    
    dd = dir(sdName);
    
    subNames = {dd.name};
    isDir = [dd.isdir];
    subDirs = subNames(find(isDir));
    subDirs = subDirs(3:end);
    
    for subs = subDirs
        oldD = pwd;
        cd(fullfile(sdName,subs{1}));
        
        D = dir;
        allNames = {D.name};
        isFile = ~[D.isdir];
        allFileNames = allNames(find(isFile));
        nFile = length(allFileNames);
        
        for fn = allFileNames
            fn_ = fn{1};
            
            %%classification
            for i = 1:nLab
                if (strfind(fn_, CWO_labels{i}))
                    CWO = i;
                end
            end
            try
                StrtIsCol = StrtIs(:,CWO);
                LsCol = Ls(:,CWO);
                realStarts = StrtIsCol(StrtIsCol~=-1);
                realLs = LsCol(LsCol ~= -1);
                
                x = dlmread(fn_);
                L = length(x);
                Ls(sdI,CWO) = L;
                
                if length(realStarts) > 0
                    %disp('a')
                    StrtIs(sdI,CWO) = sum(realLs);
                else
                    StrtIs(sdI,CWO) = 0;
                end
                
                Compiled{CWO} = [Compiled{CWO}; x];
            catch e
                disp(e)
                disp(strcat(fn{1},' does not have dlm data'));
            end
                
        end

        cd(oldD);
    end
    
end

mkdir('meta')

fspec = '%1.9e\n';

for i = 1:length(CWO_labels)
    
    lab = CWO_labels{i};
    fid = fopen(strcat('meta/',lab,'Lens.txt'),'w');
    fprintf(fid,fspec,Ls(:,i));
    fclose(fid);
    
    fid = fopen(strcat('meta/',lab,'StartI.txt'),'w');
    fprintf(fid,fspec,StrtIs(:,i));
    fclose(fid);
    
    fid = fopen(strcat('meta/',lab,'Comp.txt'),'w');
    fprintf(fid,fspec,Compiled{i});
    fclose(fid);
    
    fid = fopen(strcat('meta/',lab,'CompLen.txt'),'w');
    fprintf(fid,fspec,length(Compiled{i}));
    fclose(fid);
end

fid = fopen('meta/superDirs.txt','w');
fspec = '%s\n';
[nrows,ncols] = size(allSuperDirNames);

for row = 1:nrows
    fprintf(fid,fspec,allSuperDirNames{row,:});
end
fclose all;


clf
close all

idxOrder  = 1:sdN;
%idxOrder = order
labelOrder = [1,5,2,3,4]; %1:sdN
sampleFreq = 25e3;

lineSpecs = jet(nLab);


sdNames = allSuperDirNames;



for i = idxOrder
    %title(sdN
    figure('Name',allSuperDirNames{i})
    
    for k = 1:nLab
        j = labelOrder(k);
        subplot(nLab,1,k);
        stimCompRow = Compiled{j};
        
        if (StrtIs(i,j) ~= -1)
            is = StrtIs(i,j)+1;
            ie = is + Ls(i,j)-1;
            t = (0:(ie-is))/sampleFreq;
            plot(t,stimCompRow(is:ie),'color',lineSpecs(j,:))
        else
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
        end
        
        lab = CWO_labels{j};
        lab = lab(1:end-1); %cut trailing underscore
        ylabel(lab)
    end
    
end
function [] = PlotAllLevels(S)
% S = [et n L J F M et(w/o artificial degeneracy)]

[~,unique_inds,~] = unique(S(:,7));

for i1 = unique_inds'
    
    et = S(i1,1);
    n  = S(i1,2);
    L  = S(i1,3);
    J  = S(i1,4);
    
    if L == 0
        let = 'S';
    elseif L == 1
        let = 'P';
    elseif L == 2
        let = 'D';
    elseif L == 3
        let = 'F';
    end
    
    if J == 1/2
        frac = '1/2';
    elseif J == 3/2
        frac = '3/2';
    elseif J == 5/2
        frac = '5/2';
    elseif J == 7/2
        frac = '7/2';
    end
    
    line([L+J-0.5,L+0.5+J-0.5],[1 1]*et)
    text(L+J-0.5,et+3e13,[num2str(n) let '_{' frac '}'])
end
set(gcf,'Color','w')
set(gca,'XTickLabel','')
ylabel('Frequency (Hz)')
end

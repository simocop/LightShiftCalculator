function [] = PlotLevels(S,fh)
% S = [et n L J F M]

if nargin < 2
    fh = gcf;
end

clf(fh)
for i1 = 1:length(S(:,1));
    
    et = S(i1,1)/1e6;
%     n  = S(i1,2);
%     L  = S(i1,3);
%     J  = S(i1,4);
    F  = S(i1,5);
    M =  S(i1,6);
    
    x1 = M-0.25;
    x2 = M+0.25;
    line([x1,x2],[1 1]*et)
    text(M,et,num2str(F))
end
set(gcf,'Color','w')
ylabel('Frequency (MHz)')
xlabel('m_F')
end

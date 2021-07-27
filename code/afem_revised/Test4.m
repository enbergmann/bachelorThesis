function Test4()
tic;
initAFEM; 
B = Poisson_SquareExact(1);

% %% Case 1 (k, k-1, k+1)
% for k = 1:2
%     kU = k; kQ = k-1; kV = k+1;
%     S = AFEMRun_PoissonDPG(B, kU, kQ, kV,0);
%     generateOutput(S, 1, k);
% end
% 
% %% Case 2 (k-1, k-1, k)
% for k = 2:3
%     kU = k-1; kQ = k-1; kV = k;
%     S = AFEMRun_PoissonDPG(B, kU, kQ, kV,0);
%     generateOutput(S, 2, k);
% end
% %% Case 3 (k, k-1, k)
% for k = 1:3
%     kU = k; kQ = k-1; kV = k;
%     S = AFEMRun_PoissonDPG(B, kU, kQ, kV,0);
%     generateOutput(S, 3, k);
% end
%% Case 4 (k, k-3, k)
for k = 3
    kU = k; kQ = k-3; kV = k;
    S = AFEMRun_PoissonDPG(B, kU, kQ, kV,0);
    generateOutput(S, 4, k);
end
% %% Case 5 (k, k, k)
% for k = 1:2
%     kU = k; kQ = k; kV = k;
%     S = AFEMRun_PoissonDPG(B, kU, kQ, kV,0);
%     generateOutput(S, 5, k);
% end
toc;

end

function generateOutput(S, caseNr, k)
    name = ['../data/case',num2str(caseNr),'_',num2str(k),'.csv'];
    file = fopen(name,'w');
    fprintf(file, 'a,b,c,d,e,f,g,h\r\n');
    maxLvl = size(S.level,2);
    for lvl = 1:maxLvl
       fprintf(file, '%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%d\r\n',2^lvl,...
           S.level(lvl).errH1,S.level(lvl).rateH1,...
           S.level(lvl).errL2,S.level(lvl).rateL2,S.level(lvl).errEnergy,...
           S.level(lvl).eta,S.level(lvl).ndof);
    end
    fclose(file);
end


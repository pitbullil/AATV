function MAD = CalcMAD(ImTest,ImO)
% IdxIn = 93;%130
% IdxOut = 10;

load Arterymask256.mat
MaskIn = Bin;MaskOut = logical(1-Bin);
ImO = ImO-min(ImO(:)); ImO = ImO/max(ImO(:)); 
InShapO = mean(mean(ImO(MaskIn)));BGO = mean(mean(ImO(MaskOut)));
ImTest = ImTest-min(ImTest(:)); ImTest = ImTest/max(ImTest(:)); 
InShap = mean(mean(ImTest(MaskIn)));BG = mean(mean(ImTest(MaskOut)));
ImTest = ImTest*((InShapO - BGO)/(InShap - BG))- BG + BGO; 

MAD = sum(abs(ImTest(:) - ImO(:)))./(2*length(ImO));
end


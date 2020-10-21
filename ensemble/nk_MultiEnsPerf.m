function [hEPerf, hE, hSim] = nk_MultiEnsPerf(E, sE, L, C, G)
global MULTI RAND
hSim = []; ProbComp=[]; 
if isfield(MULTI,'ProbComp'), ProbComp=MULTI.ProbComp; end
switch MULTI.method
    case 1 % Simple One-Vs-One / One-vs-All
        [hE, hEPerf, ~, hSim] = nk_MultiDecideMaxWins(E, L, C, G, MULTI.decisiontype, 1);
    case 2 % Error-Correcting Output Codes
        if ~isfield(RAND,'Decompose'),
            decompose = 1;
        else
            decompose = RAND.Decompose;
        end
        [hE, hEPerf, ~, hSim] = nk_MultiDecideErrorCorrOutCodes(sE, L, C, G, decompose, MULTI.decoding, 1, ProbComp);
    case 3 % Directed Acyclic Graph
        [hE, hEPerf, hSim] = nk_MultiDecideDAG(E, L, C, G);
        
end


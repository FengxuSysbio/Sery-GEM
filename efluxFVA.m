function [solutionWT, solutionMU, fvaMin, fvaMax] = efluxFVA(model, expressionDataWT, expressionDataMU, fvaMin, fvaMax)
% Estimate flux distributions given a wild-type and mutant (or other condition) gene expression
% data using the method E-flux
% Used in this publication: https://journals.plos.org/ploscompbiol/article/authors?id=10.1371/journal.pcbi.1006692
% json outputs of the flux distributions for visualization in escher map
% will be saved automatically
%
% USAGE:
%    [solutionWT, solutionMU, fvaMin, fvaMax] = efluxFVA(model, expressionDataWT, expressionDataMU, fvaMin, fvaMax)
%
% INPUTS:
%    model                 COBRA model
%    expressionDataWT      wildtype gene expression data, see `mapExpressionToReactions.m` for more details
%    expressionDataMU      mutant gene expression data, see `mapExpressionToReactions.m` for more details
%  
% OPTIONAL INPUTS:
%    fvaMin                precalculated FVA range minimum. Perform FVA if not supplied
%    fvaMax                precalculated FVA range maximum. Perform FVA if not supplied
%
% OUTPUTS:
%    solutionWT            the estimated flux distribution for the wild type
%    solutionMU            the estimated flux distribution for the mutant
%    fvaMin                fvaMin computed if not supplied
%    fvaMax                fvaMax computed if not supplied

if nargin < 3
    % no FVA range given, calculate
    modelFVA = model;
    % relax any required fluxes
    modelFVA.lb(modelFVA.lb > 0) = 0;
    modelFVA.ub(modelFVA.ub < 0) = 0;
    % run simple FVA.
    [fvaMin, fvaMax] = fluxVariability(modelFVA, 0, 'max', modelFVA.rxns, 2, 1);
    % for those reactions with very large flux range, possibly in loops,
    % re-run using loopless FVA
    rLoops = fvaMin < -200 | fvaMax > 200;
    [fvaMin2, fvaMax2] = fluxVariability(modelFVA, 0, 'max', modelFVA.rxns(rLoops), 2, 0);
    fvaMin(rLoops) = fvaMin2;
    fvaMax(rLoops) = fvaMax2;
    fvaMin(abs(fvaMin) < 1e-8) = 0;
    fvaMax(abs(fvaMax) < 1e-8) = 0;
    if ~(all(fvaMin <= 0) && all(fvaMax >= 0))
        error('Problematic FVA results')
    end
end

%% Gene expression data

% mapExpressionToReactions function determines the expression data associated to each reaction present in
% the model following standards

% wild type
[expressionRxnsWT, parsedGPRwt] = mapExpressionToReactions(model, expressionDataWT);

% mutant (or different condition)
[expressionRxnsMU, parsedGPRmu] = mapExpressionToReactions(model, expressionDataMU);

%% E-flux using FVA bounds
% for each reaction, either both have expression values or none of them have expression values

[lbWT, lbMU] = deal(fvaMin);
[ubWT, ubMU] = deal(fvaMax);

rxnWtExpression = expressionRxnsWT > 0 & expressionRxnsMU > 0;

% for each reaction, get the condition with lower expression 
[~, ind] = min([expressionRxnsWT, expressionRxnsMU], [], 2);
% for each reaction, set the bounds for the condition with lower expression
% as bound * expression ratio

% reactions with lower expression in wildtype WT
rxnLowExprWT = rxnWtExpression & ind == 2;
lbWT(rxnLowExprWT) = minFlux(rxnLowExprWT) .* expressionRxns_C08(rxnLowExprWT) ./ expressionRxns_W05(rxnLowExprWT);
ubWT(rxnLowExprWT) = maxFlux(rxnLowExprWT) .* expressionRxns_C08(rxnLowExprWT) ./ expressionRxns_W05(rxnLowExprWT);

% reactions with lower expression in mutant MU
rxnLowExprMU = rxnWtExpression & ind == 1;
lbMU(rxnLowExprMU) = minFlux(rxnLowExprMU) .* expressionRxns_W05(rxnLowExprMU) ./ expressionRxns_C08(rxnLowExprMU);
ubMU(rxnLowExprMU) = maxFlux(rxnLowExprMU) .* expressionRxns_W05(rxnLowExprMU) ./ expressionRxns_C08(rxnLowExprMU);

[modelMU, modelWT] = deal(model);
[modelMU.lb, modelMU.ub, modelWT.lb, modelWT.ub] = deal(lbMU, ubMU, lbWT, ubWT);
[solutionWT, solutionMU, totalFluxDiff] = optimizeTwoCbModels(modelWT, modelMU, 'max', false, true);
 
% json output for visualization in escher map
[fWT, fMU] = deal(struct());
for j = 1:numel(model.rxns)
    fWT.(model.rxns{j}) = solutionWT.x(j);
    fMU.(model.rxns{j}) = solutionMU.x(j);
end
f = fopen('fluxWT.json','w'); 
fprintf(f, jsonencode(fWT));
fclose(f);
f = fopen('fluxMU.json','w'); 
fprintf(f, jsonencode(fMU));
fclose(f);

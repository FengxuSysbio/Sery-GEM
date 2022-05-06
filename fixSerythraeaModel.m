%% create the model from the excel file
if exist('Serythraea.mat', 'file')
    model = readCbModel('Serythraea.mat');
else
    [~, ~, rxnInfo]=xlsread('model_fixed.xlsx','Reaction List');
    [~, ~, metInfo]=xlsread('model_fixed.xlsx','Metabolite List');
    
    rxnInfo = rxnInfo(~cellfun(@(x) isempty(x) || any(isnan(x)), rxnInfo(:, 1)), 1:15);
    rxnNumericHeaders = {'Reversible','Lower bound','Upper bound','Objective','Confidence Score'};
    rxnNumericCols = ismember(rxnInfo(1, :), rxnNumericHeaders);
    for j = 1:size(rxnInfo, 2)
        if ~rxnNumericCols(j)
            rxnInfo(cellfun(@(x) any(isnan(x)), rxnInfo(:, j)), j) = {''};
        end
    end
    
    rxnInfo(:, 3) = strrep(rxnInfo(:, 3), ' => ', ' -> ');
    
    metInfo = metInfo(:, 1:12);
    metNumericHeaders = {'Charge'};
    metNumericCols = ismember(metInfo(1, :), metNumericHeaders);
    for j = 1:size(metInfo, 2)
        if ~metNumericCols(j)
            metInfo(cellfun(@(x) any(isnan(x)), metInfo(:, j)), j) = {''};
        end
    end
    metInfo = metInfo(~cellfun(@isempty, metInfo(:, 1)), :);
    model = createModel();
    
    % add metabolites
    for j = 2:size(metInfo, 1)
        model = addMetabolite(model, metInfo{j, 1}, 'metName', metInfo{j, 2}, 'metFormula', metInfo{j, 4}, ...
            'KEGGId', metInfo{j, 7}, 'ChEBIID', metInfo{j, 9}, 'PubChemID', metInfo{j, 8}, ...
            'Charge', metInfo{j, 5});
    end
    m = numel(model.mets);
    % add reactions
    [model.rxnECnumbers, model.rxnNotes, model.rxnReferences] = deal({});
    for j = 2:size(rxnInfo, 1)
        re = regexp(rxnInfo{j, 3}, '(\d*[a-zA-Z]+\w*)(?:$|\s)', 'match');
        if ~isempty(re)
            f = regexprep(rxnInfo{j, 3}, '(\d*[a-zA-Z]+\w*)(?:$|\s)', '$1\[c\] ');
            fprintf('\n\n%d\t%s\t%s\nold:  %s\nnew:  %s\n', j, rxnInfo{j, 1:3}, f);
            %         if j <= 1360
            rxnInfo{j, 3} = f;
            %         else
            %             s = input('Hit return to accept. Otherwise enter the correct formula:\n', 's');
            %             if isempty(s)
            %                 rxnInfo{j, 3} = f;
            %             else
            %                 rxnInfo{j, 3} = s;
            %             end
            %         end
        end
        model = addReaction(model, rxnInfo{j, 1}, 'reactionName', rxnInfo{j, 2}, ...
            'reactionFormula', rxnInfo{j, 3}, 'lowerBound', rxnInfo{j, 9}, 'upperBound', rxnInfo{j, 10}, ...
            'objectiveCoef', rxnInfo{j, 11}, 'subSystem', rxnInfo{j, 7}, ...
            'geneRule', rxnInfo{j, 4});
        idJ = findRxnIDs(model, rxnInfo{j, 1});
        model.rxnECnumbers{idJ, 1} = rxnInfo{j, 13};
        model.rxnNotes{idJ, 1} = rxnInfo{j, 14};
        model.rxnReferences{idJ, 1} = rxnInfo{j, 15};
        if any(strcmp(model.mets, ''))
            error('Empty mets')
        end
        %     if numel(model.mets) > m && j > 1360
        %         fprintf('\n\n%d\t%s\t%s\n%s\n%s\n', j, rxnInfo{j, 1:3}, strjoin(model.mets(m + 1:end), ', '));
        %         pause;
        %         m = numel(model.mets);
        %     end
    end
    
    if ~exist('Serythraea.mat', 'file')
        writeCbModel(model, 'format', 'mat', 'fileName', 'Serythraea.mat')
    end
end

%% Running the balancing algorithm
m = findMetIDs(model, 'dTDP-galactose[c]');
model.metFormulas{m} = 'C16H26N2O16P2';
m = findMetIDs(model, 'menaquinol[c]');
model.metFormulas{m} = '';
[metEle, elements] = getElementalComposition(model.metFormulas);
genEle = ismember(elements, {'R', 'X'});
% treat all formulae with generic groups as unknown for running the balancing algorithm
model.metFormulas(any(metEle(:, genEle), 2)) = {''};

%% Add fixes here for mass and charge balance 

%%%%%%%%%%%%%%%%% Add fixes here until 'rxnImbal' at the end is empty
%%%%%%%%%%%%%%%%%
% Unbalanced reactions should only be allowed if the reaction mechanism is
% not fully known yet, i.e., currently unbalanced in databases like MetaCyc

% example:
% r83 tagaturonate reductase: tagaturonate and altronate should switch
r = findRxnIDs(model, 'r83');  % get reaction index
m = findMetIDs(model, {'D-Tagaturonate[c]', 'D-Altronate[c]'});  % get metabolite indices
model.S(m, r) = [1; -1];  % fix the stoichiometric matrix directly

% malonyl-CoA: protonated malonyl-CoA should have 33 H
m = findMetIDs(model, 'Malonyl-CoA[c]');
model.metFormulas{m} = 'C24H33N7O19P3S';

% r117 pyruvate dehydrogenase: Ferricytochrome:Ferrocytochrome should be 2:2 to receive the correct charge? 
% ...
% ...

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

% balancing algorithm to find all hidden imbalance provided a set of correct formulae
[model2, ~, elements, metEle, rxnBal, S_fill, solInfo] = computeMetFormulae(model, 'fillMets', {'HCharge1', 'H2O'});
rxnEx = sum(model.S ~= 0, 1) <= 1;
% inactive rxns are not handled by the algorihtm
rxnActive = model.lb ~= 0 | model.ub ~= 0;
rxnImbal = model.rxns(any(abs(rxnBal) > 1e-4, 1) & ~rxnEx & rxnActive');

% automatically fill missing protons and H2O
hc = findMetIDs(model, 'H[c]');
h2o = findMetIDs(model, 'H2O[c]');
[metEleFill, elementsFill] = getElementalComposition({'HCharge1', 'H2O'}, elements);
metComp = regexp(model.mets, '\[\w\]$', 'match', 'once');
acceptAutoCorrect = false(numel(model.rxns), 1);
for j = 1:numel(model.rxns)
    % if automatic filling suggested
    if any(S_fill(:, j))
        % check whether it perfectly resolves the imbalance
        if all(abs(metEleFill' * S_fill(:, j) + metEle' * model.S(:, j)) < 1e-5)
            if all(strcmp(metComp(model.S(:, j) ~= 0), '[c]'))  % this is currently trivial since the model is not compartmentalized
                % accept
                acceptAutoCorrect(j) = true;
                model.S([hc; h2o], j) = model.S([hc; h2o], j) + S_fill(:, j);
            end
        end
    end
end

% rebalance
[model2, metFormulae, elements, metEle, rxnBal, S_fill, solInfo] = computeMetFormulae(model, 'fillMets', {'HCharge1', 'H2O'});
rxnImbal = model.rxns(any(abs(rxnBal) > 1e-4, 1) & ~rxnEx & rxnActive');
% check unbalanced reactions
printImbalance(model2, rxnImbal(1), 0, rxnBal, elements, metEle)

%%  Check ATP-generating cycles

% this is optional. For metabolites originally with generic R-, X-groups,
% model2 has all their formulae replaced by computationally defined groups,
% usually pretty close to the actual groups
model = model2;

%%%%% Add fixes to block ATP-generating cycles here (but currently seems no
%%%%% such cycle, confirm this again after finishing mass balance)


% logical index vector for exchange reactions
rxnEx = sum(model.S ~= 0, 1) <= 1;
% shut down all uptake reactions (lower bound = 0)
modelATP = changeRxnBounds(model, model.rxns(rxnEx),  0, 'l');
% non-zero ATP maintenance
modelATP = changeRxnBounds(modelATP, 'ATP maintenance',  1, 'l');
s = optimizeCbModel(modelATP, 'max', 'one');

if s.stat ~= 0
    fprintf('ATP-generating cycles!')
    surfNet(modelATP, modelATP.rxns(s.x ~= 0), 's', 0, 'f', s.x)
end

%% Call Eflux to estimate flux distributions

% format gene expression data
% ...

[solutionWT, solutionMU, fvaMin, fvaMax] = efluxFVA(model, expressionDataWT, expressionDataMU);

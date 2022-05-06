function varargout = printImbalance(modelInput, rxn, metFlagInput, rxnBalInput, elementInput, metEleInput)

persistent model
persistent metFlag
persistent rxnBal
persistent element
persistent metEle
if nargin == 1
    % query for another reaction using previous data
    if isempty(model)
        error('No previously supplied model exists. Please supply a model.')
    end
    % if only one input, it is the reaction. Use previous model and balance
    rxn = modelInput;
elseif ~isequal(model, modelInput)
    % new call
    model = modelInput;
    if nargin < 3
        metFlag = false;
    end
    if nargin < 4
        if isfield(model, 'metCharge') && ~isfield(model, 'metCharges')
            model.metCharges = model.metCharge;
        end
        if isfield(model, 'metCharges') && ~all(isnan(model.metCharges)) && ~any(contains(model.metFormulas, 'Charge'))
            % Include charges in the metFormulas field to also check their balance
            model.metFormulas = strcat(model.metFormulas, 'Charge', cellfun(@num2str, num2cell(model.metCharges), 'UniformOutPut', false));
        end
        %[rxnBal, element, metEle] = 
        [metEle, element] = getElementalComposition(model.metFormulas, [], 1);
        metUnknown = any(isnan(metEle), 2);
        rxnBal = sparse(numel(element), numel(model.rxns));
        rxnBal(:, any(model.S(metUnknown, :), 1)) = NaN;
        rxnBal(:, ~any(model.S(metUnknown, :), 1)) = metEle(~metUnknown, :)';
    elseif nargin == 6
        [rxnBal, element, metEle] = deal(rxnBalInput, elementInput, metEleInput);
    else
        error('rxnBal, element and metEle must be supplied at the same time')
    end
end


if iscell(rxn) || ischar(rxn)
    rxn = findRxnIDs(model,rxn);
end
imb = zeros(numel(element),numel(rxn));
if numel(rxn) > 1
    for j = 1:numel(rxn)
        imb(:,j) = printImbalance(rxn(j));
    end
else
    dMax = 6;
    %     printRxnFormula(model,model.rxns(rxn),1,1,metFlag);
    %     fprintf('\n');
    surfNet(model,model.rxns(rxn),metFlag);
    [r,~,e] = find(rxnBal(:,rxn));
    imb(:,1) = rxnBal(:,rxn);
    for k = 1:numel(r)
        n = e(k);
        d = 0;
        %at most keep dMax decimal places
        while n - str2double(sprintf(['%.' num2str(d) 'f'],n)) ~= 0 && d <= dMax - 1
            d = d + 1;
        end
        n = round(n, d);
        if n ~= 0
            d = 0;
            while n - str2double(sprintf(['%.' num2str(d) 'f'],n)) ~= 0 && d <= dMax - 1
                d = d + 1;
            end
        end
        fprintf('[%s]\t%', element{r(k)});
        fprintf(['%.' num2str(d) 'f\n'],n);
    end
    if isempty(r)
        fprintf('The reaction is balanced.');
    end
    fprintf('\n');
end
if nargout > 0
    varargout = {imb};
end
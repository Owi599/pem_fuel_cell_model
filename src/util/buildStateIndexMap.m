function idxMap = buildStateIndexMap(p)
    names = p.states.variableNames;
    dims  = p.states.variableDimensions;

    idxMap = struct();
    k = 1;
    for ii = 1:numel(names)
        idxMap.(names{ii}) = k:(k + dims(ii) - 1);
        k = k + dims(ii);
    end
end
function grainSize = getGrainSize(grainMap)
% getGrainSize determines grain sizes in pixels for a provided grainMap
%
% Input:
% grainMap: array of grain numbers assigned to each pixel
%
% Output:
% grainSize: a vector of grain sizes

    arguments
        grainMap double;
    end

    grainNumber = max(max(grainMap));
    grainSize = zeros(grainNumber, 1);

    for iGrain = 1:grainNumber
        grainSize(iGrain) = sum((grainMap == iGrain), 'all');
    end

end

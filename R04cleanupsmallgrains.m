function data = R04cleanupsmallgrains(data, minPoints, method)
% CLEANUPSMALLGRAINS removes small grains
% by re-assigning data points to larger neighbours
%
% Inputs:
% data: structure containing grainMap
% minPoints: integer value specyfying minimum count of points
%   for a grain to be persistent
% method: char specyfying method for reassigning points:
%   'Fast4' or 'Average8'
%
% Output:
% data: structure containing cleaned up grainMap

arguments
    data (1,1) struct;
    minPoints (1,1) double = 50;
    method (1,:) char {mustBeMember(method,{'Fast4','Average8'})} = 'Fast4';
end

% extract values from data structure
grainMap = data.grainMap;
matrices = data.orientationMatrices;
eulerAngles = data.eulerAngles;
dimensions = size(grainMap);

% direction mapping for overwrite: 1-N, 2-E, 3-S, 4-W
yMap = [-1; 0; 1; 0];
xMap = [0; 1; 0; -1];

% determine small grains to be cleaned
grainSize = getGrainSize(grainMap);
smallGrains = find(grainSize < minPoints);
largeGrains = find(grainSize >= minPoints);

% create a list of points to clean (reassign) [PTC]
PTCmap = arrayfun(@(point) any(point == smallGrains), grainMap);
[PTCy, PTCx] = find(PTCmap);

% reassign cleaned pixels to large grains
loopCount = 0;
while ~isempty(PTCy)
    sPTC = length(PTCy);
    overwriteDirection = zeros(sPTC, 1);

    % determine overwrite direction for each point
    fprintf('find overwrite directions');
    for iPoint = 1:sPTC

        % work control
        loopCount = loopCount + 1;
        if loopCount == 1000
            fprintf('cleaning points: %d', iPoint);
            loopCount = 0;
        end
            
        % create neighbors list
        nearPoints = zeros(4, 1);
        if PTCy(iPoint) > 1 % N
            nearPoints(1) = grainMap(PTCy(iPoint)-1, PTCx(iPoint));
        end
        if PTCx(iPoint) < dimensions(2) % E
            nearPoints(2) = grainMap(PTCy(iPoint), PTCx(iPoint)+1);
        end
        if PTCy(iPoint) < dimensions(1) % S
            nearPoints(3) = grainMap(PTCy(iPoint)+1, PTCx(iPoint));
        end
        if PTCx(iPoint) > 1 % W
            nearPoints(4) = grainMap(PTCy(iPoint), PTCx(iPoint)-1);
        end

        % remove grains to be cleaned from neighbors list
        uNearPoints = setdiff(nearPoints, [0; smallGrains]);
            
        % chose neighbor to overwrite
        if ~isempty(uNearPoints) % prevent the growth of false grains
            if length(uNearPoints) == 1 % single good neighbor
                newGID = uNearPoints;
            elseif length(uNearPoints) > 1 % need to choose a neighbor
                %% histcounts obcina biny na minimum: trzeba zrobic recznie binning
                histNearPoints = histcounts(uNearPoints, 'BinMethod', 'integers');
                [histGID, ~, histCount] = find(histNearPoints');
                
                %% max zwraca index tylko pierwszego maksimum
                [~, iHistCount] = max(histCount);
                    
                if numel(iHistCount) == 1 % single maximum
                    newGID = histGID(iHistCount);
                else % need to solve a draw using larger grain size
                    drawGrainSizes = grainSize(histGID);
                    [~, drawWinner] = max(drawGrainSizes);
                    newGID = histGID(drawWinner);
                end %%
            end
            overwriteDirection(iPoint) = find(nearPoints == newGID, 1);
        end
    end
        
    % update grainMap, matrices and Eulers based on overwriteDirection
    fprintf('update data');
    for iPoint = 1:sPTC
        if overwriteDirection(iPoint) > 0 % if null skip until next loop
            yNew = PTCy(iPoint) + yMap(overwriteDirection(iPoint));
            xNew = PTCx(iPoint) + xMap(overwriteDirection(iPoint));
                
            grainMap(PTCy(iPoint), PTCx(iPoint)) = grainMap(yNew, xNew);
            matrices(PTCy(iPoint), PTCx(iPoint), 1:3, 1:3) = ...
                matrices(yNew, xNew, 1:3, 1:3);
            eulerAngles(PTCy(iPoint), PTCx(iPoint), 1:3) = ...
                eulerAngles(yNew, xNew, 1:3);
        end
    end
        
    % re-create the list of points by filtering out the updated ones
    idx = overwriteDirection == 0;
    PTCy = PTCy(idx);
    PTCx = PTCx(idx);
end

% update grain IDs
for iGrain = 1:numel(largeGrains)
    grainMap(grainMap == largeGrains(iGrain)) = iGrain;
end

data.grainMap = grainMap;
data.orientationMatrices = matrices;
data.eulerAngles = eulerAngles;

finishGood;

end

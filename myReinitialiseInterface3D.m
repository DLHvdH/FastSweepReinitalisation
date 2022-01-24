function levelSet = myReinitialiseInterface3D(levelSet, h)
%{
Author: Dingeman (Danny) van der Haven
Email: dannyvdhaven@gmail.com
Last updated: 2022/01/24

Reinitialise the vales of the level set adjacent to interface by rescaling
them such that they fall between -h and h as should be. The current code
assumes that whatever distribution of level set values we have beforehand
is centered about 0.
Note: This is a hacky fix. In reality there is some probability
distribution (for the level set values in grid points next to the 
interface) that we should match. But for most purposes this fix will
suffice.

INPUT:
 levelSet (required):
    The 3D signed-distance function that is to be reinitialised.
 h (default = 1):
    The grid spacing h. The grid is assumed to be regular (i.e. the 
    spacing h is equal in all directions).

%}

% Default values
if nargin < 2
    h = 1;
end

% The exit condition set. To mark points that are next to the interface.
setQ = false(size(levelSet));

% Find the interface as indicated by a sign change in the level set.
% Compare the sign of an element to all of its neighbours, 
% this condition will be true if the sign differs.
sgnLS = sign(levelSet);
xchange = sgnLS(2:end,:,:) .* sgnLS(1:end-1,:,:) < 0;
ychange = sgnLS(:,2:end,:) .* sgnLS(:,1:end-1,:) < 0;
zchange = sgnLS(:,:,2:end) .* sgnLS(:,:,1:end-1) < 0;

% If there is any sign change, we are adjacent to an interface
setQ(1:end-1,:,:) = setQ(1:end-1,:,:) | xchange;
setQ(2:end,:,:) = setQ(2:end,:,:) | xchange;
setQ(:,1:end-1,:) = setQ(:,1:end-1,:) | ychange;
setQ(:,2:end,:) = setQ(:,2:end,:) | ychange;
setQ(:,:,1:end-1) = setQ(:,:,1:end-1) | zchange;
setQ(:,:,2:end) = setQ(:,:,2:end) | zchange;

% The vales of the level set in the cells adjacent to the interface
adjValsPos = levelSet( setQ & (sgnLS > 0) );
adjValsNeg = levelSet( setQ & (sgnLS < 0) );

% Sort these values
sortAdjValsPos = sort(adjValsPos);
sortAdjValsNeg = sort(adjValsNeg);

% Determine a cutoff based on the extreme 0.5% of the distribution.
cutoffPos = round(numel(sortAdjValsPos)/200);
cutoffNeg = max( round(numel(sortAdjValsNeg)/200), 1);

% Rescale the distrubtion such that is is uniform again.
levelSet( setQ & (sgnLS > 0) ) = ...
    adjValsPos./sortAdjValsPos(end-cutoffPos)*h;
levelSet( setQ & (sgnLS < 0) ) = ...
    -adjValsNeg./sortAdjValsNeg(cutoffNeg)*h;

% Limit the level set values to -h and +h
levelSet( setQ & (sgnLS > 0) ) = min(levelSet( setQ & (sgnLS > 0) ),h);
levelSet( setQ & (sgnLS < 0) ) = max(levelSet( setQ & (sgnLS < 0) ),-h);

%{
% Visualisation to check the result
adjValsPos = levelSet( setQ & (sgnLS > 0) );
adjValsNeg = levelSet( setQ & (sgnLS < 0) );

figure;
histogram(adjValsPos);
hold on;
histogram(adjValsNeg);
%}

end
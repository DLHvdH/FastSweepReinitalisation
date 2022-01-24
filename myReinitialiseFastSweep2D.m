function levelSet = myReinitialiseFastSweep2D(levelSet, h, maxSweeps)
%{
Author: Dingeman (Danny) van der Haven
Email: dannyvdhaven@gmail.com
Last updated: 2022/01/19

INPUT:
 levelSet (required):
    The 2D signed-distance function that is to be reinitialised.
 h (default = 1):
    The grid spacing h. The grid is assumed to be regular (i.e. the 
    spacing h is equal in all directions).
 maxSweeps (default = 8):
    The algorithm stops if the last sweep did not change the level set
    function. However, sometimes only the level set near the interface
    matters and one might want to set a maximum to the number of 
    reinitialisation sweeps. This variable sets that maximum. For no
    maximum, set maxSweeps = Inf.

REFERENCES:
 Algorithm by adam Chacon and Alexander Vladimirsk,
 SIAM Journal on Scientific computing, 2012
%}

% Default values
if nargin < 3
    maxSweeps = 8;
end
if nargin < 2
    h = 1;
end

% Dimensions of the level set field
[nx, ny] = size(levelSet); 

% The exit condition set. To mark points that are next to the interface.
setQ = false(size(levelSet));

% Find the interface as indicated by a sign change in the level set.
% Compare the sign of an element to all of its neighbours, 
% this condition will be true if the sign differs.
sgnLS = sign(levelSet);
xchange = sgnLS(2:end,:) .* sgnLS(1:end-1,:) < 0;
ychange = sgnLS(:,2:end) .* sgnLS(:,1:end-1) < 0;

% If there is any sign change, we are adjacent to an interface
setQ(1:end-1,:) = setQ(1:end-1,:) | xchange;
setQ(2:end,:) = setQ(2:end,:) | xchange;
setQ(:,1:end-1) = setQ(:,1:end-1) | ychange;
setQ(:,2:end) = setQ(:,2:end) | ychange;

% Set all elements not adjacent to the interface to sgn(levelSet) * Inf
levelSet(~setQ) =  sgnLS(~setQ).*Inf;

%figure; histogram(levelSet(setQ));

% Start and end points for the x, y, and z loops depending on the
% direction (drct) in which we want to loop.
xLoop = [1 nx; 1 nx; nx 1; nx 1];
yLoop = [1 ny; ny 1; ny 1; 1 ny];

% Initialise sweep number and a tracking to see if the level set changed
% (i.e. was updated) or not.
sweepNr = 0;
changed = true;

% Doing the reinitialisation sweeps until the level set no longer changes
% or we have reached the maximum number of sweeps.
while (changed == true) && (sweepNr < maxSweeps)
    
    % Set loop condition to false and determine direction
    changed = false;
    d = mod(sweepNr,4) + 1;

    % Apply sweep
    for i = xLoop(d,1): sign(xLoop(d,2)-xLoop(d,1)) :xLoop(d,2)
        for j = yLoop(d,1): sign(yLoop(d,2)-yLoop(d,1)) :yLoop(d,2)

            [levelSet, changed] = ...
                upwindLevelUpdate(levelSet,h,setQ,changed,i,j);

        end % For i
    end % For k

    % Update sweep number
    sweepNr = sweepNr + 1;

end % While

%disp("Total number of sweeps: "+num2str(sweepNr));

end



function [levelSet, changed] = ...
    upwindLevelUpdate(levelSet, h, setQ, changed, i, j)
%{
Do a first-order upwind update of the level set
%}
% If not an interface point
if ~setQ(i,j)

    % Select point to consider
    lev = levelSet(i,j);
    sz = size(levelSet);

    % Select neighbours (and check if we're not at a boundary)
    if i == 1
        levXpos = levelSet(i+1,j);
        levXneg = levXpos;
    elseif i == sz(1)
        levXneg = levelSet(i-1,j);
        levXpos = levXneg;
    else
        levXpos = levelSet(i+1,j);
        levXneg = levelSet(i-1,j);
    end
    
    if j == 1
        levYpos = levelSet(i,j+1);
        levYneg = levYpos;
    elseif j == sz(2)
        levYneg = levelSet(i,j-1);
        levYpos = levYneg;
    else
        levYpos = levelSet(i,j+1);
        levYneg = levelSet(i,j-1);
    end

    if (lev > 0)

        % Shortcut method, select lowest-absolute-value neighbours
        lev1 = min(levXpos,levXneg);
        lev2 = min(levYpos,levYneg);
        levHat = findRootEikonalND([lev1,lev2],h,true);

        % Update level set if distance value is lower
        if (levHat < lev)
    		levelSet(i,j) = levHat;
    		changed = true;
        end

    elseif (lev < 0)

        % Shortcut method, select lowest-absolute-value neighbours
        lev1 = max(levXpos,levXneg);
        lev2 = max(levYpos,levYneg);
        levHat = findRootEikonalND([lev1,lev2],h,false);

        % Update level set if distance value is lower
        if (levHat > lev)
    	    levelSet(i,j) = levHat;
    	    changed = true;
        end

    end

end % If ~setQ(i,j,k)

end

function levHat = findRootEikonalND(levs, h, rootChoice)
%{
The solution to the ND Eikonal equation.
If the level set value  to be updated is positive rootChoice = true.
Otherwise, for negative level set values rootChoice = false.
See: https://en.wikipedia.org/wiki/Eikonal_equation
%}

% Compute the discriminant D
N = numel(levs);
D = sum(levs)^2 + N*(h^2-sum(levs.^2));
% For a gradient with a non-unity magnitude f, substitute h^2 for h^2/f^2.

% Check the possible solutions
% U = 1/N sum(u_j) +/- 1/N * sqrt( sum(u_j)^2 + 1/N (h^2/f^2-sum(u_j^2)) )
% For a signed-distance funciton such as a level set the magnitude f of 
% the gradient is 1.
if D > 0
    % A real solution exists
    if rootChoice
        levHat = 1/N * sum(levs) + 1/N * sqrt(D);
    else
        levHat = 1/N * sum(levs) - 1/N * sqrt(D);
    end
else
    % A real solution does not exist, one of the partial derivatives is 
    % equal to zero. Perform a lower dimensional update and chose the 
    % smallest value.

    % Initalise the new estimated level value depending on the sign
    if rootChoice
        levHat = Inf;
    else
        levHat = -Inf;
    end

    % Recurvsively search all lower-dimensional solutions for the minimum
    % (or maximum if negative) level set update that can be obtained.
    for i = 1:N
        levsTemp = levs;
        levsTemp(i) = [];
        levHatTemp = findRootEikonalND(levsTemp,h,rootChoice);
        if rootChoice
            levHat = min(levHat, levHatTemp);
        else
            levHat = max(levHat, levHatTemp);
        end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    pcaMotionAnalysis
%
% FUNCTION:         pcaMotionAnalysis(positionData)
%
% DESCRIPTION:      Uses PCA to determine direction and magnitude of brain
%                   motion for easy arrow plot
%
% INPUT:
%
% VARIABLES:
%
% OUTPUT:
%
% FUNCTIONS USED:
%
% LIBARIES USED:
%
% NOTES:
%
% WRITTEN BY:       Spencer Garborg 4/16/20
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function motionVec = pcaMotionAnalysis(positionData,swapTF)

if ~exist('swapTF','var')
    swapTF = false;
end
% Mean center data
positionDataMC(:,1) = positionData(:,1) - mean(positionData(:,1));
positionDataMC(:,2) = positionData(:,2) - mean(positionData(:,2));

% Calculate covariance matrix of data for rotation information data
C = cov(positionDataMC);

% Diagonalize covariance matrix to decorrelate new variables through
% rotation. V are eigenvectors belonging to diagonalized covariance matrix
% to express correlation between old and new data. D are eigenvalues that
% describe the variance within the calculated principal components. 

[V,D] = eig(C);

% Determine the principal component with largest variance to make the
% directional unit vec of movement
if D(1,1) > D(2,2)
    if swapTF
        motionVec = [V(1,2) V(2,2)];
    else
        motionVec = [V(1,1) V(2,1)];
    end
else
    if swapTF
        motionVec = [V(1,1) V(2,1)];
    else
        motionVec = [V(1,2) V(2,2)];
    end
end

% Determine direction of principal component vector based on mean of data
% if mean(positionData(:,1)) < 0 && motionVec(1) > 0
%     motionVec = motionVec * -1;
% end
if motionVec(2) < 0
    motionVec = motionVec * -1;
end

% Calculate magnitude based on average of furthest 20% of movement data
% points from origin
numMags = ceil(size(positionData,1)*.2);
distanceVec = ((positionData(:,1).^2) + (positionData(:,2).^2)).^.5;
magnitude = mean(maxk(distanceVec,numMags));
motionVec = [motionVec * magnitude, motionVec * prctile(distanceVec,75), motionVec * prctile(distanceVec,95)];
end
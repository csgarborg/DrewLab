% How-to for running full 2P2L analysis

% 1: run single tracking on each layer 1-10 times

%           motionTracking2P2Layer(tifFileName,calibrationFileString,updateSearchTF,medFiltTF,saveOutputVideoTF,threeStepTF,compiledTifTF,tempMedFiltTF,targetAvgNum,framesPerSecond,analogSampleRate,objMag,digMag,turnabout,layer,hemisphere,commentString,tifFrameBounds)

% 2: run combination code on each layer

%           combineMotionTracking(fileName,versionVec,mode)

% 3: combine the layer combinations for final plots with skull as first one

%           plotMotionTracking2Layer(matFileName1,matFileName2)
function [tdmsDataStruct] = fixTDMSFieldNames(tdmsDataStruct)
fieldNames = fieldnames(tdmsDataStruct.Digital_Data);
if length(fieldNames) < 3
    interleavedData = tdmsDataStruct.Digital_Data.Side_Diff;
    digitalNameList = {'Right_Whisker_Diff','Left_Whisker_Diff','Side_Diff','Front_Diff','Bottom_Diff','Respiration_Sum'};
    for n = 1:6
        tdmsDataStruct.Digital_Data.(digitalNameList{n}) = interleavedData(n:6:length(interleavedData));
    end
end
end
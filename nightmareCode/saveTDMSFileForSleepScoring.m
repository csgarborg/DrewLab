function saveTDMSFileForSleepScoring(filePath)

tdmsDataStruct = openTDMS(filePath);
save([filePath(1:end-4) 'mat'],'tdmsDataStruct')
end
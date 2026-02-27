function saveTifStack(imageStack, saveFileName)
if exist(saveFileName,'file')
    delete(saveFileName)
end
for n = 1:size(imageStack,3)
    if n == 1
        imwrite(imageStack(:,:,n),saveFileName);
    else
        imwrite(imageStack(:,:,n),saveFileName,'WriteMode','append');
    end
end
end
for i = [8 9 13 15 16 22 23]
    for j = 1:20
        if i < 10
            x = ['0' num2str(i)];
        else
            x = num2str(i);
        end
        if j < 10
            y = ['0' num2str(j)];
        else
            y = num2str(j);
        end
        for k = 1:3
            if ~exist(['D:\22-08-' x '_MouseExp\2208' x '_0' y '_processed_Layer2_' num2str(k) '.mat'],'file')
                continue
            else
                load(['D:\22-08-' x '_MouseExp\2208' x '_0' y '_processed_Layer2_' num2str(k) '.mat'])
                movementData.runDate = ['2208' x];
                movementData.runNumber = num2str(j);
                save(['D:\22-08-' x '_MouseExp\2208' x '_0' y '_processed_Layer2_' num2str(k) '.mat'],'movementData');
                disp([x y ' ' movementData.mouseNumber ' ' num2str(k)])
            end
        end
    end
end
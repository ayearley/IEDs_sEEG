function[] = TWChar(ptIDs)
% TWChar will calculate the speed of all IEDs by both Euclidean 
% distance and path length
%
% inputs:   1) patient identifier

% Authors [AGY:20221226]

topDir = %%Input top directory.
Results = [];

for pt = 1:length(ptIDs)

    ptID = convertStringsToChars(ptIDs(pt));

    topDirData = [topDir ptID '/Data/'];
    topDirIED = [topDir ptID '/IEDs/'];

    IEDFiles = dir(topDirIED);
    dataFiles = dir(topDirData);

    IEDSpeedsEucD = [];
    IEDSpeedsPL = [];

    files = [];
    segments = [];
    IEDs = [];

    for fl=1:length(IEDFiles)
        fileName = IEDFiles(fl).name;
        if(~contains(fileName, ptID))
            continue;
        end
        vals = regexp(extractAfter(fileName,'file'),'\d+','match');
        files = [files vals(1)];
        segments = [segments vals(2)];
        IEDs = [IEDs vals(3)];
    end

    TotalIED = length(IEDs);

    for fl=1:length(dataFiles)

        if(~contains(dataFiles(fl).name, ptID))
            continue;
        end

        IEDseg = load([dataFiles(fl).folder '/' dataFiles(fl).name]);
        fileName = dataFiles(fl).name;
        vals = regexp(extractAfter(fileName,'file'),'\d+','match');
        for x = 1:length(segments)
            if(str2double(files(x)) == str2double(vals(1)) && str2double(segments(x)) == str2double(vals(2)))
                SigEucD = IEDseg.ptResults.sigMinModel(str2double(IEDs(x)));
                SigNumCon = IEDseg.ptResults.sigMinModel_numCons(str2double(IEDs(x)));
                SigPL = IEDseg.ptResults.sigMinModel_pathLen(str2double(IEDs(x)));

                %If the IED has been determined to be a TW by any method,
                %measure speed from euclidean distance
                if (SigEucD == 1)
                    IEDspeed = IEDseg.ptResults.minsLM(str2double(IEDs(x))).LM.Coefficients(2,1);
                    IEDspeed = 1/table2array(IEDspeed) * 1000/10; %Convert to mm/msec then to cm/s
                    IEDSpeedsEucD(x) = IEDspeed;
                end
                if (SigPL == 1)
                    IEDspeed = IEDseg.ptResults.minsLM_pathLen(str2double(IEDs(x))).LM.Coefficients(2,1);
                    IEDspeed = 1/table2array(IEDspeed) * 1000/10; %Convert to mm/msec then to cm/s
                    IEDSpeedsPL(x) = IEDspeed;
                end
            end
        end
    end
    IEDSpeedsEucD = IEDSpeedsEucD(IEDSpeedsEucD > 0);
    medianSpeedEucD = median(IEDSpeedsEucD);

    IEDSpeedsPL = IEDSpeedsPL(IEDSpeedsPL > 0);
    medianSpeedPL = median(IEDSpeedsPL);

    if(pt == 1)
        Results = {ptID, medianSpeedEucD, medianSpeedPL};
        Results = array2table(Results);
        Results.Properties.VariableNames(1:3) = {'Patient ID', 'Median Speed (Euclidean Distance)','Median Speed (Path Length)'};
    else
        Results2 = {ptID, medianSpeedEucD, medianSpeedPL};
        Results2 = array2table(Results2);
        Results2.Properties.VariableNames(1:3) = {'Patient ID', 'Median Speed (Euclidean Distance)','Median Speed (Path Length)'};
        Results = [Results; Results2]
    end
end

writetable(Results,'IED Speed.csv')



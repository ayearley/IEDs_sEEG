function[] = IEDChanOrder(ptIDs, TW_only)
% IEDChanOrder will calculate the channels in which each IED was detected 
% in sequential order
%
% inputs:   1) patient identifier
% inputs:   1) if only IED travelling waves should be included [true] or
% all IEDs [false]

% Authors [AGY:20230109]

topDir = %%Input top directory.
Results = [];

for pt = 1:length(ptIDs)

    ptID = convertStringsToChars(ptIDs(pt));

    topDirData = [topDir ptID '/Data/'];
    topDirIED = [topDir ptID '/IEDs/'];

    IEDFiles = dir(topDirIED);
    dataFiles = dir(topDirData);


    files = [];
    segments = [];
    IEDs = [];
    iIndex = 1;

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

                if(TW_only && SigPL == 0 && SigNumCon == 0 && SigEucD == 0)
                    continue
                else
                    IEDmindices = IEDseg.ptResults.IEDdata(str2double(IEDs(x))).mindices;
                    IED_chans = IEDseg.ptResults.IEDdata(str2double(IEDs(x))).chans;

                    [nums,IEDchansofMin] = sort(IED_chans.*IEDmindices);
                    check = find(nums == 0);
                    IEDchansofMin(check) = [];

                    Results.pt = ptID;
                    Results.chanOrder{iIndex,1} = IEDchansofMin;
                    iIndex = iIndex + 1;
                end
            end
        end
    end


    for i = 1:iIndex-1
        line = Results.chanOrder{i,1};
        path = sprintf('%s_results.csv', ptID);
        dlmwrite(path, line, '-append');
    end

end


function[] = CountSig(ptID)
% CountSig will count the # IEDs the # of significant models by type of TW
%
% inputs:   1) patient identifier

% Authors [AGY:20221215]

topDir = %%Input top directory.

topDirData = [topDir ptID '/Data/'];
topDirIED = [topDir ptID '/IEDs/'];

IEDFiles = dir(topDirIED);
dataFiles = dir(topDirData);

numSig = 0;
numSigEucD = 0;
numSigNumCon = 0;
numSigPL = 0;
SDmult = '';

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
            SDmult = IEDseg.ptResults.IEDdata(str2double(IEDs(x))).SDMultiplier;
            SigEucD = IEDseg.ptResults.sigMinModel(str2double(IEDs(x)));
            SigNumCon = IEDseg.ptResults.sigMinModel_numCons(str2double(IEDs(x)));
            SigPL = IEDseg.ptResults.sigMinModel_pathLen(str2double(IEDs(x)));

            if (SigEucD == 1 || SigNumCon == 1 || SigPL == 1)
                numSig = numSig + 1;
            end

            numSigEucD = numSigEucD + SigEucD;
            numSigNumCon = numSigNumCon + SigNumCon;
            numSigPL = numSigPL + SigPL;
        end
    end
end

pTW = numSig/TotalIED * 100;

Results = {ptID, SDmult, length(unique(files)), length(unique(segments)), TotalIED, numSig, pTW, numSigEucD, numSigNumCon, numSigPL};
Results = array2table(Results);
Results.Properties.VariableNames(1:10) = {'Patient ID', 'SD Threshold' 'No. Files', 'No. Segments', 'No. IEDs Detected', 'No. IEDs as TW', '%IEDs as TW', 'No. Euclidean Distance', 'No. Number Connections','No. Path Length'};
writetable(Results,'IED Summary.csv')




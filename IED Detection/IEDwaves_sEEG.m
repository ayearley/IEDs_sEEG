function [] = IEDwaves_sEEG(ptID,multiplier,visualBadChanlabels)
% IEDWAVES_SEEG looks for IED traveling waves on sEEG
%
% inputs:   1) patient identifier
%           2) SD multiplier (for magnitude of peak detection)
%           3) a cell array of bad channel labels.

% Authors [EHS:20211007, AGY:20221129]

%% Loading data

topDir = [ptID '/Cadwell/']; % Directory containing edf files 
dataFiles = subdir([topDir 'fullDurationExport.edf']);

segmentSize = 600; % in s (10 min)
Fo = 2e3; % in samples
Fs = 500; 

% Distances from stim to response electrodes
LR = '5';

% Loading electrode data
load(fullfile('Data',ptID,'Imaging','Registered','Electrodes.mat'))
load(fullfile('Data',ptID,'Imaging','Registered','ChannelMap2.mat'))

% Loading number of connections data
numConnsFile = fullfile(ptID,'Imaging','Probtrackx',sprintf('FA_%s_LR_%s',ptID,LR),'fdt_network_matrix');
waytotalFile = fullfile(ptID,'Imaging','Probtrackx',sprintf('FA_%s_LR_%s',ptID,LR),'waytotal');
load(numConnsFile)
load(waytotalFile)

% Loading path length data
pathLenFile = fullfile(ptID,'Imaging','Probtrackx',sprintf('FA_%s_LR_%s',ptID,LR),'fdt_network_matrix_lengths');
load(pathLenFile)

% Redefining channel numbers and labels in terms of edf order
chanIdcs = ismember([ElecMapRaw{:,3}],DepthElec);
[inclChans, idc] = sort(cell2mat(ElecMapRaw(chanIdcs,3)));
chanLabels = ElecMapRaw(idc,1);

% Removing NaNs from ElecXYZProj
ElecXYZProj(isnan(ElecXYZProj(:,1)),:) = [];


% Removing obviously bad channels
for tr = length(visualBadChanlabels):-1:1
    if ~isempty(inclChans(strcmp(deblank(chanLabels),visualBadChanlabels(tr))))
        bads(tr) = find(strcmp(deblank(chanLabels),visualBadChanlabels(tr)));
    end
end

% Adjusting all matrices to be without bad channels
if exist('bads','var')
    ElecXYZProj(bads,:) = [];
    chanLabels(bads) = [];
    inclChans(bads) = [];
    fdt_network_matrix(bads,:) = [];
    fdt_network_matrix(:,bads) = [];
    waytotal(bads,:) = [];
end
 
nChans = length(inclChans);

num_files = 0;
num_trueSegments = 0;

% Looping over data files
for fl = 1:length(dataFiles)

    % Determine total number of 10min segments in a file
    try
        [hdr_overall, ~] = fastEDFReadMod('File',dataFiles(fl).name,'Range',[1 2]);
        total_sg = floor(hdr_overall.records*hdr_overall.duration/segmentSize);
        num_trueSegments = num_trueSegments + total_sg;
    catch
        continue;
    end

    for sg = 1:total_sg
        clear -regexp ptResults ^IED rmsvals

        try
            if sg==1
                [hdr,D] = fastEDFReadMod('File',dataFiles(fl).name,'Range',[1 Fo*segmentSize]);
            else
                [hdr,D] = fastEDFReadMod('File',dataFiles(fl).name,'Range',[(sg-1)*Fo*segmentSize (sg)*Fo*segmentSize]);
            end
            fileFlag = true;
        catch
            fileFlag = false;
            fprintf('\n unable to load file segment %d, for file %d, for patient: %s',sg,fl,ptID)
        end

        if fileFlag
            % Determining the start time of the EDF
            fileStart = datetime(str2double(['20' hdr.startdate(7:8)]),...
                str2double(hdr.startdate(4:5)),str2double(hdr.startdate(1:2)),...
                str2double(hdr.starttime(1:2)),str2double(hdr.starttime(4:5)),...
                str2double(hdr.starttime(7:8)));

            % Resampling ECoG
            tmp = resample(double(D),Fs,Fo)';

            % Counting number of data-containing segments
            if(sg == 1 || round(mode(mean(tmp, 2))) == 1)
                num_trueSegments = num_trueSegments - 1;
            end

            % Filtering & De-meaning IED detection signal
            [b,a] = butter(4,[1 50]./(Fs/2)); %changed from [20 40] to [1 50] by ALEX

            for ch = 1:nChans
                real_chan = inclChans(ch);
                updateUser('filtering channels',ch,20,nChans)
                tmpSig = tmp(real_chan,:)-mean(tmp(real_chan,:));
                IEDsignal(ch,:) = filtfilt(b,a,tmpSig);
            end

            % Common average referencing
            IEDsignal = IEDsignal - mean(IEDsignal,1,'omitnan');

            clear tmp D                

            %% IED detection: detecting peaks
            for ch3 = nChans:-1:1
                minPkProm = multiplier*std(abs(IEDsignal(ch3,:)));
                med_val = multiplier*median(std(abs(IEDsignal),0, 2));
                if(minPkProm < med_val)
                    minPkProm = med_val;
                end
                minPkHt = 2/3*minPkProm;
                [IEDmaximaVals,IEDmaxima,IEDwidths,~] = findpeaks(smoothdata(IEDsignal(ch3,:)','gaussian',10)','MinPeakProminence',minPkProm, 'MinPeakHeight',minPkHt,'MaxPeakWidth',50,'WidthReference','halfprom','Threshold', 1e-4,'MinPeakDistance',3.01); %Consider removing filters
                IEDmaximaValsTotal{ch3} = [IEDmaximaVals];
                IEDmaximaTotal{ch3} = [IEDmaxima];
                IEDwidthsTotal{ch3} = [IEDwidths];
                IEDchansTotal{ch3} = ch3*ones(1,length(IEDwidths));

            end

            % Aggregates peak locations and widths into 1 set
            IEDmaximaTotal2 = [];
            IEDmaximaValsTotal2 = [];
            IEDwidthsTotal2 = [];
            IEDchansTotal2 = [];
            for x=1:nChans
                IEDmaximaTotal2 = [IEDmaximaTotal2 cell2mat(IEDmaximaTotal(x))];
                IEDmaximaValsTotal2 = [IEDmaximaValsTotal2 cell2mat(IEDmaximaValsTotal(x))];
                IEDwidthsTotal2 = [IEDwidthsTotal2 cell2mat(IEDwidthsTotal(x))];
                IEDchansTotal2 = [IEDchansTotal2 cell2mat(IEDchansTotal(x))];
            end

            % Sorts all IED maximas in ascending order
            sorted_IEDmaximaTotal2 = sort(IEDmaximaTotal2);

            % Consolidates overlapping peaks in multiple channels
            for x=1:length(sorted_IEDmaximaTotal2)-1
                y = x+1;
                while(sorted_IEDmaximaTotal2(y)-sorted_IEDmaximaTotal2(x)<=3)
                    sorted_IEDmaximaTotal2(y)=sorted_IEDmaximaTotal2(x);
                    if(y == length(sorted_IEDmaximaTotal2))
                        break;
                    end
                    y=y+1;
                end
                x = y-1;
            end

            % Creates a new dataset of peaks present in 5 or more channels
            U = unique(sorted_IEDmaximaTotal2);

            greater9 = find(hist(sorted_IEDmaximaTotal2,U)>=5);
            IEDmaxima = [1:length(greater9)];
            IEDmaximaVals = [1:length(greater9)];
            IEDwidths = [1:length(greater9)];

            for x=1:length(greater9)
                try
                    index = find(IEDmaximaTotal2 == U(greater9(x)),1);
                    IEDmaxima(x) = IEDmaximaTotal2(index);
                    IEDmaximaVals(x) = IEDmaximaValsTotal2(index);
                catch
                end
            end

            % Excluding IEDs within a half second of each other
            IEDmaxima((IEDmaxima-(Fs/2))<=0) = [];
            IEDmaxima((IEDmaxima+(Fs/2))>length(IEDsignal)) = [];

            [~,~,G]=uniquetol(IEDmaxima,Fs/4,'DataScale',1);
            G = G.';
            if(length(IEDmaxima)>0)
                IEDmaxima = splitapply(@max,IEDmaxima,G);
            end

            IED_chans = zeros(length(IEDmaxima),nChans);
            numC = 0;

            %% Looping over IED detections.
            for ied = length(IEDmaxima):-1:1

                % Identifying channels with IED maxima
                loc = IEDmaxima(ied);
                max_locs = find(IEDmaximaTotal2<=loc+3 & IEDmaximaTotal2>=loc);
                for x=1:length(max_locs)
                    IED_chans(ied,IEDchansTotal2(max_locs(x))) = 1;
                    numC = numC+1;
                end

                % Looping over channels to filter data
                for ch2 = nChans:-1:1

                    IEDdata(ied).data(ch2,:) = IEDsignal(ch2,IEDmaxima(ied)-(Fs/2):IEDmaxima(ied)+(Fs/2));

                    % Finding timings of IED minima and peaks
                    try
                        [IEDmins(ch2),IEDmindices(ch2)] = min(smoothdata(IEDdata(ied).data(ch2,(Fs/2):Fs),'gaussian',10));
                    catch
                        [IEDmins(ch2),IEDmindices(ch2)] = min(smoothdata(IEDdata(ied).data(ch2,250:500),'gaussian',10));
                    end
                end

                % Grouping data
                IEDdata(ied).file = dataFiles(fl).name;
                IEDdata(ied).segment = sg;
                IEDdata(ied).segmentSize = segmentSize;
                IEDdata(ied).SDMultiplier = multiplier;
                IEDdata(ied).maxindex = IEDmaxima(ied);
                IEDdata(ied).mindices = IEDmindices;
                IEDdata(ied).chans = IED_chans(ied,:);

                %% Distances
                % Look at distance from channel with earliest minima that detected an IED
                [~,sourceChan] =  min(IEDmindices);                
                [~,chansofMin] = sort(IEDmindices); %Earliest channel with a min
                for x = 1:length(chansofMin)
                    if(IED_chans(ied,chansofMin(x)) == 1)
                        sourceChan = chansofMin(x);
                        break;
                    end
                end

                % Creating a set of IED channels ordered by minima timing
                [vals,IEDchansofMin] = sort(IED_chans(ied,:).*IEDmindices);
                check = find(vals == 0);
                IEDchansofMin(check) = [];

                % Sets starting time of regression to point of first minima
                minIndex = IEDmindices(sourceChan);
                rChans = find(IEDmindices >= minIndex);
                numRChans = length(rChans);
                IEDmindices = IEDmindices(rChans);

                [~,channelMindices] = sort(IEDmins); %Channel w/ largest min

                % Euclidean distances
                D_euclidean = sqrt(sum((ElecXYZProj(rChans,:) - repmat(ElecXYZProj(sourceChan,:),numRChans,1)).^2,2));


                % Only picking included channels
                D_probCon = log10(fdt_network_matrix(sourceChan,rChans)./waytotal(rChans)');

                D_probCon(isinf(D_probCon)) = 0;

                % Only picking included channels
                D_pathLen = fdt_network_matrix_lengths(sourceChan,rChans);

                % Updating user every 50 IEDs
                updateUser('processed IED',ied,50,length(IEDmaxima))

                % Data storage structure
                ptResults.pt = ptID;
                ptResults.chanLabel{ied} = chanLabels(sourceChan);
                ptResults.chanNumber(ied) = inclChans(sourceChan); %According to edf channel number
                ptResults.fileNumber = fl;
                ptResults.segmentNumber = sg;
                ptResults.edfHeader = hdr;
                ptResults.chanOrder = IEDchansofMin;


                %% Running regression models
                % Linear  models for distance versus IED mins
                % Times 2 to convert to seconds
                minsLM = fitlm(D_euclidean,IEDmindices*2,'RobustOpts','on');

                minsLM_numCons = fitlm(D_probCon,IEDmindices*2,'RobustOpts','on');                      

                minsLM_pathLen = fitlm(D_pathLen,IEDmindices*2,'RobustOpts','on');

                %% Saving models
                % LMs
                ptResults.FileStartFromEDF = fileStart;
                ptResults.detectionTimesSec_ReFileStart = (segmentSize*(sg-1))+(IEDmaxima./Fs);
                ptResults.minsLM(ied).LM = minsLM;
                ptResults.minsLM_probCon(ied).LM = minsLM_numCons;
                ptResults.minsLM_pathLen(ied).LM = minsLM_pathLen;

                % Saving bad Channels
                ptResults.postHocBadChans_range(ied).labels = chanLabels(outliers(range(IEDdata(ied).data,2)));

                % IED data
                ptResults.IEDdata(ied) = IEDdata(ied);


                %% Visualize data for each IED
                delayMap = cool(nChans);

                if (minsLM.Coefficients{2,3}>=2.5 || (minsLM_numCons.Coefficients{2,3}<=-2.5 && minsLM_pathLen.Coefficients{2,3}>=2.5))
                    plotIED = true;
                else
                    plotIED = false;
                end

                set(0,'DefaultFigureVisible','off')

                if plotIED
                    figure(ied);
                    ax1 = subplot(4,2,1);
                    hold on

                    for ch3 = 1:nChans
                        if(IED_chans(ied, ch3)==1)
                            plot(linspace(0,1,length(IEDdata(ied).data)),smoothdata(IEDdata(ied).data(ch3,:),'gaussian',10),'color',[0 0 0],'linewidth',0.5) %'color',delayMap(channelMindices(ch3),:),channelMindices(ch3)
                        end
                    end
                    hold off
                    axis tight
                    title(sprintf('patient: %s, file: %d, segment: %d, IED number: %d',ptID,fl,sg,ied))
                    ylabel('voltage (uV)')
                    xlabel('time (s)')

                    ax3 = subplot(4,2,3);
                    mP = plotAdded(minsLM);
                    mP(1).Marker = '.';
                    mP(1).MarkerEdgeColor = 'k';
                    legend off
                    axis tight
                    ylabel('timing of IED minima')
                    xlabel('euclidean distance')
                    if minsLM.Coefficients{2,4}<=0.05 && ~isinf(minsLM.Coefficients{2,3}) && sign(minsLM.Coefficients{2,3})>0 %only looking at positive effects.
                        title(sprintf('significant: (t(%d) = %.2f, p = %.2f)',minsLM.DFE,minsLM.Coefficients{2,3},minsLM.Coefficients{2,4}))
                        ptResults.sigMinModel(ied) = true;
                        ptResults.signMinModel(ied) = sign(minsLM.Coefficients{2,3});
                    else
                        title(sprintf('NOT significant: (t(%d) = %.2f, p = %.2f)',minsLM.DFE,minsLM.Coefficients{2,3},minsLM.Coefficients{2,4}))
                        ptResults.sigMinModel(ied) = false;
                        ptResults.signMinModel(ied) = sign(minsLM.Coefficients{2,3});
                    end

                    ax5 = subplot(4,2,5);
                    mP = plotAdded(minsLM_numCons);
                    mP(1).Marker = '.';
                    mP(1).MarkerEdgeColor = 'k';
                    legend off
                    axis tight
                    ylabel('timing of IED minima')
                    xlabel('number of connections')
                    if minsLM_numCons.Coefficients{2,4}<=0.05 && ~isinf(minsLM_numCons.Coefficients{2,3}) && sign(minsLM_numCons.Coefficients{2,3})<0
                        title(sprintf('significant: (t(%d) = %.2f, p = %.2f)',minsLM_numCons.DFE,minsLM_numCons.Coefficients{2,3},minsLM_numCons.Coefficients{2,4}))
                        ptResults.sigMinModel_numCons(ied) = true;
                        ptResults.signMinModel_numCons(ied) = sign(minsLM_numCons.Coefficients{2,3});
                    else
                        title(sprintf('NOT significant: (t(%d) = %.2f, p = %.2f)',minsLM_numCons.DFE,minsLM_numCons.Coefficients{2,3},minsLM_numCons.Coefficients{2,4}))
                        ptResults.sigMinModel_numCons(ied) = false;
                        ptResults.signMinModel_numCons(ied) = sign(minsLM_numCons.Coefficients{2,3});
                    end

                    ax7 = subplot(4,2,7);
                    mP = plotAdded(minsLM_pathLen);
                    mP(1).Marker = '.';
                    mP(1).MarkerEdgeColor = 'k';
                    legend off
                    axis tight
                    ylabel('timing of IED minima')
                    xlabel('path length')
                    if minsLM_pathLen.Coefficients{2,4}<=0.05 && ~isinf(minsLM_pathLen.Coefficients{2,3}) && sign(minsLM_pathLen.Coefficients{2,3})>0
                        title(sprintf('significant: (t(%d) = %.2f, p = %.2f)',minsLM_pathLen.DFE,minsLM_pathLen.Coefficients{2,3},minsLM_pathLen.Coefficients{2,4}))
                        ptResults.sigMinModel_pathLen(ied) = true;
                        ptResults.signMinModel_pathLen(ied) = sign(minsLM_pathLen.Coefficients{2,3});
                    else
                        title(sprintf('NOT significant: (t(%d) = %.2f, p = %.2f)',minsLM_pathLen.DFE,minsLM_pathLen.Coefficients{2,3},minsLM_pathLen.Coefficients{2,4}))
                        ptResults.sigMinModel_pathLen(ied) = false;
                        ptResults.signMinModel_pathLen(ied) = sign(minsLM_pathLen.Coefficients{2,3});
                    end

                    subplot(3,2,[2 4 6])
                    hold on
                    imagesc(linspace(0,1,length(IEDdata(ied).data)),1:nChans,IEDdata(ied).data(chansofMin,:));
                    hold off
                    set(gca, 'YTick',1:nChans,'YTickLabel',chanLabels(chansofMin),'FontSize',8)
                    axis xy tight
                    ylabel('channels ordered by minima timing')
                    xlabel('time(ms)')

                    colormap(turbo)

                    % saving figure
                    halfMaximize(ied,'left')
                    set(ied,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0]) %make it smaller to save on the pdf (ALEX DID IT)
                    
                    saveas(ied,fullfile([ptID '/IEDs/'],sprintf('%s_file%d_segment%d_IED%d_distanceFromChannel%d.pdf',ptID,fl,sg,ied,sourceChan)))

                    close(ied)
                end
            end
        end
    end % looping over segments
end % looping over files.

fprintf('\nCode complete for patient %s\n',ptID)

if isfield('ptResults','sigMinModel')
    % IED minima results
    fprintf('\n%d of %d IEDs significant using timing of IED minima models (Euclidean Distance)\n',sum(ptResults.sigMinModel),length(IEDdata))
    fprintf('\n%d of %d IEDs significant using timing of IED minima models (Number of Connections)\n',sum(ptResults.sigMinModel_numCons),length(IEDdata))
    fprintf('\n%d of %d IEDs significant using timing of IED minima models (Path Length)\n',sum(ptResults.sigMinModel_pathLen),length(IEDdata))

    fprintf('\n%d of %d IEDs significant using timing of IED minima models (Euclidean Distance & Number of Connections)\n',sum(ptResults.sigMinModel),length(IEDdata))
    fprintf('\n%d of %d IEDs significant using timing of IED minima models (Euclidean Distance & Path Length)\n',sum(ptResults.sigMinModel & ptResults.sigMinModel_numCons),length(IEDdata))
    fprintf('\n%d of %d IEDs significant using timing of IED minima models (Number of Connections & Path Length)\n',sum(ptResults.sigMinModel_numCons & ptResults.sigMinModel_pathLen),length(IEDdata))
    fprintf('\n%d of %d IEDs significant using timing of IED minima models (All three dsitance metrics)\n',sum(ptResults.sigMinModel & ptResults.sigMinModel_numCons & ptResults.sigMinModel_pathLen),length(IEDdata))

    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

    % IED maxima results
    fprintf('\n%d of %d IEDs significant using timing of beta power maxima models (Euclidean Distance)\n',sum(ptResults.sigBetaModel),length(IEDdata))
    fprintf('\n%d of %d IEDs significant using timing of beta power maxima models (Number of Connections)\n',sum(ptResults.sigBetaModel_numCons),length(IEDdata))
    fprintf('\n%d of %d IEDs significant using timing of beta power maxima models (Path Length)\n',sum(ptResults.sigBetaModel_pathLen),length(IEDdata))

    fprintf('\n%d of %d IEDs significant using timing of beta power maxima models (Euclidean Distance & Number of Connections)\n',sum(ptResults.sigBetaModel),length(IEDdata))
    fprintf('\n%d of %d IEDs significant using timing of beta power maxima models (Euclidean Distance & Path Length)\n',sum(ptResults.sigBetaModel & ptResults.sigBetaModel_numCons),length(IEDdata))
    fprintf('\n%d of %d IEDs significant using timing of beta power maxima models (Number of Connections & Path Length)\n',sum(ptResults.sigBetaModel_numCons & ptResults.sigBetaModel_pathLen),length(IEDdata))
    fprintf('\n%d of %d IEDs significant using timing of beta power maxima models (All three dsitance metrics)\n',sum(ptResults.sigBetaModel & ptResults.sigBetaModel_numCons & ptResults.sigBetaModel_pathLen),length(IEDdata))

    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

else
    fprintf('\nno results for patient %s',ptID)
end

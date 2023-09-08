
function[] = SOZChanTWCalc(ptID, SOZ)

% SOZChanTWCalc will calculate the number of occurences of an IED detection 
% in each channel across all time. Can be used to quantify the number of
% IEDs passing through the SOZ.
%
% inputs:   1) patient identifier
% inputs:   1) channels that define the SOZ

% Authors [AGY:20230426]

path_tw = sprintf('Channel_Order_TWs/%s_results.csv', ptID);
path_all = sprintf('Channel_Order_All/%s_results.csv', ptID);
load(fullfile('Data',ptID,'Imaging','Registered','Electrodes.mat'))
load(fullfile('Data',ptID,'Imaging','Registered','ChannelMap2.mat'))

numChans = length(ElecXYZProj(:,1));
ChanDistance = zeros(numChans, numChans);

chanFrequency = zeros(length(ElecMapRaw(:,3)), 6);
chanFrequency(:,1) = cell2mat(ElecMapRaw(:,3));

[~, loc] = ismember(SOZ,chanFrequency(:,1));
chanFrequency(loc,5) = 1;

for x = 1:numChans
    for y = 1:numChans
        ChanDistance(x,y) = sqrt(sum((ElecXYZProj(x,:) - ElecXYZProj(y,:)).^2,2)); 
    end
end

[loc,~] = find(ChanDistance(:,SOZ)<20);
SOZRegion = unique(loc);
[~, loc] = ismember(SOZRegion,chanFrequency(:,1));
loc = loc(loc>0);
chanFrequency(loc,6) = 1;

T = csvread(path_tw);
T_dist = T;

for x = 1:length(T(:,1))
    for y = 1:length(T(1,:))
        if(T(x,y) ~= 0)
            val = find(chanFrequency(:,1) == T(x,y));
            chanFrequency(val, 2) = chanFrequency(val, 2)+1;
        end
    end
end

T = csvread(path_all);
T_dist = T;

for x = 1:length(T(:,1))
    for y = 1:length(T(1,:))
        if(T(x,y) ~= 0)
            val = find(chanFrequency(:,1) == T(x,y));
            chanFrequency(val, 3) = chanFrequency(val, 3)+1;
        end
    end
end

chanFrequency(:,4) = chanFrequency(:,2)./chanFrequency(:,3);

path = sprintf('%s_Channel_Frequencies.csv', ptID);
writematrix(chanFrequency,path);
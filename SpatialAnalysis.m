
function[] = SpatialAnalysis(ptID, SOZ)

% SOZChanTWCalc returns the frequency in which each triplet (of electrodes) 
% detected an IED passing within the SOZ (row 4) vs. outside the SOZ (row 5)
%
% inputs:   1) patient identifier
% inputs:   2) channels that define the SOZ

% Authors [AGY:20230701]

path = sprintf('%s_results.csv', ptID);
load(fullfile('Data',ptID,'Imaging','Registered','Electrodes.mat'));
load(fullfile('Data',ptID,'Imaging','Registered','ChannelMap2.mat'));

numChans = length(ElecXYZProj(:,1));

T = csvread(path);
T_dist = T;

T_vals = zeros(nchoosek(numChans,3), 5);

for z = 1:length(T(:,1))
    row = T(z,:);
    if(length(intersect(row, SOZ)) > 0)
        x = 3;
        while(x <= length(row) && T(z,x) ~= 0)
            if(ismember(row(x-2:x),T_vals(:,1:3), 'rows'))
                index = find(T_vals(:,1:3) == row(x-2:x), 1, "first");
                T_vals(index,4) = T_vals(index,4)+1;
            else
                for y = 1:length(T_vals(:,1))
                    if(T_vals(y,1) ~= 0)
                        continue
                    end
                    T_vals(y,1:3) = row(x-2:x);
                    T_vals(y,4) = 1;
                    break;
                end
            end
            x = x+1;
        end
        
    else
        x = 3;
        while(x <= length(row) && T(z,x) ~= 0)
            if(ismember(row(x-2:x),T_vals(:,1:3), 'rows'))
                index = find(T_vals(:,1:3) == row(x-2:x), 1, "first");
                T_vals(index,4) = T_vals(index,5)+1;
            else
                for y = 1:length(T_vals(:,1))
                    if(T_vals(y,1) ~= 0)
                        continue
                    end
                    T_vals(y,1:3) = row(x-2:x);
                    T_vals(y,5) = 1;
                    break;
                end
            end
            x = x+1;
        end
    end
end

T_vals = T_vals(~all(T_vals == 0, 2),:);

filename = sprintf('%s_Triplet_Frequency.csv',ptID);
csvwrite(filename,T_vals);

                



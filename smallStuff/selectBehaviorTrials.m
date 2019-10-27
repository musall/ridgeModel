function bhv = selectBehaviorTrials(bhv,trials)
% Function to select a subset of trials/settings from selected 'trials' in a 
% larger array 'bhv' that has behavioral data. 'trials' should be a vector of 
% trial numbers that can be used. 
% Usage: bhv = selectBehaviorTrials(bhv,trials)

%% get fieldnames
if isempty(bhv)
    bFields = {};
else
    bFields = fieldnames(bhv);
end
bhv.nTrials = sum(bhv.nTrials);

%% check if trials is logical index. If not create index that matches trials in bhv.
if ~islogical(trials)
    temp = false(1,bhv.nTrials);
    temp(trials) = true;
    
    if length(temp) ~= bhv.nTrials
        warning('Trial index is larger as available trials in behavioral dataset')
        trials = temp(1:bhv.nTrials);
    else
        trials = temp;
    end
else
    if length(trials) ~= bhv.nTrials
        warning('Trial index has different length as available trials in behavioral dataset')
        trials = trials(1:bhv.nTrials);
    end
end
bhv.nTrials = sum(trials);

%% cycle trough fields and carry over selected trials / sessions
for iFields = 1:size(bFields,1)
    if ~any(ismember(size(bhv.(bFields{iFields})), length(trials))) %if field does not contain single trials, it should contain session data instead
        if isstruct(bhv.(bFields{iFields})) %if field is a struct, check one layer deeper if it contains trial info
            tFields = fieldnames(bhv.(bFields{iFields}));
            if length(bhv.(bFields{iFields}).(tFields{1})) == length(trials)
                bhv.(bFields{iFields}).(tFields{1}) = bhv.(bFields{iFields}).(tFields{1})(trials);
            else
                bhv.(bFields{iFields}) =  bhv.(bFields{iFields}); %carry over complete field
            end
        else
            bhv.(bFields{iFields}) =  bhv.(bFields{iFields}); %carry over complete field
        end
    else
        if isvector(bhv.(bFields{iFields}))
            bhv.(bFields{iFields}) = bhv.(bFields{iFields})(trials); %carry over selected trials
        else %some highD matrix, find trial dimension, cut trials and reshape to match original matrix
            cIdx = find(ismember(size(bhv.(bFields{iFields})), length(trials))); %find trial dimension
            cSize = size(bhv.(bFields{iFields}));
            cSize(cIdx) = sum(trials);
            temp = reshape(bhv.(bFields{iFields}), sum(cSize(1:cIdx-1)), cSize(cIdx), []);
            if cIdx == 1
                bhv.(bFields{iFields}) = reshape(temp(trials,:),cSize);
            else
                bhv.(bFields{iFields}) = reshape(temp(:,trials,:),cSize);
            end
        end
    end 
end

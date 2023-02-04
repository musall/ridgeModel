function bhv = selectBehaviorTrials(bhv, trials, nTrials)
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

if isfield(bhv, 'nTrials')
    bhv.nTrials = sum(bhv.nTrials);
    nTrials = bhv.nTrials;
end

%% check if trials is logical index. If not create index that matches trials in bhv.
if ~islogical(trials)
    temp = false(1,nTrials);
    temp(trials) = true;
    
    if length(temp) ~= nTrials
        warning('Trial index is larger as available trials in behavioral dataset')
        trials = temp(1:nTrials);
    else
        trials = temp;
    end
else
    if length(trials) ~= nTrials
        warning('Trial index has different length as available trials in behavioral dataset')
        trials = trials(1:nTrials);
    end
end

if isfield(bhv, 'nTrials')
    bhv.nTrials = sum(trials);
end

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
            bhv.(bFields{iFields}) = arrayIndex(bhv.(bFields{iFields}), trials, cIdx); %get index from target dimension
            
        end
    end 
end

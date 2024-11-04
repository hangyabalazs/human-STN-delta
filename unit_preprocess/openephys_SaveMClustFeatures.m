function openephys_SaveMClustFeatures(TT_file_name,FeaturesToUse,ChannelValidity,TT_path_name)
%OPENEPHYS_SAVEMCLUSTFEATURES   Save feature files for MClust import.
%   OPENEPHYS_SAVEMCLUSTFEATURES(TT,FEATURESTOUSE,CHANNELVALIDITY,TT_PATH)
%   calculates and saves feature files for MClust plotting. Input
%   arguments:
%       TT - TT file name with full path or TT data; in the latter case the
%           TT_PATH argument is required
%       FEATURESTOUSE - N-by-1 cell array of feature names to calculate
%       CHANNELVALIDITY - 1-by-4 array of 0 or 1 to indicate tetrode
%           channel validity
%       TT_PATH - full path of the TT data file; required in case of numerical
%           first input argument
%
%   See also CALCULATEFEATURES and WRITE_FD_FILE.

%   Balazs Hangya
%   Laboratory of Systems Neuroscience
%   Institute of Experimental Medicine, Budapest, Hungary

% Load data
if ischar(TT_file_name)   % works both with data and file path name
    TTdata = load(TT_file_name);  % if pathname was provided, load data
    TT_path_name = TT_file_name;
else
    TTdata = TT_file_name;  % data was provided as first input argument
end

% Calculate feature data
NumFeatures = length(FeaturesToUse);   % number of features to calculate
for iF = 1:NumFeatures
    [FeatureData, FeatureNames, FeaturePar] = feval(['feature_', FeaturesToUse{iF}],...
        tsd(TTdata.TimeStamps,TTdata.WaveForms), ChannelValidity);   %#ok<ASGLU> % use MClust feature code
    if isequal(FeaturesToUse{iF},'Time')
        FeatureData = FeatureData(:);  % 'Time' is returned with flipped dimensions
    end
    
    if isempty(FeaturePar)
        FeaturePar = 'empty'; %#ok<NASGU> % no feature parameters
    end
    
    FD_av = mean(FeatureData); %#ok<NASGU>  % mean feature data
    FD_sd = std(FeatureData)+eps; %#ok<NASGU>  % feature data SD
    
    FeatureTimestamps = TTdata.TimeStamps(:) * 1e+4;   % MClust time stamp convention
    FeatureIndex = 1:length(FeatureTimestamps); %#ok<NASGU>
    
    [pth, fnm] = fileparts(TT_path_name);
    FDfname = fullfile(pth,[fnm '_' FeaturesToUse{iF} '.fd']);   % feature file name
    save(FDfname, 'FeatureIndex','FeatureTimestamps','FeatureData', 'FeaturesToUse', 'ChannelValidity', 'FeatureNames', ...
        'FeaturePar','FD_av','FD_sd', 'TT_path_name', '-mat');
    disp([  ' Wrote ' FDfname ' as a .mat formatted file']);
end
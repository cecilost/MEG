%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocessing (MEG: Cam-CAN dataset)
% uses: camcan raw fif data files
% settings: 1-45 Hz lowpass, 48-52 Hz bandstop, 1000->100 Hz downsampling,
% 2 sepochs
% outlier criteria: +-three scaled median absolute deviations
% results saved as .mat files + images for checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Filtering:
% 
% Frequency Range: The data undergoes bandpass filtering between 1 Hz and 45 Hz. 
% This means that only the frequencies within this range are retained, and frequencies outside this range are attenuated.
% 
% Filter Type: A second-order Butterworth filter should be used? 
% The Butterworth filter is a type of filter that provides a flat frequency response in the passband.
% 
% Notch Filtering:
% 
% Stop Band: A notch filter is applied to remove a specific frequency component from the data. 
% In this case, the stop band is between 48 Hz and 52 Hz. This is typically done to eliminate power line noise, which often appears at 50 or 60 Hz depending on the region.
% 
% Filter Type: The specific type of notch filter (e.g., Butterworth, Chebyshev) is not mentioned, 
% but notch filters are commonly used to remove narrow-band noise.
% 
% Down-Sampling:
% 
% Sampling Rate Reduction: The data is down-sampled from its original sampling rate to 100 Hz. 
% This reduces the number of samples per second while retaining the essential information within the new frequency range determined by the bandpass filter.
%
% It is stated that: 1000->100 Hz downsampling
% % But in the code: %% downsample
%         cfg = [];
%         cfg.resamplefs = 200;% It sets the target sampling frequency (resamplefs) to 200 Hz
%         cfg.method = 'downsample'; % It specifies the downsampling method as 'downsample'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% cd 'C:/Users/stein/OneDrive/Escritorio/norm_conn-master/camcan/meg_prep/' % change directory to script folder
addpath 'C:/Program Files/fieldtrip-20230118'
ft_defaults
addpath('C:/Users/stein/OneDrive/Escritorio/norm_conn-master/camcan');
addpath('C:/Users/stein/OneDrive/Escritorio/norm_conn-master/utils/');
addpath('Matlab_FHZ_latex_library/');
addpath('Matlab_FHZ_latex_library/latex_library');
    
% fid1 = fopen('/datasabzi/results/CamCAN_new/preproc_errors.txt', 'w'); % open or create a file named 'preproc_errors.txt' in the specified directory.
fid1 = fopen('C:/Users/stein/OneDrive/Escritorio/norm_conn-master/results/CamCAN_new/preproc_errors.txt', 'w'); % open or create a file named 'preproc_errors.txt' in the specified directory

setup
% run('C:/Users/stein/OneDrive/Escritorio/norm_conn-master/camcan/setup.m'); % Use this if setup throws an error

%% Identify the subject you want to process
target_subject = 'sub-CC110069'; % One subject for this example

% For more than one subject, for example, create a list of subjects
% subject_list = {'sub-CC110069', 'sub-CC110070', 'sub-CC110071', ...}; % Add all subjects   

for nsub = 1%:length(target_subject)
% for nsub = 1:length(sub_info) % This loop iterates over each subject in the sub_info variable (631 x 3 double)
% for nsub = 1:length(sub_list) % This loop iterates over each subject in the subject_list   
    tic % starts a stopwatch timer. The 'toc' command is used later to measure how much time has passed during the execution of the script.
    clear info % clears the variable info. If info was previously defined or used in the script, this ensures a fresh start for each subject.
    
    %% Load the data for one subject
    % Set subj directly to target_subject
    subj = target_subject; 
    % %% load the data for one subject
    cfg = []; % Variable Initialization: A structure variable cfg is initialized, for configuration settings.
    % subj = ['sub-CC' num2str(sub_info(nsub))]; % Variable: This creates a string variable subj containing the subject identifier, like 'sub-CCxxx'.
    info.subject_id = subj; % Structure Update: The subject_id field in the info structure is updated with the value of subj.
    data_dir = [meg_dir subj '/'];
    save_prep = [preprocessed_dir];
    save_report = [report_dir subj '/'];
    % check if folder exists and report folder hasn't been created yet
    if (isfolder(data_dir)) %&& (~isfolder(save_report))
        filename = [data_dir '/ses-rest/meg/' subj '_ses-rest_task-rest_proc-sss.fif']; % use processed rest
        % Configuration setting:
        cfg.dataset = filename; % for fif data % 'cfg.dataset' is set to the 'filename', specifying the input data file for further processing.
        % raw_data = ft_preprocessing(cfg); % this was commented out in the original script
        % Preprocessing Configuration (cfg) Settings:
        cfg.continuous='yes'; % data is continuous
        cfg.channel='all'; % all channels should be included
        cfg.demean='yes'; % data should be demeaned
        cfg.bpfilter = 'yes'; % bandpass filtering will be applied
        cfg.bpfreq = [1 45]; % for 100 hz resolution. Lowpass: bandpass frequency range (1 to 45 Hz in this case). It  means: only the frequencies within this range are retained
        % cfg.bpfiltord = 2; % Should we add this? 'cfg.bpfiltord' should specify the filter order: second-order Butterworth (as this filter is mentioned in the draft but was not specified here)???        
        cfg.bsfilter = 'yes'; % a bandstop filter will be applied
        cfg.bsfreq = [48 52]; % bandstop frequency range (48 to 52 Hz in this case)
        
        % Data preprocessing:
        data=ft_preprocessing(cfg); % applies the preprocessing steps to the data using the FieldTrip 'ft_preprocessing' function
        
        % update the info structure with information about the preprocessing:
        info.freq = cfg.bpfreq; % frequency range. 'cfg.bpfreq' contains the bandpass frequency range specified in the preprocessing
        info.bs = cfg.bsfreq; % this was commented out in the original script % It would assign the bandstop frequency range ('cfg.bsfreq') to 'info.bs'
        info.orig_sf = data.fsample; % original sampling frequency. 'data.fsample' contains the original sampling frequency of the preprocessed MEG data
        info.nchan = length(data.label); % number of channels. 'data.label' contains the labels of all channels in the preprocessed data
        info.minutes = (size(data.time{1,1},2))/data.fsample/60; % duration of the data in minutes
      
             
        %% layout % This part of the script creates layouts for different sensor types (magnetometers & planar gradiometers) 
        grad = ft_read_sens(filename, 'senstype', 'meg'); % reads the sensor information from the MEG data file (filename) using the ft_read_sens function and stores it in the grad variable
        cfg = []; % sets up a configuration (cfg) structure for layout preparation
        cfg.grad = grad; % the sensor information (grad) is specified
        cfg.projection = 'polar'; % a projection method ('polar') is specified
	    % cfg.channel = 'MEG*1'; % these are the magnetometers 
        % layout.mag = ft_prepare_layout(cfg); % ERROR: Unable to perform assignment because dot indexing is not supported for variables of this type.
        % cfg.channel = 'MEG*2'; % planar gradiometers
        % layout.grad1 = ft_prepare_layout(cfg); % ERROR
        % cfg.channel = 'MEG*3'; % planar gradiometers
        % layout.grad2 = ft_prepare_layout(cfg); % ERROR
        layout = {}; % Initialization: an empty cell array (layout) is initialized to store the layout information for different sensor types
	    channel_names = {'MEG*1','MEG*2','MEG*3'}; % 'MEG*1' are the magnetometers, 'MEG*2' & 'MEG*3' are planar gradiometers
	    names = {'mag','grad1','grad2'};
        for c = 1:3
            cfg.channel = channel_names{c}; % these are the magnetometers % that is what said the original script, BUT??: 'MEG*1' are the magnetometers, 'MEG*2' & 'MEG*3' are planar gradiometers
            layout{c} = ft_prepare_layout(cfg);
        end % It loops over three sensor types (magnetometers and two types of planar gradiometers). For each type, it sets the channels to be considered in the layout using 'cfg.channel'. The layout information is then prepared using 'ft_prepare_layout' and stored in the 'layout' cell array.
        % get chanlabels for things outside FT
        eogs = find(~cellfun(@isempty,regexp(data.label,'^EOG')));
        ecg = find(~cellfun(@isempty,regexp(data.label,'^ECG')));
        mag = find(~cellfun(@isempty,regexp(data.label,'MEG\d*1$')));
        grad1 = find(~cellfun(@isempty,regexp(data.label,'MEG\d*2$')));
        grad2 = find(~cellfun(@isempty,regexp(data.label,'MEG\d*3$')));
        chantypes = [mag, grad1, grad2]; % It creates a combined list of indices for magnetometer, grad1, and grad2 channels in the chantypes variable
       
        %% downsample
        cfg = [];
        cfg.resamplefs = 200; % It sets the target sampling frequency (resamplefs) to 200 Hz. % IT SHOULD BE = 100 INSTEAD?
        cfg.method = 'downsample'; % It specifies the downsampling method as 'downsample'
        data_resampled=ft_resampledata(cfg,data); % It applies the downsampling configuration to the MEG data using the ft_resampledata function         
        info.fs = cfg.resamplefs; % It stores the updated sampling frequency in the info structure for later reference
      
        %% plot channels after loading, filtering and downsampling % visualize the MEG data quality after preprocessing
        % plot_data(save_report, data, layout.mag, 'MEG*1', 'after_import'); % ERROR: Dot indexing is not supported for variables of this type
      for c= 1:1 % This loop is set to run once (1:1): it will iterate only for one channel type, the 1st one:'MEG*1'
          plot_data(save_report, data, layout{c}, channel_names{c}, ['after_import_' names{c}]); % in plot_data() The data were split into epochs of 2 second length with no overlap
      end 
       
      %% regress EOG and ECG data % regressing out the EOG (electrooculogram) and ECG (electrocardiogram) data from the resampled EEG data
        % This step is commonly performed to remove artifacts from EEG data caused by eye movements (EOG) and heartbeats (ECG)
        eeg_data = data_resampled.trial{1,1}'; % Extract EEG data from the resampled MEG data
        eog_data = data_resampled.trial{1,1}(eogs,:)'; % Extract EOG data from the resampled MEG data
        ecg_data = data_resampled.trial{1,1}(ecg,:)'; % Extract ECG data from the resampled MEG data
        eeg_data = eeg_data - eog_data*(eog_data\eeg_data); % Perform linear regression to remove the contribution of EOG from the EEG data
        eeg_data = eeg_data - ecg_data*(ecg_data\eeg_data); % Perform linear regression to remove the contribution of ECG from the EEG data
        data_resampled.trial{1,1} = eeg_data'; % Update the resampled MEG data with the regressed EEG data
        clearvars eeg_data eog_data ecg_data % Clear intermediate variables to free up memory
        info.eog = eogs; % Update the info structure with information about the EOG channels.
        info.ecg = ecg; % Update the info structure with information about the ECG channels.
       
        %% check after regression % checking the data quality after the regression step
        % plot_data(save_report, data, layout.mag, 'MEG*1', 'after_phys_regr'); % ERROR: Dot indexing is not supported for variables of this type
    for c = 1:1
        plot_data(save_report, data_resampled, layout{c}, channel_names{c}, ['after_regression_' names{c}]);
	end
        
    %% reject outlying channels % identify and interpole outlying channels in the MEG data for different sensor types
        % this part detects outlying channels for each sensor type 
        % (outliers: 3 scaled median deviations from the median)
        % and interpolates them using a weighted average of neighbors
        chanlist = 1:length(data.label);
        all_bad_channels = []; 
        % all_bad_ch = []; % I changed the two mentions to 'all_bad_ch' to 'all_bad_channels'. 1st mention here
        data_cleaned = data_resampled; % fieldtrip structure % HERE I CHANGED 'data' to 'data_resampled'; ths fixed the ERROR with 'trial'
        for ch_type = 1:size(chantypes,2) % iterates over channel types (mag and grad)
            % returns the interpolated data and a list of bad channels (bc):
            [data_interpolated{ch_type}, bc] = interpolate_bad_channels(save_report, chantypes(:,ch_type), data_resampled, layout{ch_type}, channel_names{ch_type}); % I added 'save_report', ex 'data' now 'data_resampled'
            % put into FieldTrip struct
            data_cleaned.trial{1,1}(chantypes(:,ch_type),:) = data_interpolated{ch_type}; % Updates the cleaned data with the interpolated data for the current channel type
            all_bad_channels = [all_bad_channels;bc]; % Appends the information about bad channels for the current type to the overall list of bad channels
        end

        %% I commented out the lines below because 'data_cleaned = data'; try this instead: 'data_cleaned = data_resampled' (as in preprocess_remote)
        % With the code below: Error using interpolate_bad_channels. The input data does not have the expected 'trial' field.
        % chanlist = 1:length(data.label);
        % all_bad_ch = [];
        % data_cleaned = data; % fieldtrip structure
        %% I commented out the lines below because 'channel_types' does not exist; I refer to 'chantypes' instead
        % for ch_type = 1:size(channel_types,2) % iterate over type (mag and grad)
        %     [data_interpolated{ch_type}, bc] = interpolate_bad_channels(channel_types(:,ch_type), data, layout{ch_type}, channel_names{ch_type});
        %     % put into FieldTrip struct
        %     data_cleaned.trial{1,1}(channel_types(:,ch_type),:) = data_interpolated{ch_type};
        %     all_bad_channels = [all_bad_channels;bc];
        % end
        
        clearvars bc
        info.badchans = all_bad_channels; % Assigns the variable all_bad_channels to the field badchans in the structure info
        % info.badchans = all_bad_ch; I changed the two mentions to 'all_bad_ch' to 'all_bad_channels'. 2nd mention here
        
        %% segment the data into 2s % and update the info structure with relevant information about the segmentation
        cfg = [];
        cfg.length = 2; % Configures the segment length to be 2 seconds
        cfg.overlap = 0; % Specifies that there is no overlap between consecutive segments
        data_seg = ft_redefinetrial(cfg,data_cleaned); % Redefines the trials of the cleaned data (data_cleaned) into non-overlapping segments of 2 seconds each. This function is part of the FieldTrip toolbox and is used for trial-based analysis
        % data_seg = ft_redefinetrial(cfg,interdata); ERROR: Unrecognized function or variable 'interdata'. % I changed 'interdata' to 'data_cleaned' 
        info.sam_length = cfg.length; % Stores the segment length in the info structure
        info.nsam = size(data_seg.time{1,1},2); % Retrieves the number of samples in each segment and stores it in the info structure
        info.ntri = size(data_seg.trial,2); % Retrieves the number of trials (segments) and stores it in the info structure
        
        %% reject outlier trials
        all_bad_tri = [];
        for chtype = 1:size(chantypes,2) % a loop iterating over different channel types (2 specifies that we are interested in the number of columns). Iterates over the range from 1 to the number of columns in chantypes.
            bt = detect_bad_trials(save_report, chantypes(:,chtype), data_seg); % calls a function detect_bad_trials to identify bad trials for the current channel type
            all_bad_tri = [all_bad_tri;bt]; %  appends the bad trials detected for the current channel type to the overall list of bad trials
        end
        all_bad_tri = unique(all_bad_tri); %  removes duplicate entries from the list of bad trials
        disp('bad trials:');
        fprintf(1, '%d \n', all_bad_tri);
        info.badtrials = all_bad_tri;
        info.prcbadtri = num2str((length(all_bad_tri)/length(data_seg.trial))*100); % Updates the info structure with information about bad trials and the percentage of trials that are considered bad
        % reject trials here and select only meg channels for saving
        cfg = [];
        cfg.trials = 1:length(data_seg.trial);
        cfg.trials(all_bad_tri) = [];
        cfg.channel = 'meg';
        data_seg = ft_selectdata(cfg,data_seg); % rejects the bad trials from the data and selects only MEG channels for further processing, using the FieldTrip ft_selectdata function
       
        %% plot again for report
        % plot_data(save_report, data_seg, layout.mag,'MEG*1', 'after_rej'); % ERROR
    for c= 1:1
        plot_data(save_report, data_seg, layout{c}, channel_names{c}, ['after_rej_' names{c}]);
	end
        timer = toc; % (it started with 'tic'). 'toc': to measure how much time has passed during the execution of the script.
        info.time = timer;
      
              
        %% save
        disp(['This data is saved in ' [save_prep subj] ' under the name' subj]); % displays a message in the command window
        ft_write_data([save_prep subj], data_seg, 'dataformat', 'matlab'); % uses the FieldTrip function ft_write_data to write the data (data_seg) to a file. The data is saved in MATLAB MAT-file format
        mkdir(save_report); % creates a directory with the name specified by save_report if it doesn't already exist
        save([save_report '/info.mat'],'info'); % saves the MATLAB variable info into a MAT-file named info.mat in the directory specified by save_report
        % mkdir(save_report); % this line should ideally come before the attempt to save the 'info.mat' file
        % preproc_report(info,save_report); % ERROR: Unrecognized function or variable 'preproc_report' 
        
        % Create a PDF report with subject identifier
        reportFileName = [save_report '/report_' subj '.pdf'];
        fidReport = fopen(reportFileName, 'w');

     if fidReport ~= -1
            % Open the PDF document for writing
            fprintf(fidReport, '\\documentclass{article}\n\\begin{document}\n');

            % Write information from the info struct to the PDF
            fprintf(fidReport, 'Report for subject %s:\n', subj);
            fprintf(fidReport, 'Subject: %s\n', info.subject_id); % this is the same info than in the line above
            fprintf(fidReport, 'Frequency range (bandpass): %s\n', num2str(info.freq));
            fprintf(fidReport, 'Frequency range (bandstop or notch filter): %s\n', num2str(info.bs));
            fprintf(fidReport, 'Original Sampling Frequency: %d\n', info.orig_sf);
            fprintf(fidReport, 'Number of Channels: %d\n', info.nchan);
            fprintf(fidReport, 'Duration of Data (in minutes): %.2f\n', info.minutes);
            fprintf(fidReport, 'Sampling frequency: %d\n', info.fs);
            fprintf(fidReport, 'EOG channels: %d, %d\n', info.eog);
            fprintf(fidReport, 'ECG channels: %d\n', info.ecg);
            fprintf(fidReport, 'Bad channels: %d\n', info.badchans);
            fprintf(fidReport, 'Segment length: %d\n', info.sam_length);
            fprintf(fidReport, 'Number of samples in each segment: %d\n', info.nsam);
            fprintf(fidReport, 'Number of trials (segments): %d\n', info.ntri);
            fprintf(fidReport, 'Bad trials: %d\n', info.badtrials);
            
            % Close the PDF document
            fprintf(fidReport, '\\end{document}');
            fclose(fidReport);

            % Convert tex to pdf
            command = ['cd ', save_report, '; ', 'pdflatex report.tex;'];
            system(command);

        % % convert tex to pdf
        % command = ['cd ', save_report,'; ', 'yes " " | /usr/bin/pdflatex  ',...
        % 'report.tex;']; % changes the current working directory to save_report and then executes a shell command (/usr/bin/pdflatex) to convert a LaTeX file (report.tex) to a PDF file
        % system(command);
    else
        fprintf(fid1,'%s \n',['No MEG data found for subject ' subj]);
        continue
     end
    end
end
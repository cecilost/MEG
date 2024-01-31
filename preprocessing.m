%add fieldtrip to path
%replace your fieldtrip path 
addpath /home/erfan/fieldtrip-20230118/
ft_defaults
clear info

report_dir = './camcan_results/';
preprocessed_dir = report_dir;
save_prep = preprocessed_dir;

%dir for the sample data
data_dir = './MEG_sample/camcan_sample/';
subjects = dir(data_dir);
subjects = subjects([subjects.isdir]);  % keep only directories
subjects = subjects(~ismember({subjects.name}, {'.', '..'}));  % remove '.' and '..'

%% loop through each subject 
for i = 1:1 %size(subjects) %in case for all the present subj
    subj = subjects(i).name;
     
    % search for .fif files in any subfolder
    fif_files = dir(fullfile(data_dir, subj, '**', '*.fif'));

    save_report = [report_dir subj '/'];
    if ~isfolder(save_report)
        mkdir(save_report)
    end
    % create the full path to the MEG data file for the current subject
    filename = [fif_files.folder '/' fif_files.name];
%% Read and plot raw data 
    % call ft_preprocessing and et the channel length--not used for further step-- 
    cfg=[];
    cfg.dataset= filename;  
    raw_data=ft_preprocessing(cfg);
    nchans = length(raw_data.label);
    
    %magnetometters layout
    grad = ft_read_sens(filename, 'senstype', 'meg');
    cfg.grad = grad;
    cfg.projection = 'polar';
    layout = {};
    channel_names = 'MEG*1';
    names = 'mag';
    
    cfg.channel = channel_names; % these are the magnetometers 
    layout = ft_prepare_layout(cfg);
    
    plot_data(save_report, raw_data, layout, channel_names, ['raw_after_import_' names]);
    
    clearvars raw_data layout
%% Read and apply filter channesl by channel
    % loading each channel at at a time (for memory constrains)
    % filteriing, and resampleing to 200
    for i=1:nchans
        
        cfgp= [];
        cfgp.dataset = filename;
        cfgp.channel = i;
        %baseline correction
        cfgp.demean='yes';
        cfgp.bpfilter = 'yes';
        % default bpfilter type is Butterworth 
        % from ft_preproc_bandpassfilter, low level preprocessin function:
        % filter order, default is 4 (but) or dependent on frequency band and data length (fir/firls)
        cfgp.bpfreq = [1 45]; %for 100 hz resolution
        cfgp.bsfilter = 'yes';
        cfgp.bsfreq = [48 52];
        datp= ft_preprocessing (cfgp);
        %lp hp ??
        
        % cfgp.lpfilter="yes";
        % cfgp.hpfilter="yes";
        % cfgp.lpfreq=48;
        % cfgp.hpfreq=52;

        %saving the info
        info.freq = cfgp.bpfreq;
        info.bs = cfgp.bsfreq;
        info.orig_sf = datp.fsample;
        info.nchan = length(datp.label);
        info.minutes = (size(datp.time{1,1},2))/datp.fsample/60;
        
        %Downsampleing 
             
        cfgr= [];
        cfgr.resamplefs = 200;
        cfgr.method = 'downsample';
        info.fs = cfgr.resamplefs;
        datr{i}= ft_resampledata(cfgr, datp);
        
        clearvars datp
    end
    
    % apend the channels
    cfg = [];
    resampled_filtered_data = ft_appenddata(cfg, datr{:});
    %% layout
    % https://www.fieldtriptoolbox.org/faq/why_are_there_multiple_neighbour_templates_for_the_neuromag306_system/
    % https://imaging.mrc-cbu.cam.ac.uk/meg/VectorviewDescription#Magsgrads
    %  a combination of magnetometers and two planar gradiometers at each sensor location. planargradiometers (fT/cm) and magnetometeres (fT)
    % https://www.fieldtriptoolbox.org/template/layout/#neuromagelektamegin-system
    
    grad = ft_read_sens(filename, 'senstype', 'meg');
    cfg = [];
    cfg.grad = grad;
    cfg.projection = 'polar';
    layout = {};
    channel_names = {'MEG*1','MEG*2','MEG*3'};
    names = {'mag','grad1','grad2'};
    
    for c = 1:3
        cfg.channel = channel_names{c}; % these are the magnetometers 
        %https://github.com/fieldtrip/fieldtrip/blob/master/ft_prepare_layout.m
        layout{c} = ft_prepare_layout(cfg);
    end
        %get chanlabels for things outside FT
        eogs = find(~cellfun(@isempty,regexp(resampled_filtered_data.label,'^EOG')));
        ecg = find(~cellfun(@isempty,regexp(resampled_filtered_data.label,'^ECG')));
        mag = find(~cellfun(@isempty,regexp(resampled_filtered_data.label,'MEG\d*1$')));
        grad1 = find(~cellfun(@isempty,regexp(resampled_filtered_data.label,'MEG\d*2$')));
        grad2 = find(~cellfun(@isempty,regexp(resampled_filtered_data.label,'MEG\d*3$')));
        channel_types = [mag, grad1, grad2];
      
    for c= 1:1
            plot_data(save_report, resampled_filtered_data, layout{c}, channel_names{c}, ['after_import_' names{c}]);
    end
    
    %% regress EOG and ECG data
    eeg_data = resampled_filtered_data.trial{1,1}';
    eog_data = resampled_filtered_data.trial{1,1}(eogs,:)';
    ecg_data = resampled_filtered_data.trial{1,1}(ecg,:)';
    eeg_data = eeg_data - eog_data*(eog_data\eeg_data);
    eeg_data = eeg_data - ecg_data*(ecg_data\eeg_data);
    resampled_filtered_data.trial{1,1} = eeg_data';
    clearvars eeg_data eog_data ecg_data
    
    info.eog = eogs;
    info.ecg = ecg;
    
    for c = 1:1
        plot_data(save_report, resampled_filtered_data, layout{c}, channel_names{c}, ['after_regression_' names{c}]);
    end
    
    chanlist = 1:length(resampled_filtered_data.label);
    all_bad_channels = [];
    data_cleaned = resampled_filtered_data; % fieldtrip structure
    for ch_type = 1:size(channel_types,2) % iterate over type (mag and grad)
        [data_interpolated{ch_type}, bc] = interpolate_bad_channels(save_report, channel_types(:,ch_type), resampled_filtered_data, layout{ch_type}, channel_names{ch_type});
        % put into FieldTrip struct
        data_cleaned.trial{1,1}(channel_types(:,ch_type),:) = data_interpolated{ch_type};
        all_bad_channels = [all_bad_channels;bc];
    end
    clearvars bc
    info.badchans = all_bad_channels;
    
    %% segment the data into 2s
    cfg = [];
    cfg.length = 2;
    cfg.overlap = 0;
    data_seg = ft_redefinetrial(cfg,data_cleaned);
    info.sam_length = cfg.length;
    info.nsam = size(data_seg.time{1,1},2);
    info.ntri = size(data_seg.trial,2);
    %% reject outlier trials
    all_bad_trials = [];
    for ch_type = 1:size(channel_types,2)
        bt = detect_bad_trials(save_report, channel_types(:,ch_type), data_seg);
        all_bad_trials = [all_bad_trials;bt];
    end
    all_bad_tri = unique(all_bad_trials);
    disp('bad trials:');
    fprintf(1, '%d \n', all_bad_trials);
    info.badtrials = all_bad_tri;
    info.prcbadtri = num2str((length(all_bad_tri)/length(data_seg.trial))*100);
    % reject trials here and select only meg channels for saving
    cfg = [];
    cfg.trials = 1:length(data_seg.trial);
    cfg.trials(all_bad_tri) = [];
    cfg.channel = 'meg';
    data_seg = ft_selectdata(cfg,data_seg);
    %% plot again for report
    for c= 1:1
        plot_data(save_report, data_seg, layout{c}, channel_names{c}, ['after_rej_' names{c}]);
    end
    timer = toc;
    info.time = timer;
     %% save
    info.subject_id= subj;
    disp(['This data is saved in ' [save_prep subj] ' under the name' subj]);
    ft_write_data([save_prep subj], data_seg, 'dataformat', 'matlab');
    save([save_report '/info.mat'],'info');
    % mkdir(save_report);
    preproc_report(info,save_report);
    % convert tex to pdf
    command = ['cd ', save_report,'; ', 'yes " " | /usr/bin/pdflatex  ',...
    'report.tex;'];
    system(command);

    clearvars resampled_filtered_data data_seg data_cleaned
end

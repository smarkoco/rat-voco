%% clear and close figures, load data, and initialize key variables
clear; close all
audio_filepath_parent = "C:\Users\Public\Documents\Avisoft Bioacoustics\data\";
disp("Choose an audio recording.")
[audio_filename, audio_filepath_parent]= uigetfile(strcat(audio_filepath_parent,"*wav"));
dbstop if error
% read audio file
full_audio_filepath = strcat(audio_filepath_parent, audio_filename);
int16_wavdata = audioread(full_audio_filepath,'native');
audio_fs = 250000; % Hz, sampling frequency
time_axis_for_audio_sec = (0:1/audio_fs:(length(int16_wavdata)-1)/audio_fs)'; % make time axis a column array

% convert to binary and extract LSB (contains SYNC)
binary_wavdata = dec2bin(int16_wavdata,16);
LSB_array = binary_wavdata(:,16)=='1';
TTL_rising_edge = find(LSB_array,1,'first');
TTL_falling_edge = find(LSB_array,1,'last');
Avisoft_data = double(int16_wavdata);

% read and create time axis for ephys data
disp("Now choose the corresponding Intan file.")
Intan_parent_folder = "C:\Users\Public\Documents\Intan\";
read_Intan_RHD2000_file(Intan_parent_folder)
ephys_channel_to_view = 2;
ephys_fs = 30000; % Hz, sampling frequency
Intan_amplifier_data = amplifier_data';
time_axis_for_ephys_sec = 0:1/ephys_fs:(length(Intan_amplifier_data(:, ephys_channel_to_view))-1)/ephys_fs;

% align audio and ephys data, accounting for sampling rate differences
TTL_rising_edge_in_time = TTL_rising_edge/audio_fs;
TTL_falling_edge_in_time = TTL_falling_edge/audio_fs;
time_difference_of_last_point = time_axis_for_audio_sec(end) - TTL_falling_edge_in_time;
% shifted_time_axis_for_audio_sec = time_axis_for_audio_sec-TTL_rising_edge_in_time;
SYNC_received_on_Avisoft_mV = double(LSB_array)*5000;
SYNC_received_on_Intan_mV = (board_dig_in_data*5000)';
pretime = 1.0;
posttime = 1.0;

%% plot Intan SYNC pulses recorded on both Intan and Avisoft
% figure(1); hold on
% plot(time_axis_for_audio_sec,LSB_array)
% plot(time_axis_for_ephys_aligned_to_TTL_rise,board_dig_in_data)
% title('TTL "Mark Out" Signal from Intan RHD Controller')
% xlabel('Time (s)')
% legend(["250kHz Audio Data","30kHz Ephys Data"])
% hold off

% count number of pulses and stretch Audio data
% by the ratio of time passed for ephys and audio
% (using the times between 1st and Nth pulse)
[~,all_Intan_pulses,~,~] = findpeaks(SYNC_received_on_Intan_mV,'MinPeakWidth',ephys_fs*(1/100/3)); % divide by 3 to ensure less than half-cycle time of arduino square (100Hz)
[~,all_Avisoft_pulses,~,~] = findpeaks(SYNC_received_on_Avisoft_mV,'MinPeakWidth',audio_fs*(1/100/3));
% [~,all_big_Intan_pulses,~,~] = findpeaks(SYNC_received_on_Intan_mV,'MinPeakWidth',ephys_fs*(1/100)); % wide pulses mark 1 second
% [~,all_big_Avisoft_pulses,~,~] = findpeaks(SYNC_received_on_Avisoft_mV,'MinPeakWidth',audio_fs*(1/100));
[~,all_Avisoft_falling_edge_idxs,~,~] = findpeaks(-SYNC_received_on_Avisoft_mV);
idx_of_Avisoft_falling_edges_which_are_wide = find(diff(all_Avisoft_falling_edge_idxs)>0.016*audio_fs)+1;
all_Avisoft_big_falling_edge_idxs = all_Avisoft_falling_edge_idxs(idx_of_Avisoft_falling_edges_which_are_wide);
[~,all_Intan_falling_edge_idxs,~,~] = findpeaks(-SYNC_received_on_Intan_mV);
idx_of_Intan_falling_edges_which_are_wide = find(diff(all_Intan_falling_edge_idxs)>0.016*ephys_fs)+1;
all_Intan_big_falling_edge_idxs = all_Intan_falling_edge_idxs(idx_of_Intan_falling_edges_which_are_wide);

% % validate Avisoft big pulse findings
% figure(2)
% plot(SYNC_received_on_Avisoft_mV); hold on
% scatter(all_Avisoft_big_falling_edge_idxs,ones(size(all_Avisoft_big_falling_edge_idxs))*5000,'ro'); hold off
% % validate Intan big pulse findings
% figure(3)
% plot(SYNC_received_on_Intan_mV); hold on
% scatter(all_Intan_big_falling_edge_idxs,ones(size(all_Intan_big_falling_edge_idxs))*5000,'ro'); hold off

%% alignment and stretching of Avisoft (Intan serves as the reference)
time_axis_for_audio_aligned_to_TTL_rise = time_axis_for_audio_sec - (all_Avisoft_pulses(1)-1)/audio_fs;
time_axis_for_ephys_aligned_to_TTL_rise = time_axis_for_ephys_sec - (all_Intan_pulses(1)-1)/ephys_fs;
time_axis_for_ephys_aligned_to_TTL_fall = time_axis_for_ephys_sec + time_difference_of_last_point;
ephys_duration = length(time_axis_for_ephys_aligned_to_TTL_rise)/ephys_fs;
num_Intan_pulses = length(all_Intan_pulses);
num_Avisoft_pulses = length(all_Avisoft_pulses);
index_of_1st_Avisoft_pulse = all_Avisoft_pulses(1);
index_of_1st_Intan_pulse_in_Avisoft_rate = round(all_Intan_pulses(1)*audio_fs/ephys_fs);
N = num_Intan_pulses; % set N=num_Intan_pulses for most accurate result
index_of_Nth_Avisoft_pulse = all_Avisoft_pulses(N);
index_of_Nth_Intan_pulse_in_Avisoft_rate = round(all_Intan_pulses(N)*audio_fs/ephys_fs);
factor_to_stretch_Avisoft_by = (index_of_Nth_Intan_pulse_in_Avisoft_rate-index_of_1st_Intan_pulse_in_Avisoft_rate)/(index_of_Nth_Avisoft_pulse-index_of_1st_Avisoft_pulse);
stretched_time_axis_for_audio = (time_axis_for_audio_aligned_to_TTL_rise)*factor_to_stretch_Avisoft_by;
relative_audio_fs = audio_fs*factor_to_stretch_Avisoft_by;
factor_Avisoft_stretched_by = factor_to_stretch_Avisoft_by;

%% plot aligned data
% figure(4); hold on
% plot(time_axis_for_audio_aligned_to_TTL_rise,SYNC_received_on_Avisoft_mV,'-+')
% plot(time_axis_for_ephys_aligned_to_TTL_rise,SYNC_received_on_Intan_mV,'-+')
% plot(time_axis_for_audio_aligned_to_TTL_rise,int16_wavdata)
% scatter(time_axis_for_ephys_aligned_to_TTL_rise, Intan_amplifier_data(ephys_channel_to_view,:),2,'k.')
% plot(-pretime:pretime:0, [5000, 5000], 'r', 'LineWidth', 2); % check that alignment matches the expected 1.0 sec pre-trigger time window
% plot(time_axis_for_audio_aligned_to_TTL_rise(end)-posttime:posttime:time_axis_for_audio_aligned_to_TTL_rise(end), [5000, 5000], 'r', 'LineWidth', 2) % check expected posttime seconds "hold time" window
% title('Aligned but Unstretched Audio and Ephys Data')
% xlabel('Time (s)')
% legend(["SYNC recv'd on Avisoft","SYNC recv'd on Intan","250kHz Audio Data","30kHz Ephys Data","Pre-rec Time","Post-rec Time"])
% hold off

%%% segment the data into minute-long chunks, and store in cell array, for
%%% saving later
% % audio data
% rounded_minutes_of_audio_data = round(length(stretched_time_axis_for_audio)/relative_audio_sample_rate/60);
% whole_minutes_containing_audio_data = ceil(length(stretched_time_axis_for_audio)/relative_audio_sample_rate/60);
% audio_data_cell_array = cell(whole_minutes_containing_audio_data,1);
% audio_file_boundaries = 1:relative_audio_sample_rate:whole_minutes_containing_audio_data*relative_audio_sample_rate;
% % ephys data
% whole_minutes_containing_ephys_data = ceil(length(time_axis_for_ephys_aligned_to_TTL_rise)/ephys_fs/60);
% ephys_data_cell_array = cell(whole_minutes_containing_ephys_data,1);

%% plot aligned and stretched data
% figure(5); hold on
% plot(stretched_time_axis_for_audio, SYNC_received_on_Avisoft_mV,'-+')
% plot(time_axis_for_ephys_aligned_to_TTL_rise, SYNC_received_on_Intan_mV,'-+')
% plot(stretched_time_axis_for_audio,int16_wavdata)
% scatter(time_axis_for_ephys_aligned_to_TTL_rise, Intan_amplifier_data(ephys_channel_to_view,:),2,'k.')
% plot(-pretime:pretime:0, [5000, 5000], 'r', 'LineWidth', 2); % check that alignment matches the expected 1.0 sec pre-trigger time window
% plot(stretched_time_axis_for_audio(end)-posttime:posttime:stretched_time_axis_for_audio(end), [5000, 5000], 'r', 'LineWidth', 2) % check expected 1.0 sec "hold time" window
% scatter(time_axis_for_ephys_aligned_to_TTL_rise(all_Intan_pulses(100:99:end)),5000*ones(size(all_Intan_pulses(100:99:end)))) % plot 100th pulses
% title('Aligned and Stretched Audio and Ephys Data')
% xlabel('Time (s)')
% legend(["SYNC recv'd on Avisoft","SYNC recv'd on Intan","250kHz Audio Data","30kHz Ephys Data","Pre-rec Time","Post-rec Time","Every 100th TTL"])
% hold off

% histogram
% histogram(diff(all_Intan_pulses/ephys_fs),'facecolor','k','binwidth',1/relative_audio_fs); hold on
% histogram(diff(all_Avisoft_pulses/relative_audio_fs),'facecolor','r','binwidth',1/relative_audio_fs)
% title('Comparison of Inter-Pulse Timing for SYNC Across Intan and Avisoft Devices')
% xlabel('Time (s)')
% ylabel('Pulse Count')
% legend(["Intan SYNC","Avisoft SYNC"])

%% chop data into 1 minute chunks; also FYI, there is an extra 5ms at beginning due to Arduino timer starting HIGH for half a 100Hz cycle
ephys_data_1sec_delimiters = [(all_Intan_pulses(1)-1); all_Intan_big_falling_edge_idxs];
audio_data_1sec_delimiters = [(all_Avisoft_pulses(1)-1); round(all_Avisoft_big_falling_edge_idxs*factor_Avisoft_stretched_by)];
ephys_data_1min_delimiters = ephys_data_1sec_delimiters(1:60:end);
audio_data_1min_delimiters = audio_data_1sec_delimiters(1:60:end);
num_chunks = ceil(length(all_Intan_big_falling_edge_idxs)/60); % minute long chunks
% preallocate variables
sessionID = string(audio_filename).split(".");
sessionID = sessionID(1);
synced_struct(num_chunks) = struct();
% load chunked data into appropriate struct locations
for iChunk=1:num_chunks
    ephys_start_idx = ephys_data_1min_delimiters(iChunk);
    audio_start_idx = audio_data_1min_delimiters(iChunk);
    if iChunk ~= num_chunks
        ephys_end_idx = ephys_data_1min_delimiters(iChunk+1);
        audio_end_idx = audio_data_1min_delimiters(iChunk+1);
    else
        ephys_end_idx = ephys_data_1sec_delimiters(end);
        audio_end_idx = audio_data_1sec_delimiters(end);
    end
    synced_struct(iChunk).ephys.time = time_axis_for_ephys_aligned_to_TTL_rise(ephys_start_idx:ephys_end_idx)';
    synced_struct(iChunk).ephys.data = Intan_amplifier_data(ephys_start_idx:ephys_end_idx,:);
    synced_struct(iChunk).ephys.sync = SYNC_received_on_Intan_mV(ephys_start_idx:ephys_end_idx);
    synced_struct(iChunk).audio.time = stretched_time_axis_for_audio(audio_start_idx:audio_end_idx);
    synced_struct(iChunk).audio.data = int16_wavdata(audio_start_idx:audio_end_idx);
    synced_struct(iChunk).audio.sync = SYNC_received_on_Avisoft_mV(audio_start_idx:audio_end_idx);
    synced_struct(iChunk).info.factor_Avisoft_stretched_by = factor_Avisoft_stretched_by;
    synced_struct(iChunk).info.sessionID = sessionID;
    synced_struct(iChunk).info.chunkID = iChunk;
end

%% ask user whether to save data into .mat files, does not save by default
% TODO: add saving modes, to store data by full blocks or minute chunks
user_answer = input("Save audio and SYNC data into .mat? (y/[n])",'s');
if any(char(user_answer) == 'y') || any(char(user_answer) == 'Y')
    save_mat_path_arr = split(full_audio_filepath,'.');
    save_mat_path = save_mat_path_arr{1};
%     save(save_mat_path, ...
%         "stretched_time_axis_for_audio","Avisoft_data","SYNC_received_on_Avisoft_mV","factor_Avisoft_stretched_by", ...
%         "time_axis_for_ephys_aligned_to_TTL_rise","Intan_amplifier_data","SYNC_received_on_Intan_mV", '-v7.3')
    for iChunk=1:num_chunks
        iChunk_str = string(iChunk);
        synced = synced_struct(iChunk);
        save_mat_path_chunk = strcat(save_mat_path,"-",iChunk_str,".mat");
        save(save_mat_path_chunk, "synced", "-v7.3")
    end
    parts = strsplit(path, filesep);
    parent_path = strjoin(parts(1:end-1), filesep);
    fprintf("Saved files to %s\n", parent_path)
elseif any(isempty(user_answer))
    disp("Did not save any data arrays into .mat files.")
else
    disp("Did not save any data arrays into .mat files.")
end

%% validate synced_struct chunking by plotting sequentially
f6 = figure(6); hold on
title("Synchonized, Chunked Data (written to disk as separate files)")
xlabel("Time (s)")
ylabel("Various signal amplitudes")
good_chan_idxs = [1,2];%,4,6,7,8,14,15];
disp("Using channels: ")
disp(good_chan_idxs)
disp("Number of chunks in this recording:")
disp(num_chunks)

for iChunk=1:num_chunks
    disp(strcat("chunk: ",string(iChunk)))
    user_answer2 = input("See this chunk? y/[n]","s");
    if any(char(user_answer2) == 'y') || any(char(user_answer2) == 'Y')
        plot(synced_struct(iChunk).audio.time, synced_struct(iChunk).audio.sync)
        plot(synced_struct(iChunk).ephys.time, synced_struct(iChunk).ephys.sync)
        plot(synced_struct(iChunk).audio.time, synced_struct(iChunk).audio.data)
        plot(synced_struct(iChunk).ephys.time, synced_struct(iChunk).ephys.data(:,good_chan_idxs))

        legend(["ephys sync","audio sync","ephys data","audio data"])
%     elseif any(isempty(user_answer))
%         disp("Quitting. Have a good day :)")
%         close(f6)
%         return
%     else
%         disp("Quitting. Have a good day :)")
%         close(f6)
%         return
%     end
    end
end
hold off
loaded_file = load('B1_230302121440-1');
audio = single(loaded_file.synced.audio.data);

% set spectral variables
fs = 250000;
points_per_segment = 2500; % 10ms bins
segment_overlap_pts = 1000; 
window = bartlett(points_per_segment);
Ndft = 4096; %max(256,2^nextpow2(points_per_segment));
disp("Ndft")
disp(Ndft)

[ywinhat,f,t,Pow] = spectrogram(audio,window,segment_overlap_pts,Ndft,fs);
S = Pow*fs/Ndft;
SdB = 10*log10(S); % Convert spectral power to decibel scale
figure(1)
% plot audio spectrogram
ax1 = subplot('Position',[.1 .5 .8 .4]);
imagesc(loaded_file.synced.audio.time,f/1000,SdB)
% uses the range of values of SdB to make the color scale.
ylabel("Frequency (kHz)")
axis xy % Puts low frequencies at the bottom
colormap parula
c1 = colorbar;
c1.Label.String = 'Power/Frequency (dB/Hz)';
caxis([-10 60])

ax2 = subplot('Position',[.1 .1 .8 .4]);
% ax2 = subplot();
% plot ephys underneath spectrogram
channels_to_show = [1];
plot(loaded_file.synced.ephys.time, loaded_file.synced.ephys.data(:,channels_to_show),'k')
xlabel('Time (s)')
ylabel('Voltage (uV)')
c2 = colorbar('Visible','off');
% axis linked
linkaxes([ax1,ax2],'x');
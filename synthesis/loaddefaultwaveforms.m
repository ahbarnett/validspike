function wf = loaddefaultwaveforms
% LOADDEFAULTWAVEFORMS - load a default waveform object for testing
%
% wf = loaddefaultwaveforms creates a waveform struct with fields:
%      W - (M*T*K) upsampled waveforms
%      d - EC data struct from which they were extracted (d.samplefreq has
%          sampling rate)
%      fac - upsampling factor
%      freqs - firing counts from the datasets they were extracted

wf = load('data/waveforms_2.5ms_e_fac3_tv_dbk.mat');

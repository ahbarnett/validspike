function image_spike_shapes(X)
% IMAGE_SPIKE_SHAPES - show image of single events matrix X
%
% image_spike_shapes(X) stacks the M waveforms down the y axis, events are cols
%
% Inputs:
%  X - M*Nt*Ns array of waveform data or events
%
% Barnett 12/22/14

[M T Ns] = size(X);
figure; imagesc(reshape(permute(X,[2 1 3]),[M*T Ns])); colorbar;

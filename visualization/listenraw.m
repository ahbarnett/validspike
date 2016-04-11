function p = listenraw(d,ch)
% LISTENRAW - dump two channels of raw data to WAV file and play it
%
% p = listenraw(d,ch) where d is an extracellular data struct and ch = [m1 m2]
% is a vector of two channel numbers (1-indexed) plays in real time as
% stereo audio. Also writes to raw.wav (overwriting any previous such file).
% Outputs:
% p - an audioplayer object that can, eg, be stopped via stop(p)
%
% Example: see test_listenraw
%
% Barnett 11/14/14

if nargin<1, test_listenraw; return; end

y = d.A(ch,:);         % extract 2 chans. todo: sanity check
y = y/max(abs(y(:)));  % prevent clipping
p = audioplayer(y', d.samplefreq); % note time is downwards
play(p);

% todo: make this optional
%wavwrite(y',d.samplefreq,'/tmp/raw.wav');  % deprecated and replaced by (R2014b)...
audiowrite('raw.wav',y',d.samplefreq);  % post-R2014b...

%%%%
function test_listenraw
d = loaddemodata; p = listenraw(d,[1 2]); pause(2); stop(p);

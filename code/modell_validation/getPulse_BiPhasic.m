function [t, sig] = getPulse_BiPhasic(t, sig, pulseWidthUs, pauseWidthUs, currentMa, delayStartUs_, delayEndUs_, debugPlot_)
%GETPULSE_BIPHASIC returns the time and current vector for a biphasic pulse

% defaults for optional parameters
if ~exist('delayStartUs_', 'var'), delayStartUs_ = 0; end
if ~exist('delayEndUs_', 'var'),   delayEndUs_ = 0;  end
if ~exist('debugPlot_', 'var'),    debugPlot_ = false;  end

% sample time
if (length(t) < 2)
    ts = 8-3;   % 1Mhz / 1us
else
    ts = t(2) -t(1);
end

%% build signal

% start delay
if (delayStartUs_ > 0)
    sig = addSteadyState(sig, ts, delayStartUs_/1000/1000, 0);
end

% first part of the pulse
sig = addSteadyState(sig, ts, pulseWidthUs/1000/1000, currentMa/1000);
% 100 us pause
if (pauseWidthUs > 0)
    sig = addSteadyState(sig, ts, pauseWidthUs/1000/1000, 0);
end
% second part of the pulse
sig = addSteadyState(sig, ts, pulseWidthUs/1000/1000, currentMa/1000*-1);

% end delay
if (delayEndUs_ > 0)
    sig = addSteadyState(sig, ts, delayEndUs_/1000/1000, 0);
end


%% update time
n = length(sig) - length(t);
if isempty(t)
    t = 0 : ts : ts*(n-1);
else
    t(end+1 : end+n) = t(end)+ts : ts : t(end)+ts*n;
end


%% control plot
if debugPlot_
    figure(); hold on;
    plot(t*1000*1000, sig);
    legend('signal / current');
    xlabel('time [us]');
    ylabel('current [A]');
    grid on;
end


%% Helper Functions
function nSamples = nSamples4Time(sampleRate, duration)
nSamples = round(duration / sampleRate);
end

function [vec, iStart, iStop] = addLinearRise(vec, ts, duration, vStart, vEnd)
% linear rise
nSamples = nSamples4Time(ts, duration);
iStart = length(vec) +1;
iStop  = iStart -1 +nSamples;
vec(iStart:iStop) = ((1:nSamples)/nSamples *(vEnd -vStart)) +vStart;  % linear rise + normalisation to 1 + normalisation 2 final value
end

function [vec, iStart, iStop] = addSteadyState(vec, ts, duration, vEnd)
% steady for x sec
nSamples = nSamples4Time(ts, duration);
iStart = length(vec) +1;
iStop  = iStart -1 +nSamples;
vec(iStart:iStop) = vEnd;  % final value
end

end
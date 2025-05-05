function out = findPulses(in)

%%%
%
%   Copyright MIT (c) 2017 Arthur Morrissette
%
%%%

inDiff = diff(in);
[pulseStarts] = find(inDiff == 1);
[pulseEnds] = find(inDiff == -1);

if isempty(pulseEnds) || isempty(pulseStarts)
    pulseEnds = [];
    pulseStarts = [];
end

if length(pulseStarts) > length(pulseEnds)
    pulseStarts = pulseStarts(1:end-1);
elseif length(pulseStarts) < length(pulseEnds)
    if pulseStarts(1) < pulseEnds(1)
        pulseEnds = pulseEnds(1:end-1);
    else
        pulseEnds = pulseEnds(2:end);
    end
end

% Add removal of artifact pulses (equal to or shorter than 2 samples or longer than 20 samples)
pulseWidths = pulseEnds-pulseStarts;
artifacts = pulseWidths<=2;
pulseStarts(artifacts) = [];
pulseEnds(artifacts) = [];
% artifacts2 = pulseWidths>20;
% pulseStarts(artifacts2) = [];
% pulseEnds(artifacts2) = [];

% Save as output
out.starts = pulseStarts;
out.ends = pulseEnds;

end
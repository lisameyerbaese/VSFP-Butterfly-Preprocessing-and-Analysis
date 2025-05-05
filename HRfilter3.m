function HRwin = HRfilter3(imgD, fDate, fNum, Fs, mouseID, auto)

ROIs = 0;

if ROIs == 1
    for x = 1:4
            region = roiSelect(fDate, fNum, mouseID);
            xLoc(x) = region(1);
            yLoc(x) = region(3);
    end
else
    xLoc = [25,25,75,75];
    yLoc = [25,25,75,75];
end

sZ = size(imgD, 3);
Df = ones(length(xLoc),sZ);

ind = 1;
temp = imgD(1,1,ind);
while temp == 0
    ind = ind +1;
    temp = imgD(1,1,ind);
end

DfTest = (squeeze(imgD(1,1,:))-imgD(1,1,ind))/imgD(1,1,ind);
[pxxTest,~] = pwelch(DfTest,[],[],[],Fs);
pxx = ones(length(xLoc), length(pxxTest));
F = ones(1, length(pxxTest));

for x = 1:length(xLoc)
        xL = xLoc(x);
        yL = yLoc(x);
        Df(x, :) =(squeeze(imgD(xL,yL,:))-imgD(xL,yL,1))/imgD(xL,yL,1);
        [pxx(x,:),F] = pwelch(Df(x,:),[],[],[],Fs);
end

% Crop freq range to analyze from 8-15 Hz
f_inds = find(F<15 & F>8);
F_crop = F(f_inds);
mean_pxx = mean(pxx);
smooth_pxx = smooth(mean_pxx);
pxx_crop = smooth_pxx(f_inds);

%Remove linear trend in data
[~, locs, ~, p] = findpeaks(pxx_crop);

peakLoc = locs((p == max(p)));
HRfreq = F_crop(peakLoc);

HRwin = [HRfreq - 1.5, HRfreq + 1.5];
end
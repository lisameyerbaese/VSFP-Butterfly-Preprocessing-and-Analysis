% this function saves figures for you as both 
% a .fig and a .png file 

function saveFig(figs, fname, FolderName)

%     temp = [fname, '.fig'];
%     savefig(figs, fullfile(FolderName, temp));
    temp = [fname, '.png'];
    saveas(figs, fullfile(FolderName, temp));

end
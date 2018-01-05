clear;
close all;

w = warning ('off','all');

wsizes = {50, 100, 200, 400, 600, 800, 1000, 1250, 1275, 1300};
sigmas = {[0.1], [0.2], [0.225], [0.25], [0.275], [0.5], [1]};

nsig = length(sigmas);

filter_size = [1,3,6,10];
overlap_ratio = [50]; % Adaptative
ms = [3]; % Inlier range
nPeaks = [3];

ratesac = {};
ratehst = {};
ratepks = {};
ratekms = {};
ratecor = {};

afile = fopen (['SIM_LAUNCH nPeaks: ', num2str(nPeaks), '; OR: ', num2str(overlap_ratio),'.txt'], 'w');

datasetdir = '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/'
% datasetdir = '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/'

for fi=1:length(filter_size)
    for or=1:length(overlap_ratio)
        for np=1:length(nPeaks);
            for mi=1:length(ms)
                [ratesac{fi,np,or,mi}, ratehst{fi,np,or,mi}, ratepks{fi,np,or,mi}, ratekms{fi,np,or,mi}, ratecor{fi,np,or,mi}, clevel_pks{fi,np,or,mi}, clevel_hst{fi,np,or,mi}, clevel_sac{fi,np,or,mi}, clevel_km{fi,np,or,mi}] = SOURCE_MULTICAM_DATASET_F (datasetdir, wsizes, sigmas, filter_size(fi), overlap_ratio(or), ms(mi), nPeaks(np));

    %             figure;
    %             set(gcf,'name',['inlrRange = ', num2str( ms(mi)), ' | OR = ', num2str(overlap_ratio(or)), 'nPeaks = ', num2str(nPeaks(np))]);
    %     
    %             subplot(2,2,1);
    %             mesh([wsizes{:}],[sigmas{:}],ratehst{fi,np,or,mi}(:,:));
    %             title (['Probability rate CCC HST']);
    %             xlabel ('window size');
    %             ylabel ('sigmas');
    %     
    %             subplot(2,2,2);
    %             mesh([wsizes{:}],[sigmas{:}],ratepks{fi,np,or,mi}(:,:));
    %             title (['Probability rate PKS']);
    %             xlabel ('window size');
    %             ylabel ('sigmas');
    %     
    %             subplot(2,2,3);
    %             mesh([wsizes{:}],[sigmas{:}],ratesac{fi,np,or,mi}(:,:));
    %             title (['Probability rate SAC']);
    %             xlabel ('window size');
    %             ylabel ('sigmas');
    %     
    %             subplot(2,2,4);
    %             mesh([wsizes{:}],[sigmas{:}],ratekms{fi,np,or,mi}(:,:));
    %             title (['Probability rate KMS']);
    %             xlabel ('window size');
    %             ylabel ('sigmas');



                fprintf (afile, '\n\nResults. Filter: %d, Overlap: %d, inlrRange: %d, nPeaks: %d \n', filter_size(fi), overlap_ratio(or), ms(mi), nPeaks(np));

                ccc_max = max(ratesac{fi,np,or,mi}(:));
                ccc_avg = mean(ratesac{fi,np,or,mi}(:));
                ccc_min = min(ratesac{fi,np,or,mi}(:));
                fprintf (afile, 'SAC CCC min|avg|max: %.2f | %.2f | %.2f \n', ccc_min, ccc_avg, ccc_max);
                fprintf (afile, 'SAC CLV min|avg|max: %.2f | %.2f | %.2f \n', min(clevel_sac{fi,np,or,mi}(:)), mean(clevel_sac{fi,np,or,mi}(:)), max(clevel_sac{fi,np,or,mi}(:)));
                fprintf (afile, '300: %.2f | 800: %.2f | 1250: %.2f \n', max(ratesac{fi,np,or,mi}(:,3)), max(ratesac{fi,np,or,mi}(:,6)), max(ratesac{fi,np,or,mi}(:,8)));

                ccc_max = max(ratekms{fi,np,or,mi}(:));
                ccc_avg = mean(ratekms{fi,np,or,mi}(:));
                ccc_min = min(ratekms{fi,np,or,mi}(:));
                fprintf (afile, 'KMS CCC min|avg|max: %.2f | %.2f | %.2f \n', ccc_min, ccc_avg, ccc_max);
                fprintf (afile, 'KMS CLV min|avg|max: %.2f | %.2f | %.2f \n', min(clevel_km{fi,np,or,mi}(:)), mean(clevel_km{fi,np,or,mi}(:)), max(clevel_km{fi,np,or,mi}(:)));
                fprintf (afile, '300: %.2f | 800: %.2f | 1250: %.2f \n', max(ratekms{fi,np,or,mi}(:,3)), max(ratekms{fi,np,or,mi}(:,6)), max(ratekms{fi,np,or,mi}(:,8)));

                ccc_max = max(ratehst{fi,np,or,mi}(:));
                ccc_avg = mean(ratehst{fi,np,or,mi}(:));
                ccc_min = min(ratehst{fi,np,or,mi}(:));
                fprintf (afile, 'HST CCC min|avg|max: %.2f | %.2f | %.2f \n', ccc_min, ccc_avg, ccc_max);
                fprintf (afile, 'HST CLV min|avg|max: %.2f | %.2f | %.2f \n', min(clevel_hst{fi,np,or,mi}(:)), mean(clevel_hst{fi,np,or,mi}(:)), max(clevel_hst{fi,np,or,mi}(:)));
                fprintf (afile, '300: %.2f | 800: %.2f | 1250: %.2f \n', max(ratehst{fi,np,or,mi}(:,3)), max(ratehst{fi,np,or,mi}(:,6)), max(ratehst{fi,np,or,mi}(:,8)));

                ccc_max = max(ratepks{fi,np,or,mi}(:));
                ccc_avg = mean(ratepks{fi,np,or,mi}(:));
                ccc_min = min(ratepks{fi,np,or,mi}(:));
                fprintf (afile, 'PKS CCC min|avg|max: %.2f | %.2f | %.2f \n', ccc_min, ccc_avg, ccc_max);
                fprintf (afile, 'PKS CLV min|avg|max: %.2f | %.2f | %.2f \n', min(clevel_pks{fi,np,or,mi}(:)), mean(clevel_pks{fi,np,or,mi}(:)), max(clevel_pks{fi,np,or,mi}(:)));
                fprintf (afile, '300: %.2f | 800: %.2f | 1250: %.2f \n', max(ratepks{fi,np,or,mi}(:,3)), max(ratepks{fi,np,or,mi}(:,6)), max(ratepks{fi,np,or,mi}(:,8)));

                ccc_max = max(ratecor{fi,np,or,mi}(:));
                ccc_avg = mean(ratecor{fi,np,or,mi}(:));
                ccc_min = min(ratecor{fi,np,or,mi}(:));
                fprintf (afile, 'COR min|avg|max: %.2f | %.2f | %.2f \n', ccc_min, ccc_avg, ccc_max);

            end
        end
    end
end















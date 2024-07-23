function Ind_outlier = DetectOutlier_TRUST(dif_all, control_all, option)
% Based on Structural Correlation-Based Outlier Rejection (SCORE)
% Algorithm  for Arterial Spin Labeling Time Series, JMRI, 2017;45:1786–1797
if nargin < 3
    option = 'whole';
end
brain_thresh_scale = 0.5;
[Nrow, Ncol, Ndyn] = size(dif_all);
dif_all = reshape(dif_all, [Nrow, Ncol, Ndyn]);
control_all = reshape(control_all, [Nrow, Ncol, Ndyn]);
mean_diff_ratio = zeros(Ndyn, 1);
for iDyn = 1:Ndyn
    cnt_control = control_all(:,:,iDyn);
    cnt_dif = dif_all(:,:,iDyn);
    tmpmat = cnt_control(23:42,23:42);
    tmpvox = tmpmat(tmpmat>0);
    brain_thresh = mean(tmpvox)*brain_thresh_scale;
    brain_mask = cnt_control > brain_thresh;
    if strcmp(option, 'half') || strcmp(option, 'quarter')
        [rowInd, colInd] = find(brain_mask);
        centroid_row = round(mean(rowInd));
        brain_mask(1:centroid_row, :) = 0; % only focusing on the lower half of the brain
        if strcmp(option, 'quarter')
            [rowInd, colInd] = find(brain_mask);
            centroid_row = round(mean(rowInd));
            brain_mask(1:centroid_row, :) = 0; % only focusing on the lower quarter of the brain
        end
    end
    mean_diff_ratio(iDyn) = mean(cnt_dif(brain_mask))./mean(cnt_control(brain_mask));
end
median_diff_ratio = median(mean_diff_ratio);
med_abs_dev = median(abs(mean_diff_ratio - median_diff_ratio));
robsd = 1.4826*med_abs_dev; % robust standard deviation
thre_robsd = 2.5;
figure,plot(mean_diff_ratio);
AvgLine = refline([0, median_diff_ratio]);AvgLine.Color = 'r';
StdLineUp = refline([0, median_diff_ratio + thre_robsd*robsd]);StdLineUp.LineStyle = '--';
StdLineDown = refline([0, median_diff_ratio - thre_robsd*robsd]);StdLineDown.LineStyle = '--';

Ind_outlier = find(mean_diff_ratio > median_diff_ratio + thre_robsd*robsd | mean_diff_ratio < median_diff_ratio - thre_robsd*robsd);
hold on;plot(Ind_outlier, mean_diff_ratio(Ind_outlier), 'r*');hold off;
set(gca, 'fontsize',15);
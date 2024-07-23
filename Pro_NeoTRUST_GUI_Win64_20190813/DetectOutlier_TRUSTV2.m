function Ind_outlier = DetectOutlier_TRUSTV2(dif_all, control_all, option)
    % Function to do outlier images detection based on L1 ratio norm method

    rep_num = 3; eTE_num = 4;
    brain_thresh_scale = 0.5;
    L1_factor = [1 1.0863 1.1670 1.2329];
    columnnum = size(dif_all,1);
    rownum = size(dif_all,2);
    L1 = zeros(rep_num*eTE_num,1);

    for idyn = 1:rep_num*eTE_num
        cnt_control = control_all(:,:,idyn);
        cnt_diff = abs(dif_all(:,:,idyn));
        tmpmat = cnt_control(round(columnnum/3):round(2*columnnum/3),round(rownum/3):round(2*rownum/3));
        tmpvox = tmpmat(tmpmat>0);
        brain_thresh = mean(tmpvox)*brain_thresh_scale;
        brain_mask = cnt_control > brain_thresh;

        if strcmp(option, 'whole')
            L1(idyn) = mean(cnt_diff(brain_mask));
        elseif strcmp(option, 'bottom')
            % get the bottom part of brain to calculate L1
            brain_mask_bottom = brain_mask;
            [rowInd, ~] = find(brain_mask);
            centroid_row = round(mean(rowInd));
            brain_mask_bottom(1:centroid_row, :) = 0;
            L1(idyn) = mean(cnt_diff(brain_mask_bottom));
        end
    end

    for ieTE =  1:eTE_num
        L1_normed((ieTE-1)*rep_num+1 : ieTE*rep_num,:) = L1((ieTE-1)*rep_num+1 : ieTE*rep_num,:)*L1_factor(ieTE);
    end
    
    if strcmp(option, 'whole')
        fprintf('Using whole brain to calculate L1\n');
    elseif strcmp(option, 'bottom')
        fprintf('Using bottom brain to calculate L1\n');
    end

    Ind_outlier = find(L1_normed > 26);
end
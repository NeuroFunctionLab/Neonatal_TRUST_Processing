function Ind = GetRefInd(varargin)
    % Function to calculate best reference image using different methods
    % based on maxmizing total similarity between reference images

    I = varargin{1};

    if strcmp(varargin{2}, 'whole')
        bApplyMask = 0;
    elseif strcmp(varargin{2}, 'mask')
        bApplyMask = 1;
        MaskThresholdFactor = 0.4;
        fprintf('Applying eTE0 mask to all eTE to find optimal reference images\n');
    else
        error('Input incorrect');
    end
    
    if strcmp(varargin{3}, 'MI')    
        fprintf('Using MI to caluclate similarity\n');
    elseif strcmp(varargin{3}, 'Normed MI')
        fprintf('Using normalized MI to caluclate similarity\n');
    elseif strcmp(varargin{3}, 'normxcorr2')
        fprintf('Using norm cross correlation to caluclate similarity\n');
    end

    TotalSimi = zeros(size(I,3),1);
    Comb = zeros(size(I,3),4);
    i = 1;
    for ieTE0 = 1:3
        for ieTE40 = 1:3
            for ieTE80 = 1:3
                for ieTE160 = 1:3
                    Ind1 = GetInd(1,2,ieTE0); I1 = I(:,:,Ind1);
                    Ind2 = GetInd(2,2,ieTE40); I2 = I(:,:,Ind2);
                    Ind3 = GetInd(3,2,ieTE80); I3 = I(:,:,Ind3);
                    Ind4 = GetInd(4,2,ieTE160); I4 = I(:,:,Ind4);
                    
                    if bApplyMask
                        tmpmat = I1(23:42,23:42);
                        tmpvox = tmpmat(tmpmat>0);
                        brain_thresh = mean(tmpvox) * MaskThresholdFactor;
                        MaskUsed = I1 > brain_thresh;
                        MaskUsed = imfill(MaskUsed, 'holes');
                        MaskUsed = logical(MaskUsed);
                        if ~strcmp(varargin{3}, 'normxcorr2')
                            I1 = I1.*MaskUsed; I2 = I2.*MaskUsed; I3 = I3.*MaskUsed; I4 = I4.*MaskUsed;
                        end
                    end

                    if strcmp(varargin{3}, 'MI')
                        MI12 = mi(I1,I2);
                        MI13 = mi(I1,I3);
                        MI14 = mi(I1,I4);
                        TotalSimi(i) = MI12 + MI13 + MI14;
                    elseif strcmp(varargin{3}, 'Normed MI')
                        NormMI12 = MI2(I1,I2,'Normalized');
                        NormMI13 = MI2(I1,I3,'Normalized');
                        NormMI14 = MI2(I1,I4,'Normalized');
                        TotalSimi(i) = NormMI12 + NormMI13 + NormMI14;
                    elseif strcmp(varargin{3}, 'normxcorr2')
                        if bApplyMask
                            corr12 = normxcorr2(I1(MaskUsed),I2(MaskUsed));
                            corr13 = normxcorr2(I1(MaskUsed),I3(MaskUsed));
                            corr14 = normxcorr2(I1(MaskUsed),I4(MaskUsed));
                            len = length(I1(MaskUsed));
                            TotalSimi(i) = corr12(len) + corr13(len) + corr14(len);
                        else
                            corr12 = normxcorr2(I1,I2);
                            corr13 = normxcorr2(I1,I3);
                            corr14 = normxcorr2(I1,I4);
                            TotalSimi(i) = corr12(64,64) + corr13(64,64) + corr14(64,64);
                        end
                    end

                    Comb(i,1) = Ind1; Comb(i,2) = Ind2; Comb(i,3) = Ind3; Comb(i,4) = Ind4;
                    i = i+1;
                end
            end
        end
    end

    MaxInd = find(TotalSimi == max(TotalSimi));
    Ind = Comb(MaxInd,:);
end
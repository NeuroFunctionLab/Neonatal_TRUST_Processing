function [t2,ci,Yv, LabelEff]=trustcode_pl_GUI(filename1,matrix,seq,hct)
% TRUST acquired with any eTE order
% Input:
% filename: .REC file directory
% nlist: pixel number of interest
% roiflag: dynamic drawing roi for different TE set '1' or static drawing roi from ...
%  first TE set '0'
% datastructure: data structure N,N,repeat*2 (include label and control, # of total dynamics each TE weighted scan has)
% nrep: n out of 4 repeat measurements want to use

% Output:
% T2 value in milisecond
% CI 95%
% Yv
% LabelEff: labelling efficiency
% Author: Peiying Liu, 08/23/2011
% Modified by Dengrong Jiang, 5/14/2018
bSelectRef = 1; % use normalized cross correlation or mutual information to select the best set of reference images
LCflag = 1; % use label or control image as target image for realignment Label: LCflag=1, Control: LCflag =2
close all;
warning off;
[fname_path, fname_body, fname_ext] = fileparts(filename1);
filename = [fname_path filesep fname_body];
etetemp=unique(seq);
etelen=length(etetemp);
te=sort(etetemp,'ascend')';
array_temp = cell(etelen, 1);
len_temp = zeros(etelen, 1);
for iTE = 1:etelen
    array_temp{iTE} = find(seq==te(iTE));
    len_temp(iTE) = length(array_temp{iTE});
end
roiflag=0;
bloodt1=1624;
ti=1020;
rowspacing=2.5;
columnspacing=2.5;
slicespacing=5;
nlist=4; % number of voxels for averaging

columnnum=matrix(1);
rownum=matrix(2);
totaldyn=matrix(3);
if strcmpi(fname_ext, '.REC') % read .REC for Philips data
    noalign=read_images_vms(strcat(filename,'.REC'),rownum,columnnum, 'int16',totaldyn);
% Dengrong Jiang, 5/7/2018: Add codes to read DICOM images and reorder them
elseif strcmpi(fname_ext, '.IMA') || strcmpi(fname_ext, '.DCM') % read .IMA or .DCM for Siemens data
    if strcmpi(fname_ext, '.IMA') % read .IMA for Siemens data
        PathList = dir([fname_path,filesep,'*.IMA']);
    elseif strcmpi(fname_ext, '.DCM')
        PathList = dir([fname_path,filesep,'*.DCM']);
    end    
    Filename_list = extractfield(PathList,'name');
%     Ind_1st_file = find(strcmp(Filename_list, strcat(fname_body, fname_ext)));
    Ind_1st_file = find(contains(Filename_list, fname_ext), 1);
    noalign = zeros(matrix);
    acq_ind = zeros(totaldyn, 1);
    % read the 1st image
    [cntIm] = dicomread(strcat(fname_path, filesep, Filename_list{Ind_1st_file}));
    cntInfo = dicominfo(strcat(fname_path, filesep, Filename_list{Ind_1st_file}));
    if isfield(cntInfo, 'NumberOfFrames') % multi-frame enhanced DICOM
        noalign = squeeze(cntIm);
    else % single-frame DICOM, read each dynamic one by one
        noalign(:,:,1) = cntIm;
        acq_ind(1) = cntInfo.AcquisitionNumber;
        for iIm = 2:totaldyn
            [cntIm] = dicomread(strcat(fname_path, filesep, Filename_list{Ind_1st_file + iIm - 1}));
            cntInfo = dicominfo(strcat(fname_path, filesep, Filename_list{Ind_1st_file + iIm - 1}));
            noalign(:,:,iIm) = cntIm;
            acq_ind(iIm) = cntInfo.AcquisitionNumber;
        end
        % reorder the images
        [~, acq_ind_sort] = sort(acq_ind, 'ascend');
        noalign = noalign(:,:,acq_ind_sort);
    end
end
    
clear con_all
clear lab_all

%--------------- Comment out the SPM realignment code in adult TRUST pipeline ---------------
% for i=1:etelen
%     delete_if_exist(strcat(filename,'_',int2str(te(i))));
% end
% 
% % regroup different eTEs ------------------
% for i=1:etelen
%     tempstr=strcat(filename,'_',int2str(te(i)));
%     mkdir(tempstr);
%     dynnum=len_temp(i);
%     dyntemp=array_temp{i};
%     for j=1:dynnum
%         if j<10
%             a=strcat('00',int2str(j));
%         elseif j<100
%             a=strcat('0',int2str(j));
%         else
%             a=int2str(j);
%         end
%         tempstr1=strcat(tempstr,filesep,'img_', a);
%         ndyn=dyntemp(j);
%         write_ANALYZE(noalign(:,:,ndyn),tempstr1,[rownum columnnum 1], [rowspacing columnspacing slicespacing], 1,16);
%     end
% end
% 
% % realign within different eTEs to its fist control image--------------------
% for i=1:etelen
%     clear img_te P;
%     tempstr=strcat(filename,'_',int2str(te(i)));
%     target=strcat(strcat(tempstr,filesep,'img_002.img'));
%     P{1} = spm_select('list',tempstr,'^img.*\.img'); %CR original datafiles
% %     P{1} = spm_get('Files',tempstr,'img*.img'); %CR original datafiles
%     [dynnum,nname]=size(P{1});
%     dyn(i)=dynnum;
%     for count=1:dynnum
%         if count<10
%             a=strcat('00',int2str(count));
%         elseif count<100
%             a=strcat('0',int2str(count));
%         else
%             a=int2str(count);
%         end
%         source=strcat(strcat(tempstr,filesep,'img_', a, '.img'));
%         delete_if_exist(['img_' a '.mat']);
%         other=source;
%         spm_defaults;
%         defs=defaults.realign;
%         FlagsC = struct('quality',defs.estimate.quality,'fwhm',5,'rtm',0);
%         spm_realign([target;source],FlagsC);
%         FlagsR = struct('interp',defs.write.interp,...
%             'wrap',defs.write.wrap,...
%             'mask',defs.write.mask,...
%             'which',2,'mean',1);
%         FlagsR.which = 2; FlagsR.mean = 0;
%         spm_reslice([target;source],FlagsR);
%             
%         thisimg=loadimage(strcat(tempstr,filesep,'rimg_', a, '.img'), 16);
%         img_te(:,:,count)=thisimg(:,:);
%         img_all{i}=img_te;
%     end
% end

%--------------- Copy the SPM realignment code in Peiying's TRUST pipeline ---------------
spm_defaults;
global defaults;
trustfilename='img_';

for i=1:etelen
    delete_if_exist(strcat(fname_path, filesep,int2str(te(i))));
end
% regroup different eTEs ------------------
for i=1:etelen
    tempstr=strcat(strcat(fname_path, filesep,int2str(te(i))));
    mkdir(tempstr);
    dynnum=len_temp(i);
    dyntemp=array_temp{i};
    for j=1:dynnum
        if j<10
            a=strcat('00',int2str(j));
        elseif j<100
            a=strcat('0',int2str(j));
        else
            a=int2str(j);
        end
        tempstr1=strcat(tempstr,filesep,trustfilename, a);
        ndyn=dyntemp(j);
        write_ANALYZE(noalign(:,:,ndyn),tempstr1,[rownum columnnum 1], [rowspacing columnspacing slicespacing], 1,16);
    end
end

% ------ realign within each eTE ------
if bSelectRef
    RefInd = GetRefInd(noalign,'mask','normxcorr2', LCflag);
    RepTime4eTE = ceil(RefInd/etelen/2);
    fprintf('Indices of reference images: %d\n', RepTime4eTE);
end

dyn = zeros(etelen, 1);
motion_vec = cell(etelen, 1);
for i=1:etelen
    motion_vec{i} = zeros(dynnum, 6);
    tempstr=strcat(fname_path, filesep,int2str(te(i)),filesep,trustfilename);
    if bSelectRef
        target_dyn = (RepTime4eTE(i)-1)*2 + LCflag;
    else
        target_dyn = LCflag;
    end
    if target_dyn < 10
        target_fn = strcat(trustfilename, '00',int2str(target_dyn));
    elseif target_dyn < 100
        target_fn = strcat(trustfilename, '0',int2str(target_dyn));
    else
        target_fn = [trustfilename, int2str(target_dyn)];
    end
    targetfile = spm_get('Files', strcat(fname_path, filesep, int2str(te(i))), [target_fn, '.img']);

    FlagsC = struct('quality',defaults.realign.estimate.quality,'fwhm',5,'rtm',0);
    dyn(i)=dynnum;
    for j=1:dynnum
        if j<10
            a=strcat('00',int2str(j));
        elseif j<100
            a=strcat('0',int2str(j));
        else
            a=int2str(j);
        end
        delete_if_exist([tempstr,a, '.mat']);
        source=strcat(strcat(tempstr, a, '.img'));
        spm_realign([targetfile;source],FlagsC);
        % read motion vector
        cnt_motion_vec = importdata([fname_path, filesep, int2str(te(i)), filesep, 'rp_', target_fn, '.txt']);
        motion_vec{i}(j, :) = cnt_motion_vec(2, :);
        which_writerealign = 2;
        mean_writerealign = 0;
        FlagsR = struct('interp',defaults.realign.write.interp,...
            'wrap',defaults.realign.write.wrap,...
            'mask',defaults.realign.write.mask,...
            'which',which_writerealign,'mean',mean_writerealign); 
        spm_reslice([targetfile;source],FlagsR);
    end
end

imgall=zeros(rownum,columnnum,etelen*dynnum);
ncon=0;
for i=1:etelen
    for count=1:dynnum
        if count<10
            a=strcat('00',int2str(count));
        elseif count<100
            a=strcat('0',int2str(count));
        else
        a=int2str(count);
        end
        P = spm_get('Files',strcat(fname_path, filesep,int2str(te(i))),['r',trustfilename '*.img']); %CR original datafiles
        thisimg=loadimage(P(count,:), 16);
        ncon=ncon+1;
        imgall(:,:,ncon)=thisimg(:,:);
    end
end
write_ANALYZE(imgall,[fname_path filesep 'img_all.img'],[rownum, columnnum, totaldyn],[rowspacing columnspacing slicespacing],1,16);
save([fname_path filesep 'motion_vec.mat'], 'motion_vec', 'RefInd', 'RepTime4eTE');

img_all = cell(etelen, 1);
for i=1:etelen
    img_all{i} = imgall(:,:,(i-1)*dynnum+1:i*dynnum);
end
%--------------- End of SPM realignment code ---------------

con_mean=zeros(rownum,columnnum,etelen);
lab_mean=zeros(rownum,columnnum,etelen);
dif_mean=zeros(rownum,columnnum,etelen);
for i=1:etelen
    clear tempimg
    tempimg=img_all{i};
    con_all{i}=tempimg(:,:,2:2:end);
    lab_all{i}=tempimg(:,:,1:2:end);
    dif_all{i}=con_all{i}-lab_all{i};
    
    con_mean(:,:,i)=mean(con_all{i},3);
    lab_mean(:,:,i)=mean(lab_all{i},3);
    dif_mean(:,:,i)=mean(dif_all{i},3);
end

mask=zeros(rownum,columnnum,etelen);

% Detect outliers --------------------------------------------------------
control_all_mat = imgall(:,:,2:2:end);
diff_all_mat = imgall(:,:,2:2:end) - imgall(:,:,1:2:end);
close all;
show_imgs_sc(control_all_mat,1:totaldyn/2, 0,1000);
saveas(gcf, [fname_path, filesep, 'control_all.fig']);
show_imgs_sc(abs(diff_all_mat),1:totaldyn/2, 0,100);
saveas(gcf, [fname_path, filesep, 'diff_all.fig']);
% ind2exclude = DetectOutlier_TRUSTV2(diff_all_mat, control_all_mat, 'bottom');
% ind2exclude
ind2exclude = input('Manual exclusion: ');
% check if any excluded dynamics are eTE0
ind2include_eTE0 = [1:dynnum/2];
ind2exclude_eTE0 = ind2exclude(ind2exclude <= dynnum/2);
ind2include_eTE0(ind2exclude_eTE0) = [];

% draw ROI --------------------------------------------------------
for t=1:etelen
    % draw ROI
    anat=mean(dif_all{1}(:,:,ind2include_eTE0), 3);
    
    droi=zeros(rownum, columnnum);
    figure;
    himage=imshow(anat,[min(anat(:)) max(anat(:))]);
    temp=roipoly();colorbar;
    droi=droi | temp;
    close;
    clear mylist;
    mycount=0;        
    for i=1:rownum
        for j=1:columnnum
            if droi(i,j)==1
                mycount=mycount+1;
                mylist(mycount,1)=i;
                mylist(mycount,2)=j;
                mylist(mycount,3)=anat(i,j,1);
            end
        end
    end
    [Y,I]=sort(mylist(:,3),'descend');
    if roiflag ==0
        mask=zeros(rownum,columnnum);
        for tt=1:nlist
            mask(mylist(I(tt),1),mylist(I(tt),2))=1;
        end
        break;
    else
        for tt=1:nlist
            mask(mylist(I(tt),1),mylist(I(tt),2),t)=1;
        end
    end
end
display_fmri(anat,mask);

% fitting for T2 ----------------------------------------
nte=etelen;
for i=1:nte
    clear temp3 temp5;
    temp3=dif_all{i};
    temp5=con_all{i};
    if roiflag ==0
        temp4=average_actpixel1(1,mask,temp3);
        b_all{i}=temp4(:);
        temp6=average_actpixel1(1,mask,temp5);
        c_all{i}=temp6(:);
    else
        temp4=average_actpixel1(1,mask(:,:,i),temp3);
        b_all{i}=temp4(:);
        temp6=average_actpixel1(1,mask(:,:,i),temp5);
        c_all{i}=temp6(:);
    end   
end

%%% T2 fitting by all b points
tetemp=seq(1:2:end);
te_all=sort(tetemp,'ascend');
b_all_vector = [];
c_all_vector = [];
for i=1:nte
    b_all_vector = [b_all_vector; b_all{i}];
    c_all_vector = [c_all_vector; c_all{i}];
end

% Dengrong Jiang, 5/14/2018: Compute labelling efficiency
Ind4LabelEff = find(te_all == te(1));
LabelEff = mean(b_all_vector(Ind4LabelEff))/mean(c_all_vector(Ind4LabelEff));

dyn0=0;
for i=1:etelen-1
    b_all_vector(dyn0+1:dyn0+dyn(i)/2)=exp(-(ti-te(end))/bloodt1)/exp(-(ti-te(i))/bloodt1)*b_all_vector(dyn0+1:dyn0+dyn(i)/2);
    c_all_vector(dyn0+1:dyn0+dyn(i)/2)=exp(-(ti-te(end))/bloodt1)/exp(-(ti-te(i))/bloodt1)*c_all_vector(dyn0+1:dyn0+dyn(i)/2);
    dyn0=dyn0+dyn(i)/2;
end

% Exclude dynamics with strong motion artifacts
ind2fit = [1:length(te_all)];
ind2fit(ind2exclude) = [];

[temp1,resid, jacob]=nlinfit_hlu(te_all(ind2fit)/1000, b_all_vector(ind2fit),'monexp_model',[100,10]);
t2=1000./temp1(2);

conintval=nlparci(temp1,resid, jacob); %95% confidence interval for estimates
aa=1000./conintval;
ci=aa(2,:);

[temp1_c,resid_c, jacob_c]=nlinfit_hlu(te_all(ind2fit)/1000, c_all_vector(ind2fit),'monexp_model',[100,10]);
conintval_c=nlparci(temp1_c,resid_c, jacob_c); %95% confidence interval for estimates
t2_c=1000./temp1_c(2);
aa_c=1000./conintval_c;
ci_c=aa_c(2,:);
Yv=100*neoT2toY(t2,hct);

figure,plot(transpose(reshape(te_all, [3,4])), transpose(reshape(b_all_vector,[3,4])),'*');
legend({'rep1','rep2','rep3'});
set(gca, 'fontsize',15);
saveas(gcf, [fname_path, filesep, 'ROI_signal_all.fig']);

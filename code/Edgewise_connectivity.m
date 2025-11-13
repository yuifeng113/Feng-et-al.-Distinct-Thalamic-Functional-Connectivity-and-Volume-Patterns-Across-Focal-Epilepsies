%% Step1 mapping source data
clear
clc
%addpath 'C:\Users\skgtxfe\OneDrive - University College London\Downloads\fdr_bh'

% Load subjects
%infotable = 'C:\Users\skgtxfe\OneDrive - University College London\Data_Code\ICHdatabase\Focal_Con_270125.xlsx';
infotable = '/Users/fxy/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Data_Code/ICHdatabase/Focal_Con_180825.xlsx';
%infotable = 'C:\Users\skgtxfe\OneDrive - University College London\Data_Code\ICHdatabase\TLE_Con_working_090824.xlsx';
datapre = readtable(infotable);% close the excel if unreadable!
%data=datapre(datapre.Include_fmri==1,:);
data=datapre(datapre.Include_vol==1,:);% for tiv adjustment
subjects = data.conn_id;%subjects = importdata('\sbj.txt'); % List of your subjects

% Controls (all 1.5T)
idx{1} = find(strcmp(data.Focal,'N')); name{1} = 'Con';
%idx{1} = find(data.TLE == 0);name{1} = 'Con';% for 'TLE_Con_working_090824.xlsx'
% % ETLE
% %idx{2} = find(strcmp(data.Focal,'F') | strcmp(data.Focal,'P')); name{2} = 'ETLE';
% % % % TLE
% idx{2} = find(strcmp(data.Focal,'T')); name{2} = 'TLE';
% % % % all Focal
% % % %idx{2} = find(~strcmp(data.Focal,'N')); name{2} = 'all Focal';
% % % % FLE
% idx{3} = find(strcmp(data.Focal,'F')); name{3} = 'FLE';
% % % % PQE
% idx{4} = find(strcmp(data.Focal,'P')); name{4} = 'PQE';

% % Mts vs. no mts
idx{2} = find(data.mts == 1 & strcmp(data.Focal,'T'));name{2} = 'TLE-HS';
idx{3} = find(data.mts == 0 & strcmp(data.Focal,'T')); name{3} = 'TLE-other';

% %1.5T patients
% idx{2} = find((data.TLE == 1 & ~strcmp(data.Scanner, 'Prisma')));
% name{2} = '1.5TP';
% % 3T patients
% idx{2} = find((data.TLE == 1 & strcmp(data.Scanner, 'Prisma')));
% name{2} = '3TP';

% Overall patients(1.5T+3T)
%idx{2} = find(data.TLE == 1);name{2} = '1.5T+3T P';
% idx{2} = find((data.TLE == 1& strcmp(data.mri_side, 'left')));name{2} = '1.5T+3T left TLE';
% idx{2} = find((data.TLE == 1& strcmp(data.mri_side, 'right')));name{2} = '1.5T+3T right TLE';

%  No sf vs. sf
%idx{3} = find(strcmp(data.sf_latest,'0'));name{3} = 'nosf';
%idx{2} = find(strcmp(data.sf_latest,'1'));name{2} = 'sf';
%idx{3} = find(strcmp(data.sf_latest,'0')&strcmp(data.Focal,'T'));name{3} = 'nosf';
%idx{2} = find(strcmp(data.sf_latest,'1')&strcmp(data.Focal,'T'));name{2} = 'sf';
% idx{3} = find(strcmp(data.sf_latest,'0')&strcmp(data.Focal,'F'));name{3} = 'nosf';
% idx{2} = find(strcmp(data.sf_latest,'1')&strcmp(data.Focal,'F'));name{2} = 'sf';
% idx{3} = find(strcmp(data.sf_latest,'0')&data.mts == 1);name{3} = 'nosf';
% idx{2} = find(strcmp(data.sf_latest,'1')&data.mts == 1);name{2} = 'sf';

% %%No sf vs. sf, 1.5T, ATLR
% idx{2} = find((data.sf_latest == 0 & data.ATLR==1 & ~strcmp(data.Scanner, 'Prisma')));
% name{2} = 'nosf-ATLR';
% idx{3} = find((data.sf_latest == 1 & data.ATLR==1 & ~strcmp(data.Scanner, 'Prisma')));
% name{3} = 'sf-ATLR';

% %fbtcs vs. not-fbtcs
% idx{2} = find((strcmp(data.types_focal_to_bilateral_tonic_clonic,'yes')));
% name{2} = 'FBCTS';
% idx{3} = find((strcmp(data.types_focal_to_bilateral_tonic_clonic,'no')));
% name{3} = 'Not-FBCTS';
%fbtcs vs. not-fbtcs,1.5T
% idx{2} = find((strcmp(data.types_focal_to_bilateral_tonic_clonic,'yes')& ~strcmp(data.Scanner, 'Prisma')));
% name{2} = 'FBCTS';
% idx{3} = find((strcmp(data.types_focal_to_bilateral_tonic_clonic,'no')& ~strcmp(data.Scanner, 'Prisma')));
% name{3} = 'Not-FBCTS';


for g = 1:numel(name)
    nosub{g} = length(idx{g}); 
    age{g} = data.age_fmri_years(idx{g});
    sex{g} = data.sex(idx{g});
    epilepsy_duration{g} = data.duration_epilepsy_fmri_years(idx{g});
    %tiv{g}=data.tiv_fs(idx{g});
    %l_thalamus_vol{g}=data.l_THALAMUS(idx{g});
end

disp([name nosub]);
disp('------------');

% Load graph metrics (local measure)
clear all_grp* roi BG thal mtemp ltemp
%cd 'C:\Users\skgtxfe\OneDrive - University College London\Study2MultimodelThal\Focal_Con\hubness\pre-processed_data_TS\filtered_with_thalamic_subregion\Harv-ox_thal8nuclei\';
cd '/Users/fxy/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Study2MultimodelThal/Focal_Con/hubness/pre-processed_data_TS/filtered_with_thalamic_subregion/Harv-ox_thalnuclei_native/';
%load Harv-oxthalnuclei_all_wei_bin_harmo.mat % for 'TLE_Con_working_090824.xlsx'
load Harv-oxthalnuclei_all_wei_bin_harmo_TLE_ETLE_Con.mat
nw='Harv-oxthalnuclei_harmo';

load brainregionFeng.mat
load brainregionFengidx.mat
roi=brainregionFeng;

den_switch=2;
%den_thr = 10:50; % num_step=41
den_thr = 5:50; % num_step=46
%den_thr = 1:99; % num_step=99

% Load ROIs
clear hub_roi titles
% % ROIs=Whole brain nodes 
% titles = roi(:);
% hub_roi=1:length(roi);  

% ROIs=predefined nodes % NOTE:use correct brainregionidx of the atlas!
% % AAL
% thal=[91 92];
% thalnuclei=[91 92 93 94 95 96 97 98];
% mtemp=[37 38 39 40 41 42 89 90];
% ltemp=[79 80 81 82 83 84 85 86 87 88];
% BG=[71 72 73 74 75 76];

% % Harv-ox
thal=[107 108];
mtemp=[62 63 64 65 98 99 100 101 105 106];
ltemp=[17 18 19 20 21 22 23 24 27 28 29 30 15 16];
BG=[92 93 94 95 96 97 102 103]; %19/08/24: remove 104 brainstem
thalnuclei=length(roi) - 7:length(roi); % 4 nuclei groups
%thalnuclei=107:122; % 8nuclei

 %hub_roi=[thal mtemp BG];
 hub_roi=thalnuclei;
 titles = roi(hub_roi,:);

%% Step 2.1 uncorrected thalamic nuclei mapping
clear all_grp_nw all_grp_nw_z
% 1.extract 50% density nw (den_idx=46)
for g = 1:numel(name)
    for ii = 1:length(den_thr)
       % % Binary graph   
       %  graph='bin';
       %  all_grp_nw{g}(:,:) = squeeze(mean(W_ot_bin(:,:,idx{g},46),3)); 

        %Weight graph   
        graph='wei';
        all_grp_nw{g}(:,:) = squeeze(mean(W_ot(:,:,idx{g},46),3)); % get the group mean connectiivty matrix
    end
end

% % 2.1 z-score per nuclei
% for g = 1:numel(name)
%     tmp_z=zscore(all_grp_nw{g}(:,:),[],2);
%     all_grp_nw_z{g}(:,:)=tmp_z(hub_roi(:),:)';%selected nuclei
% end

% 2.2 Normalize per nuclei without changing sign(perverse + & - direction
% using max absolute value scaling to [-1,1]
for g = 1:numel(name)
    tmp_max_abs = max(abs(all_grp_nw{g}), [], 2); % Get max absolute value per row
    tmp_max_abs(tmp_max_abs == 0) = 1; % Avoid division by zero
    tmp_norm = all_grp_nw{g} ./ tmp_max_abs; % Scale by max absolute value
    all_grp_nw_z{g}(:,:) = tmp_norm(hub_roi(:), :)'; % Selected nuclei
end

% Control: mean of left and right
% % i+1
% for i=1:2:length(hub_roi) 
%     clear tmp*
%     tmp_1 = squeeze(all_grp_nw_z{1}(:,i)); 
%     tmp_2 = squeeze(all_grp_nw_z{1}(:,i+1)); 
%     tmp_mean = mean([tmp_2,tmp_1],2);
%     all_grp_nw_z{1}(:,i)=tmp_mean;
%     all_grp_nw_z{1}(:,i+1)=tmp_mean; 
% end  
% OR i+8
for i=1:length(hub_roi)/2 
    clear tmp*
    tmp_1 = squeeze(all_grp_nw_z{1}(:,i)); 
    tmp_2 = squeeze(all_grp_nw_z{1}(:,i+8)); %i+1 or i+8
    tmp_mean = mean([tmp_2,tmp_1],2);
    all_grp_nw_z{1}(:,i)=tmp_mean;
    all_grp_nw_z{1}(:,i+8)=tmp_mean; 
end 

%% Step 2.2.1 thalamic nuclei mapping with age, sex corrected, using Method 2
%all_grp_nw{g}: 3D matrix for group g (other ROIs × hub ROI × subjects).
%age{g}, sex{g}: Vectors of age and sex for each subject in group g.
%X_g * beta: Predicts the adjusted mean connectivity for group g at reference age and sex.

clear all_grp_n* Z* t_grp

% 1. Extract 50% density network (den_idx = 46)
for g = 1:numel(name)
    % Extract weight graph: other ROIs * hub ROI * subjects, at density index 46
    all_grp_nw{g} = squeeze(W_ot(setdiff(1:length(roi), hub_roi), hub_roi, idx{g}, 46));

    %unadjusted mean
    all_grp_nw_unadj{g}(:,:) = squeeze(mean(W_ot(setdiff(1:length(roi), hub_roi), hub_roi,idx{g},46),3)); % get the group mean connectiivty matrix

end

% Initialize cell array 
for g = 1:numel(name)
    t_grp{g} = nan(size(all_grp_nw{1}, 1), size(all_grp_nw{1}, 2));
    all_grp_nw_adj{g}=nan(size(all_grp_nw{1}, 1), size(all_grp_nw{1}, 2));
end

% 2. Adjust for age, sex. Get the group mean

% % Regressor for Method 1  
% sex_categorical_c = categorical(sex{1});% {1} = controls, {2} = patients
% sex_categorical_p = categorical(sex{2});

% Regressor for Method 2
age_all = [];
sex_all = [];
group_all = [];

for g = 1:numel(name)
    n_subj = length(idx{g});
    age_all = [age_all; age{g}];
    sex_all = [sex_all; sex{g}];
    group_all = [group_all; repmat(name{g}, n_subj, 1)];
end

sex_categorical_all = categorical(sex_all);
group_categorical_all = categorical(cellstr(group_all));

mean_age_ref = mean(age_all);

mean_sex_ref = 0.5;% assuming half M & half F
%mean_sex_ref = mean(double(categorical(sex_all)) - 1); % the proportion of M in the overal sample(con + pat) (Assumes binary coding: 0 = F, 1 = M)

% Loop through each region 
for j = 1:size(all_grp_nw{1}, 1) %j=other ROIs 
    for x = 1:size(all_grp_nw{1}, 2) %x=hub ROI
  
        % % Method 1: Fit GLM using controls only 
        %Y_c = squeeze(all_grp_nw{1}(j, x, :)); 
        %Y_p = squeeze(all_grp_nw{2}(j, x, :));  

        % control_table = table(Y_c, age{1}, sex_categorical_c, 'VariableNames', {'Y', 'Age', 'Sex'});
        % control_glm = fitglm(control_table, 'Y ~ Age + Sex', 'Distribution', 'normal');      
        % 
        % % Get residuals (adjusted connectivity) for both controls and patients
        % X_c = [ones(length(age{1}), 1), age{1}, double(sex_categorical_c)-1]; % Controls
        % X_p = [ones(length(age{2}), 1), age{2}, double(sex_categorical_p)-1]; % Patients
        % 
        % all_grp_nw_adj{1}(j, x, :) = Y_c - X_c * control_glm.Coefficients.Estimate; 
        % all_grp_nw_adj{2}(j, x, :) = Y_p - X_p * control_glm.Coefficients.Estimate;
        % % Methods 1 cannot get group mean because the control's residuals are designed to have a mean =0

        % Method 2: Fit combined model
        % Y for all groups
        Y_all=[];
        for g = 1:numel(name)
            Y_g = squeeze(all_grp_nw{g}(j, x, :)); 
            Y_all = [Y_all; Y_g];  % concatenate across subjects
        end

        tbl = table(Y_all, age_all, sex_categorical_all, group_categorical_all, ...
            'VariableNames', {'Y', 'Age', 'Sex', 'Group'});
        glm = fitglm(tbl, 'Y ~ Age + Sex + Group', 'Distribution', 'normal');
        beta = glm.Coefficients.Estimate;
        coef_names = glm.CoefficientNames; % check glm.CoefficientNames

        % For each group, create X vector
        for g = 1:numel(name)
            X_g = zeros(1, length(beta));
            X_g(strcmp(coef_names, '(Intercept)')) = 1;
            X_g(strcmp(coef_names, 'Age')) = mean_age_ref;
            X_g(strcmp(coef_names, 'Sex_M')) = mean_sex_ref;      
            
            % Group dummy variables
            grp_label = ['Group_' name{g}]; % e.g., Group_TLE

            if ismember(grp_label, coef_names)
                idx_coef = find(strcmp(coef_names, grp_label));
                X_g(idx_coef) = 1;

                % Extract t-statistics for group difference per connection(positive means patient > control, negative means patient < control)
                t_grp{g}(j, x) = glm.Coefficients.tStat(idx_coef);
            else
                t_grp{g}(j, x) = NaN; % Controls as the reference group
            end
            
            %  adjusted group mean
            all_grp_nw_adj{g}(j, x) = X_g * beta;

        end
    end
end

%% Step 2.2.2 thalamic nuclei mapping with age, sex corrected, using Method 1&2 conbimed
% adjust using beta for control; residual for patients
% independent t tests for patient-control difference
%07072025 Feng: currently used for the paper

clear all_grp_n* Z* t_grp

% 1. Extract 50% density network (den_idx = 46)
for g = 1:numel(name)
    % Extract weight graph: other ROIs * hub ROI * subjects, at density index 46
    all_grp_nw{g} = squeeze(W_ot(setdiff(1:length(roi), hub_roi), hub_roi, idx{g}, 46));

    %unadjusted mean
    all_grp_nw_unadj{g}(:,:) = squeeze(mean(W_ot(setdiff(1:length(roi), hub_roi), hub_roi,idx{g},46),3)); % get the group mean connectiivty matrix

end

% Initialize cell array 
for g = 1:numel(name)
    all_grp_nw_adj{g}=nan(size(all_grp_nw{1}, 1), size(all_grp_nw{1}, 2));
end

% 2. Adjust for age, sex. Get the group mean

% % Regressor for Method 1  
sex_categorical_c = categorical(sex{1});% {1} = controls

% Regressor for Method 2
age_all = [];
for g = 1:numel(name)
    age_all = [age_all; age{g}];
end
mean_age_ref = mean(age_all);

mean_sex_ref = 0.5;% assuming half M & half F
%mean_sex_ref = mean(double(categorical(sex_all)) - 1); % the proportion of M in the overal sample(con + pat) (Assumes binary coding: 0 = F, 1 = M)

% Loop through each region 
for j = 1:size(all_grp_nw{1}, 1) %j=other ROIs 
    for x = 1:size(all_grp_nw{1}, 2) %x=hub ROI
  
        % Fit GLM using controls only (Method 1)
        Y_c = squeeze(all_grp_nw{1}(j, x, :));         
        control_table = table(Y_c, age{1}, sex_categorical_c, 'VariableNames', {'Y', 'Age', 'Sex'});
        control_glm = fitglm(control_table, 'Y ~ Age + Sex', 'Distribution', 'normal');  
        beta = control_glm.Coefficients.Estimate;
        coef_names = control_glm.CoefficientNames; % check glm.CoefficientNames

        % Use residuals to test group difference
        X_c = [ones(length(age{1}), 1), age{1}, double(sex_categorical_c)-1]; % Controls
        residuals_c = Y_c - X_c * beta; 

        for g = 2:numel(name)   % each patient group
            clear residuals_p Y_p X_p stats
            Y_p = squeeze(all_grp_nw{g}(j, x, :));  
            X_p = [ones(length(age{g}), 1), age{g}, double(categorical(sex{g}))-1]; 
            residuals_p = Y_p - X_p * beta;

            % % Option1: compare with Control
            % [~, p, ~, stats] = ttest2(residuals_p, residuals_c, 'Vartype', 'unequal'); % Welch's t-tests for unequal variance
            % all_grp_nw_adj{g}(j, x) = stats.tstat; % t>0: patient>control. t<0: patient<control
            % %all_grp_nw_adj{g}(j, x) = p; % save uncorrect p value

            % Option2: compare two patient grps g2&g3 (tmp(2)&tmp(3))
            tmp{g}=residuals_p;
        end
        % Continue Option2: compare two patient grps
        [~, p, ~, stats] = ttest2(tmp{2}, tmp{3}, 'Vartype', 'unequal'); % Welch's t-tests for unequal variance
        all_grp_nw_adj{2}(j, x) = stats.tstat; all_grp_nw_adj{3}(j, x) = stats.tstat; % t>0: g2>g3. t<0: g2<g3, results saved in both g2 and g3 location
        %all_grp_nw_adj{2}(j, x) = p; all_grp_nw_adj{3}(j, x) = p; % uncorrect p value

        % Use beta for normative control connectivity profile (Method 2)
        X_g = zeros(1, length(beta));
        X_g(strcmp(coef_names, '(Intercept)')) = 1;
        X_g(strcmp(coef_names, 'Age')) = mean_age_ref;
        X_g(strcmp(coef_names, 'Sex_M')) = mean_sex_ref;      
                  
        all_grp_nw_adj{1}(j, x) = X_g * beta;%  adjusted control group mean
   
    end
end

% all_grp_nw_adj{1} is normative control connectivity; all_grp_nw_adj{2:end} is patient minus control t value

%% fdr correction of p value
for g = 2:numel(name)   % each patient group
    for x = 1:size(all_grp_nw{1}, 2) %x=hub ROI
    P = all_grp_nw_adj{g}(:, x); % uncorrect p value
    fdrP = mafdr(P', 'BHFDR', true);
    all_grp_nw_adj{g}(:, x)=fdrP'; 
    disp(min(fdrP))
    end
end

%% follow Step 2.2.1 or 2.2.2 

% For controls: average L/R 

% remaining regions without left and right subparts
special_idx = [49 52 55 56 57 104]; % idx for nodes without left and right subparts
% remaining region average left and right
pair_idx = setdiff(1:size(all_grp_nw{1}, 1), special_idx);
for i=1:2:length(pair_idx) % use mean of left and right
    x=pair_idx(i);
    y=pair_idx(i+1);

    tmp_1 = all_grp_nw_adj{1}(x, :);
    tmp_2 = all_grp_nw_adj{1}(y, :);
    tmp_mean = mean([tmp_2;tmp_1]);

    all_grp_nw_adj{1}(x, :) = tmp_mean;
    all_grp_nw_adj{1}(y, :) = tmp_mean;
    clear tmp*
end  

% average across left and right thalamus 
% for 4 nuclei group
for i=1:2:length(hub_roi) 
    clear tmp*
    tmp_1 = squeeze(all_grp_nw_adj{1}(:,i)); 
    tmp_2 = squeeze(all_grp_nw_adj{1}(:,i+1)); 
    tmp_mean = mean([tmp_2,tmp_1],2);
    all_grp_nw_adj{1}(:,i)=tmp_mean;
    all_grp_nw_adj{1}(:,i+1)=tmp_mean;  
end  
% %OR for 8 nuclei
% for i=1:length(hub_roi)/2 
%     clear tmp*
%     tmp_1 = squeeze(all_grp_nw_adj{1}(:,i)); 
%     tmp_2 = squeeze(all_grp_nw_adj{1}(:,i+8)); %i+1 or i+8
%     tmp_mean = mean([tmp_2,tmp_1],2);
%     all_grp_nw_adj{1}(:,i)=tmp_mean;
%     all_grp_nw_adj{1}(:,i+8)=tmp_mean; 
% end 

% For patient: select ipsi/contra side
hub_roi_l=[1:2:length(hub_roi)];%idx for ipsi/left nuclei

for g = 1:numel(name)
    tmp_r=[]; tmp_l=[];tmp_r_t=[]; tmp_l_t=[];
    for j = 1:size(all_grp_nw{1}, 1) % loop through brain regions

        if endsWith(roi{j}, ' r') % controlateral(right hemisphere) % double check if region names end with ' r'!!
          tmp_r =[tmp_r;all_grp_nw_adj{g}(j,hub_roi_l+1)];
          %tmp_r_t =[tmp_r_t;t_grp{g}(j,hub_roi_l+1)]; % only use it when applying Method 2 solely
        end 

        if endsWith(roi{j}, ' l') % ipsilateral(left hemisphere) % double check if region names end with ' l'!!
          tmp_l =[tmp_l;all_grp_nw_adj{g}(j,hub_roi_l)]; 
          %tmp_l_t =[tmp_l_t;t_grp{g}(j,hub_roi_l)];
        end
    end 
    all_grp_nw_adj_r{g}=[tmp_r;all_grp_nw_adj{g}(special_idx,hub_roi_l+1)];% At the end, add regions without left and right subparts
    all_grp_nw_adj_l{g}=[tmp_l;all_grp_nw_adj{g}(special_idx,hub_roi_l)];

    %t_grp_r{g}=[tmp_r_t;t_grp{g}(special_idx,hub_roi_l+1)];
    %t_grp_l{g}=[tmp_l_t;t_grp{g}(special_idx,hub_roi_l)];
end

% copy all_grp_nw_adj_r/r=l to 'thalamicmapping_Focal_Con.xlsx', 'Sheet',{'controls_adj', 'TLE_adj', 'FLE_adj', 'PQE_adj'}
% copy t_grp_r/r=l{2}{3}{4} to 'thalamicmapping_Focal_Con.xlsx', 'Sheet',{'TLE_adj', 'FLE_adj', 'PQE_adj'} e.g., Ant l T

%% Option 1: manually set color and node size(may need to normalise) for Brainnewtviewer

%% Option 2: Simple brain plot, See below to manually match with my labels with Lausanne scale 1
%% Match Harv-ox to Lausanne and write to multiple sheets
clear
cd /Users/fxy/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Study2MultimodelThal/Feng_code/mapping

% Define group sheet names
groupNames = {'controls_adj', 'TLE_adj', 'FLE_adj', 'PQE_adj','HS_adj', 'NHS_adj','HS-NHS_adj'};

% Read Lausanne idx to match
Tabeltitle='thalamicmapping_Focal_Con_native_flip.xlsx';
roiList = readtable(Tabeltitle, 'Sheet', 'Lausanne_match_template');
roi_idx_list = roiList.roi_idx;

% Loop through each group
for g = 1:length(groupNames)
    groupName = groupNames{g};

    % Read original roi_idx for the current group
    mainTable = readtable(Tabeltitle, 'Sheet', groupName);
    main_idx = mainTable.roi_idx;

    matchedTable = table();
    
    % For each roi_idx in roi_idx_list, find the corresponding row in mainTable
    for i = 1:length(roi_idx_list)
        roi_id = roi_idx_list(i);

        % Find the first occurrence (or handle if multiple in mainTable)
        idx = find(main_idx == roi_id, 1, 'first');

        if ~isempty(idx)
            matchedTable = [matchedTable; mainTable(idx, :)];
        else
            % If not found, append a row of NaNs or missing data
            emptyRow = mainTable(1, :);
            emptyRow{:,:} = {NaN};  % Fill entire row with NaN
            matchedTable = [matchedTable; emptyRow];
        end
    end
    
    % Write to the corresponding sheet in the output Excel file
    writetable(matchedTable,'matched_roi.xlsx', 'Sheet', groupName);
end
% copy 'matched_roi.xlsx' to 'thalamicmapping_Focal_Con.xlsx', 'Sheet', 'Lausanne_match'

%% calculate average for patially matched Lausanne label
clear
cd /Users/fxy/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Study2MultimodelThal/Feng_code/mapping

% Read the main data table
Tabeltitle='thalamicmapping_Focal_Con_native_flip.xlsx';
mainTable = readtable(Tabeltitle, 'Sheet', 'Lausanne_match');

% Extract Lausanne_idx
lausanne_idx = mainTable.Lausanne_idx;

% Define the range of columns (R to AO)
varNames = mainTable.Properties.VariableNames;
%colStart = find(strcmp(varNames, 'AntControl'));
colStart = find(strcmp(varNames, 'AntLTHS_NHS'));
colEnd = find(strcmp(varNames, 'PostRT_5'));

% Get the data to average
dataToAverage = mainTable{:, colStart:colEnd};

% Get unique Lausanne_idx
uniqueIDs = unique(lausanne_idx);

% Pre-allocate result table
avgData = array2table(nan(length(uniqueIDs), colEnd - colStart + 1), ...
    'VariableNames', varNames(colStart:colEnd));
avgData.Lausanne_idx = uniqueIDs;

% Loop through each unique Lausanne_idx and calculate mean per column
for i = 1:length(uniqueIDs)
    id = uniqueIDs(i);
    rows = lausanne_idx == id;
    
    % For each column, compute mean across rows (omit NaN)
    avgValues = mean(dataToAverage(rows, :), 1, 'omitnan');
    
    avgData{i, 1:(colEnd - colStart + 1)} = avgValues;
end

% Move Lausanne_idx to the first column
avgData = movevars(avgData, 'Lausanne_idx', 'Before', 1);

% Write to a new Excel file
writetable(avgData, 'average_by_lausanne.xlsx');
% copy 'average_by_lausanne.xlsx' to 'thalamicmapping_Focal_Con.xlsx', 'Sheet', 'allgrp_Lausanne'

%% Use SimpleBrainPlot to visualise: https://github.com/dutchconnectomelab/Simple-Brain-Plot

%% only plot Controls
clear
cd /Users/fxy/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Study2MultimodelThal/Feng_code/mapping

T = readtable('thalamicmapping_Focal_Con_native_flip.xlsx', 'Sheet', 'allgrp_adj_Lausanne');

labels = T.LausanneLabel;

% Define columns for plotting
varNames = T.Properties.VariableNames; disp(varNames);

% colStart = find(strcmp(varNames, 'AntControl'));
% colEnd = find(strcmp(varNames, 'PostR_2'));
% colRange = colStart:colEnd;
colRange = 28:31;

% Color scheme 
cm = [1 1 1; 0.6 0.0 0.0];%Dark Red
cm = interp1(cm, 1:0.01:size(cm, 1));

%cm = flipud(gray(256)); % 256-step grayscale colormap (black to white)

% colors = [0 0.4470 0.7410;  % Blue
%           1 1 1;           % White
%           1 0 0];  % Red
% colors = [0.2 0.4 0.8;  % Blue
%           0.9 0.9 0.9;  % Light gray
%           0.9 0.3 0.1]; % Orange
% colors = [0.1 0.2 0.6;  % Dark Blue
%           1.0 1.0 1.0;  % White
%           0.6 0.0 0.0]; % Dark Red
% cm = interp1([0 0.5 1], colors, linspace(0, 1, 256));

cd simple_brain_plot/4nucleigrp_age_sex_adj_native_flip/allredp10p90
%cd simple_brain_plot/4nucleigrp_age_sex_adj/con/red_white_blue_con

% Loop over columns and plot
for c = colRange
    % Extract mapplot for this column
    mapplot = T{:, c};
    
    % Set limits 
    % Option 1:10th–90th percentile p10-p90, all red colorscale
    lims = prctile(mapplot, [10 90], 'all');

    % % Option 2:include both pos & neg, blue-red colorscale
    % lim_val = max(abs(prctile(mapplot, [10 90]))); %Symmetric lim to make sure white in the middle(zero)
    % %lim_val = max(abs(mapplot));
    % lims = [-lim_val, lim_val];
    
    plotBrain(labels, mapplot, cm, ...
        'atlas', 'lausanne120_aseg', ...
        'savePath', [varNames{c}], ...
        'limits', lims);
end

%% plot Patients & controls on same scale per nuclei
clear
cd /Users/fxy/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Study2MultimodelThal/Feng_code/mapping
T = readtable('thalamicmapping_Focal_Con_native_flip.xlsx', 'Sheet', 'allgrp_adj_Lausanne');

labels = T.LausanneLabel;
varNames = T.Properties.VariableNames;
disp(varNames);

% Color scheme
cm = [1 1 1; 0 0.4470 0.7410];
cm = interp1(cm, 1:0.01:size(cm, 1));

cd simple_brain_plot/4nucleigrp_age_sex_adj_native_flip
nucleusGroups = {'Ant', 'Lat', 'Med', 'Post'};

% Loop through each nucleus group
for n = 1:length(nucleusGroups)
    nucleus = nucleusGroups{n};

    % Find columns that start with this nucleus name
    cols_nucleus = find(startsWith(varNames, nucleus));

    % Combine data across all columns of this nucleus group
    combinedVals = [];
    for c = cols_nucleus
        combinedVals = [combinedVals; T{:, c}];
    end

    % Compute 10th–90th percentile across combined data
    lims = prctile(combinedVals, [10 90], 'all');

    % Loop over each column and plot with the same limits
    for c = cols_nucleus
        mapplot = T{:, c};

       % Plot using SimpleBrainPlot
        plotBrain(labels, mapplot, cm, ...
            'atlas', 'lausanne120_aseg', ...
            'savePath', [varNames{c}], ...
            'limits', lims);
    end
end

%% Plot patient-control difference
clear
cd /Users/fxy/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Study2MultimodelThal/Feng_code/mapping
T = readtable('thalamicmapping_Focal_Con_native_flip.xlsx', 'Sheet', 'allgrp_adj_Lausanne');

labels = T.LausanneLabel;
varNames = T.Properties.VariableNames;

% Custom colormap: Blue–white–red
colors = [0 0.4470 0.7410; 1 1 1; 1 0 0];  % blue;white;red
cm = interp1([0 0.5 1], colors, linspace(0, 1, 256));

%% --- 1. Plot patient-control difference's t-stats 
% Define nucleus keywords (as prefixes)
nucleusPrefixes = {'Ant', 'Lat', 'Med', 'Post'};

% Compute color bar limits across all columns
combinedVals = [];
for n = 1:length(nucleusPrefixes)
    prefix = nucleusPrefixes{n};

    % Get columns that contain both nucleus prefix and 'tstats'
    cols = find(contains(varNames, prefix) & contains(varNames, 'Tstats'));

    for c = cols(:)'
        vals = T{:, c};
        combinedVals = [combinedVals; vals];
    end
end
%lim_val = max(abs(prctile(combinedVals, [5 95])));% show 5-95%
% lim_val = max(abs(combinedVals));% show full range
% lims = [-lim_val, lim_val];

lims = [-3.5, 3.5];% t stat >= 3.5 to show approx significant p<0.001

%cd simple_brain_plot/4nucleigrp_age_sex_adj/pat_con_difference/t-stats/t-stats3/nsf_sf/
cd simple_brain_plot/4nucleigrp_age_sex_adj_native_flip

for n = 1:length(nucleusPrefixes)
    prefix = nucleusPrefixes{n};
    %cols = find(contains(varNames, prefix) & contains(varNames, 'Tstats'));
    %cols = find(contains(varNames, prefix) & contains(varNames, 'sf_nsf'));
    cols = find(contains(varNames, prefix) & contains(varNames, 'HS_NHSTstats')); % plot noHS-HS

    for c = cols(:)'
        mapplot = T{:, c};

        %mapplot = mapplot*-1;% plot reversed test results
              
        % % Mask non-significant results: keep only |t| ≥ 2 
        % mapplot(abs(mapplot) < 2) = 0;     
        % if all(mapplot == 0) % Check if all values are zero
        %     custom_cm = repmat([1 1 1], 256, 1);  % All white colormap
        % else
        %     custom_cm = cm; % Use default (e.g., blue-white-red)
        % end

        custom_cm = cm; % Use default (e.g., blue-white-red)

        plotBrain(labels, mapplot, custom_cm, ...
            'atlas', 'lausanne120_aseg', ...
            'savePath', varNames{c}, ...
            'limits', lims);
    end
end

%% --- 2. Plot patient-control difference in group mean
cd simple_brain_plot/4nucleigrp_age_sex_adj/pat_con_difference/mean/
nucleusGroups = {'Ant', 'Lat', 'Med', 'Post'};

% Compute 5th–95th percentile across all nuclei to plot with the same limits
combinedVals = [];
for n = 1:length(nucleusGroups) % loop through each nuclei
    nucleus = nucleusGroups{n};

    % Find columns that start with this nucleus name
    cols_nucleus = find(startsWith(varNames, nucleus));

    % Combine data across all columns of nucleus groups
    for c = cols_nucleus(2:end)
        mapplot = T{:, c}-T{:, cols_nucleus(1)};% Difference=patient MINUS control
        combinedVals = [combinedVals; mapplot];
    end
end
% Symmetric percentile limits to make sure white in the middle(zero)
lim_val = max(abs(prctile(combinedVals, [10 90])));
lims = [-lim_val, lim_val];

% Loop through each nucleus group
for n = 1:length(nucleusGroups)
    nucleus = nucleusGroups{n};

    % Find columns that start with this nucleus name
    cols_nucleus = find(startsWith(varNames, nucleus));

    % % Compute 10th–90th percentile per nuclei to plot with the same limits
    % combinedVals = [];
    % for c = cols_nucleus(2:end)
    %     mapplot = T{:, c}-T{:, cols_nucleus(1)};% Difference=patient MINUS control
    %     combinedVals = [combinedVals; mapplot];
    % end
    % 
    % %lims = prctile(combinedVals, [10 90], 'all');
    % % Symmetric percentile limits
    % lim_val = max(abs(prctile(combinedVals, [10 90])));
    % lims = [-lim_val, lim_val];

    for c = cols_nucleus(2:end)
        mapplot = T{:, c}-T{:, cols_nucleus(1)}; % Difference=patient MINUS control

       % Plot using SimpleBrainPlot
        plotBrain(labels, mapplot, cm, ...
            'atlas', 'lausanne120_aseg', ...
            'savePath', [varNames{c}], ...
            'limits', lims);
    end
end

%% Step 2.3 TBC struc &func mapping for 8 nuclei ("AV""VA""VLa""VLP""VPL""Pul""CM""MD_Pf")
% 13 May Feng: select LEFT nuclei
% corrected for age, sex and tiv
% normalise to make these two comparable

% 13 May Feng: func subj all have a larger tiv than struc subj
% -> unable to adjust tiv through regression model due to collinearity between tiv and Group('struc', 'func')
% ? alternatively, use left whole thalamus volume=sum of vol("AV""VA""VLa""VLP""VPL""Pul""CM""MD_Pf")

%% Load control's functional data
clear labels all_grp* 
cd('/Users/fxy/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Study2MultimodelThal/Feng_code/mapping/simple_brain_plot/');

% Extract connectome: other ROIs * hub ROI(LEFT) * subjects, at density index 46
all_grp_nw{1} = squeeze(W_ot(setdiff(1:length(roi), hub_roi), hub_roi(1:8), idx{1}, 46));

% Transfer Harv-ox atlas to Lausanne scale 1

% load Lausanne match table
lausanne_match = readtable('LausanneScale1_new_thomas_labels_thalamus.xlsx', Sheet='Lausanne_match_template');
lausanne_idx = lausanne_match.Lausanne_idx_rory;
unique_lausanne = unique(lausanne_idx); 

harvox_idx = lausanne_match.roi_idx;

all_grp_func{1} = zeros(length(unique_lausanne), size(all_grp_nw{1}, 2), size(all_grp_nw{1}, 3)); % unique_lausanne * num_hubROIs * num_subjects

% Loop through each unique Lausanne index
for i = 1:length(unique_lausanne)

    % Find all harvox_idx corresponding to this Lausanne index
    matching_rows = find(lausanne_idx == unique_lausanne(i));
    source_rows = harvox_idx(matching_rows);

    % Compute mean across these rows
    all_grp_func{1}(i, :, :) = mean(all_grp_nw{1}(source_rows, :, :), 1);
end

%% load control's structural data 
load('structural_RP/feng.mat');% from feng.mat fron RP

clear map maps
for b = 1:length(thalamus_generic_labels)
    clear map
    for a = 1:length(controls)
        map(a,:) = controls(a).connectome(selected_nodes(b+8),setdiff(1:length(labels), selected_nodes));% % select LEFT nuclei 
    end
    % Extract connectome: hub ROI(LEFT) * subjects * other ROIs
    maps(b,:,:) = map; 
end

% reshape to be in same dimension as all_grp_func
all_grp_struc{1} = permute(maps, [3, 1, 2]); % other ROIs * hub ROI(LEFT) * subjects

%% Adjust for age, sex, tiv. Get the group mean

% set X & Y for regression models
age_func=age{1};
age_struc = [controls.age]';
age_all = [age_func; age_struc]; % concatenate across subjects

% tiv_func=single(tiv{1});
% tiv_struc=[controls.tiv_sum]';

% Use left whole thalamus volume instead of tiv 
tiv_func= sum([data.l_AV(idx{1}),data.l_VA(idx{1}),data.l_VLa(idx{1}),data.l_VLP(idx{1}),data.l_VPL(idx{1}),data.l_Pul(idx{1}),data.l_CM(idx{1}),data.l_MD_Pf(idx{1})],2);
for i = 1:length(controls)
    tiv_struc(i) = sum(controls(i).vols_1(selected_nodes(9:end)));
end
tiv_struc = tiv_struc';
tiv_all = [tiv_func; tiv_struc];

sex_func=sex{1};
sex_struc=cellstr([controls.sex]');
sex_all = [sex_func; sex_struc];
sex_categorical_all = categorical(sex_all);

group_func=cellstr(repmat('func_con', size(all_grp_func{1},3), 1));
group_struc=cellstr(repmat('struc_con', size(all_grp_struc{1},3), 1));
group_all = [group_func; group_struc];
group_categorical_all = categorical(cellstr(group_all));

Y_all = cat(3, squeeze(all_grp_func{1}), squeeze(all_grp_struc{1})); % (roi x roi x subjects) stacking along 3rd dim

% set average co-var for group mean
mean_age_ref = mean(age_all);
mean_sex_ref = 0.5;% assuming half M & half F
%mean_sex_ref = mean(double(categorical(sex_all)) - 1); % the proportion of M in the overal sample
mean_tiv_ref = mean(tiv_all);

% adjust using regression models
% model 1
% for j = 1:size(all_grp_func{1}, 1) %j=other ROIs 
%     for x = 1:size(all_grp_func{1}, 2) %x=thalamic ROI
% 
%         % Method 2: Fit combined model
%         Y_tmp=squeeze(Y_all(j, x, :));
%         tbl = table(Y_tmp, age_all, tiv_all, sex_categorical_all, group_categorical_all, ...
%             'VariableNames', {'Y', 'Age', 'TIV', 'Sex', 'Group'});
%         glm = fitglm(tbl, 'Y ~ Age + TIV + Sex + Group', 'Distribution', 'normal');
%         beta = glm.Coefficients.Estimate;
%         % check if glm.CoefficientNames in correct order for the following X
%  
%         % For each group, create X vector
%         % {'Group_struc_con'} is the reference group in coef_names, so:
%         X_func=[1,mean_age_ref,mean_tiv_ref,mean_sex_ref,0]; % 'Group'=0
%         X_struc=[1,mean_age_ref,mean_tiv_ref,mean_sex_ref,1];% 'Group'=1
% 
%         %  Get adjusted group mean
%         all_grp_func_adj{1}(j, x) = X_func * beta;
%         all_grp_struc_adj{1}(j, x) = X_struc * beta;
%     end
% end

% model 2: add interaction terms for different age/sex/tiv effect in each group
for j = 1:size(all_grp_func{1}, 1) % j = other ROIs
    for x = 1:size(all_grp_func{1}, 2) % x = thalamic ROI

        % Extract the connectivity values for this ROI and thalamic ROI
        Y_tmp = squeeze(Y_all(j, x, :));

        % Create a table with the predictors and outcome variable
        tbl = table(Y_tmp, age_all, tiv_all, sex_categorical_all, group_categorical_all, ...
            'VariableNames', {'Y', 'Age', 'TIV', 'Sex', 'Group'});

        % Fit the GLM with interaction terms
        glm = fitglm(tbl, 'Y ~ Age * Group + TIV * Group + Sex * Group', 'Distribution', 'normal');
        
        % Get the coefficient estimates
        beta = glm.Coefficients.Estimate;
        
        % Check the coefficient names 
        %disp(glm.CoefficientNames); 

        % For each group, create X vectors
        % {'Group_struc_con'} is the reference group in coef_names, so: For 'func_con' (Group = 0)
        X_func = [1, mean_age_ref, mean_tiv_ref, mean_sex_ref, 0, mean_age_ref*0, mean_tiv_ref*0, mean_sex_ref*0]; % Interactions with 0
        
        % For 'struc_con' (Group = 1)
        X_struc = [1, mean_age_ref, mean_tiv_ref, mean_sex_ref, 1, mean_age_ref*1, mean_tiv_ref*1, mean_sex_ref*1]; % Interactions with 1

        % Get adjusted group mean for func and struc
        all_grp_func_adj{1}(j, x) = X_func * beta;
        all_grp_struc_adj{1}(j, x) = X_struc * beta;

    end
end

%% plot struc & func seperatly

% Color scheme
cm = [1 1 1; 0 0.4470 0.7410];
cm = interp1(cm, 1:0.01:size(cm, 1));

% output figure folder
savefolder='8nuclei_adj/';
mkdir savefolder

% output figure names
Names = selected_nodes_names(9:end); % left nuclei

for b = 1:size(all_grp_func_adj{1},2)

    % %plot func
    mapplot_func = all_grp_func_adj{1}(:,b);
    
    lims_func = prctile(mapplot_func,[10 90],"all"); %percentile 10-90%

    plotBrain(labels(setdiff(1:length(labels), selected_nodes)), mapplot_func, cm, ...
            'atlas', 'lausanne120_aseg', ...
            'savePath', [savefolder Names{b} '_func'], ...
            'limits', lims_func);

    %plot struc
    mapplot_struc = all_grp_struc_adj{1}(:,b);

    lims_struc = prctile(mapplot_struc,[10 90],"all"); 

    plotBrain(labels(setdiff(1:length(labels), selected_nodes)), mapplot_struc, cm, ...
            'atlas', 'lausanne120_aseg', ...
            'savePath', [savefolder Names{b} '_struc'], ...
            'limits', lims_struc);
end

%% plot struc-func consistent connectivity patterns
% 1. ?? Select only ipsilateral connections (e.g., left side nodes for left thalamic nuclei): 
% 2. For each dataset, divide regions into 4 subgroups (100/75/50/25%) based on quantiles of connectivity strength.
% 3. For regions in the same subgroup for both struc & func, assign a value: 4/3/2/1 for 100/75/50/25% subgroup.
% 4. Assign 0 for regions not in the same subgroup.

% Color scheme: 
cm_red = [1 1 1; 1 0.6 0.6]; %light red
cm_red = interp1(cm_red, 1:0.01:size(cm_red, 1));

savefolder = '8nuclei_adj/consistency/';

Names = selected_nodes_names(9:end); % left nuclei

% TBC Select Only ipsilateral nodes

for b = 1:size(all_grp_func_adj{1},2)

    mapplot_func = all_grp_func_adj{1}(:,b);
    mapplot_struc = all_grp_struc_adj{1}(:,b);

    % Quantile thresholds (3 quantiles: 25, 50, 75%)
    q_func = quantile(mapplot_func, [0.25 0.5 0.75]);
    q_struc = quantile(mapplot_struc, [0.25 0.5 0.75]);

    % Assign subgroup for each region (1=lowest, 4=highest)
    q_func_group = 1 + (mapplot_func > q_func(1)) + (mapplot_func > q_func(2)) + (mapplot_func > q_func(3));
    q_struc_group = 1 + (mapplot_struc > q_struc(1)) + (mapplot_struc > q_struc(2)) + (mapplot_struc > q_struc(3));

    % Find regions in the same subgroup for both datasets
    same_group = (q_func_group == q_struc_group);

    % Representative value: 4/3/2/1 for 100/75/50/25% subgroup, else 0
    rep_conn = zeros(size(mapplot_func));
    rep_conn(same_group) = q_func_group(same_group);

    plotBrain(labels(setdiff(1:length(labels), selected_nodes)), rep_conn, cm_red, ...
        'atlas', 'lausanne120_aseg', ...
        'savePath', [savefolder Names{b}], ...
        'limits', [0 4]);
end

%% plot struc-func consistent regions, intrinsic connectivity patterns of each group not perserved
% 1. Min-max normalization: For each dataset, subtract the minimum value and divide by the range (max-min), scaling all values between 0 and 1.

% 2. Consistency=1−∣Si−Fi∣where Si and Fi are the normalized structural and functional values for region i 
% This index will be highest =1 when the values are identical and decrease as the difference increases.

% Color scheme: 
cm = [1 1 1; 0 0.4470 0.7410]; %blue
cm = interp1(cm, 1:0.01:size(cm, 1));

cm_orange = [1 1 1; 0.91 0.41 0.17];% orange
cm_orange = interp1(cm_orange, 1:0.01:size(cm_orange, 1));

% Output figure folder
savefolder = '8nuclei_adj/';

% Output figure names
Names = selected_nodes_names(9:end); % left nuclei

% Preallocate for consistency index
consistency = zeros(size(all_grp_func_adj{1}));

% Min-max normalization for each column (nuclei) independently
norm_func = normalize(all_grp_func_adj{1},"range");   % [subjects x regions]
norm_struc = normalize(all_grp_struc_adj{1},"range"); % [subjects x regions]

for b = 1:size(all_grp_func_adj{1},2)

    mapplot_func = norm_func(:,b);
    mapplot_struc = norm_struc(:,b);

    % Consistency index 
    mapplot_consistency = ones(size(mapplot_func)) - abs(mapplot_func - mapplot_struc);

    % Plot normalized func
    lims_func = prctile(mapplot_func, [10 90], "all");
    plotBrain(labels(setdiff(1:length(labels), selected_nodes)), mapplot_func, cm, ...
        'atlas', 'lausanne120_aseg', ...
        'savePath', [savefolder 'norm01/' Names{b} '_func_norm01'], ...
        'limits', lims_func);

    % Plot normalized struc
    lims_struc = prctile(mapplot_struc, [10 90], "all");
    plotBrain(labels(setdiff(1:length(labels), selected_nodes)), mapplot_struc, cm, ...
        'atlas', 'lausanne120_aseg', ...
        'savePath', [savefolder 'norm01/' Names{b} '_struc_norm01'], ...
        'limits', lims_struc);

    % Plot consistency
    lims_consistency = prctile(mapplot_consistency, [10 90], "all");
    plotBrain(labels(setdiff(1:length(labels), selected_nodes)), mapplot_consistency, cm_orange, ...
        'atlas', 'lausanne120_aseg', ...
        'savePath', [savefolder 'consistency_func-struc/' Names{b} '_consistency'], ...
        'limits', lims_consistency);
end

%% Option: check collinearity 
% Ensure numeric format
age_numeric = double(age_all);
tiv_numeric = double(tiv_all);
sex_numeric = double(categorical(sex_all)) - 1; % F=0, M=1
group_numeric = double(categorical(group_all)) - 1; % struc_con=1 (reference), func_con=0

% Combine into matrix
X = [age_numeric, tiv_numeric, sex_numeric, group_numeric];

% Compute correlation matrix
clear corr
corr_matrix = corr(X);

% Display with labels
labels = {'Age', 'TIV', 'Sex(M=1)', 'Group(pat=1)'};
disp(array2table(corr_matrix, 'VariableNames', labels, 'RowNames', labels));

% % Reults:
% % TIV= tiv
%                          Age         TIV      Sex(M=1)     Group(struc=1)
%                       _________    _______    _________    ______________
% 
%     Age                       1    0.15592    -0.075446       -0.14668   
%     TIV                 0.15592          1      0.35149        -0.8566   
%     Sex(M=1)          -0.075446    0.35149            1       -0.30107   
%     Group(struc=1)     -0.14668    -0.8566     -0.30107              1 
% %tiv= left nuclei sum volume
%                          Age         TIV       Sex(M=1)     Group(struc=1)
%                       _________    ________    _________    ______________
% 
%     Age                       1    -0.18567    -0.075446       -0.14668   
%     TIV                -0.18567           1      0.14491       -0.10532   
%     Sex(M=1)          -0.075446     0.14491            1       -0.30107   
%     Group(struc=1)     -0.14668    -0.10532     -0.30107              1   

%% Option: plot only struc connectivity (rory's data)
% Color scheme
cm = [1    1    1; 0 0.4470 0.7410];
cm = interp1(cm, 1:0.01:size(cm,1));

% output figure folder
savefolder='structural_RP/8nuclei_unadj/';

% output figure names
Names = selected_nodes_names; % left nuclei

% Develop mean values across the controls
clear map maps
for b = 1:length(thalamus_generic_labels)
    clear map
    for a = 1:length(controls)
        map(a,:) = controls(a).connectome(selected_nodes(b),:);
    end
    maps(b,:) = mean(map,1);
end
maps_c_R = maps; % right nuclei

clear map maps
for b = 1:length(thalamus_generic_labels)
    clear map
    for a = 1:length(controls)
        map(a,:) = controls(a).connectome(selected_nodes(b+8),:);
    end
    maps(b,:) = mean(map,1);
end
maps_c_L = maps; % left nuclei

% Use SimpleBrainPlot to visualise: https://github.com/dutchconnectomelab/Simple-Brain-Plot
for b = 1:length(thalamus_generic_labels)

    %plot R
    mapplot_R = maps_c_R(b,:)';

    lims_R = prctile(mapplot_R,[10 90],"all");
 
    % remove
    mapplot_R(mapplot_R<0)=0;

    plotBrain(labels, mapplot_R, cm,'atlas', 'lausanne120_aseg', 'savePath', [savefolder Names{b}], 'limits',lims_R);

    %plot L
    mapplot_L = maps_c_L(b,:)';

    lims_L = prctile(mapplot_L,[10 90],"all");
 
    % remove
    mapplot_L(mapplot_L<0)=0;

    plotBrain(labels, mapplot_L, cm,'atlas', 'lausanne120_aseg', 'savePath', [savefolder Names{b+8}], 'limits',lims_L);
end
%%12/08/20 Feng: harmonise ind corr matrix before mst
%%090724 Feng:Summary of the Script:s Workflow: 1)Original Weighted Matrix; 2)Minimal Spanning Tree (MST); 3)Remaining Links
% 4)Thresholding; 5)Conversion to Connection-Length Matrix; 6)Graph Theory Measures
% MST is constructed to include the strongest (most positive) correlations, 
% and the rest of the links are handled separately to analyze their incremental addition to the network.
%%25/07/24 Feng: add option--build graph solely by negative and positive corr

addpath 'C:\Users\skgtxfe\OneDrive - University College London\Downloads\2019_03_03_BCT' % add brainconnectivitytoolbox to matlab path
clear
% Load subjects
%infotable = '/Users/fxy/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Data_Code/ICHdatabase/Focal_Con_270125.xlsx';
infotable = 'C:\Users\skgtxfe\OneDrive - University College London\Data_Code\ICHdatabase\Focal_Con_221025.xlsx';
%infotable = 'C:\Users\skgtxfe\OneDrive - University College London\Data_Code\ICHdatabase\TLE_Con_working_090824.xlsx';
datapre = readtable(infotable);% close the excel if unreadable!
%data=datapre(datapre.Include_fmri==1,:);
%data=datapre(datapre.Include_fmri~=0,:);
data=datapre(datapre.Include_vol==1,:);% use include_vol for fmri method v2
subjects = data.conn_id;%subjects = importdata('\sbj.txt'); % List of your subjects

% Load atlas
cd 'C:\Users\skgtxfe\OneDrive - University College London\Study2MultimodelThal\Focal_Con\hubness\pre-processed_data_TS\filtered_with_thalamic_subregion\Harv-ox_thalnuclei_native\' %'/Users/fxy/OneDrive - University College London/data/Piper2021/pre-processed_data_TS/filtered_with_thalamic_subregion/'
ts_path='C:\Users\skgtxfe\OneDrive - University College London\Study2MultimodelThal\Focal_Con\hubness\pre-processed_data_TS\filtered_with_thalamic_subregion\Harv-ox_thalnuclei_native\ind_ts\';
%cd '/Users/fxy/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Study2MultimodelThal/Focal_Con/hubness/pre-processed_data_TS/filtered_with_thalamic_subregion/Harv-ox_thal/'
%ts_path='/Users/fxy/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Study2MultimodelThal/Focal_Con/hubness/pre-processed_data_TS/filtered_with_thalamic_subregion/Harv-ox_thal/ind_ts/';

atlas = load('brainregionFeng.mat').brainregionFeng; % load atlas where ROI are selected
load brainregionFeng.mat % load ROIs
N = length(brainregionFeng); % Number of ROIs
idx_sub = find(ismember(atlas,brainregionFeng));

% Set network density
den_switch = 2;
if den_switch == 1 
% Generate number of links for each density From 89/4005 = 2.22% to 3% 4% ... 100%
% For N=90, the total number of possible edges is N*(N-1)=90×89=4005; minimum num of edges linking all nodes=89 edges,
% To contrsucture a network including all nodes, density should begin > 89/4005(=2.22%), so 3% are used here
    en = zeros(1,99);
    num_steps=99;
    for ii=3:100 %full density: minimum(2.22%)~100%
        en(ii-1)=round(N*(N-1)*ii/(100*2))-(N-1);
    end
    disp('densities 1 to 99')
elseif den_switch == 2
    start_density = 5;
    end_density = 50;
    step_density = 1;
    num_steps = (end_density - start_density) / step_density + 1;
    en = zeros(1, num_steps);
    for ii = start_density:step_density:end_density
        en(ii - start_density + 1) = round(N * (N - 1) * ii / (100 * 2)) - (N - 1);
        if en(ii - start_density + 1) < 0
            en(ii - start_density + 1) = 0; % Ensure no negative values
        end
    end    
    disp('densities 5 to 50')
end

% Preallocating Measurements
% Weighted
W_ot = zeros(N,N,length(subjects),num_steps);
W_ot_bin = zeros(N,N,length(subjects),num_steps);
corr = zeros(N,N,length(subjects));
Egw = zeros(length(subjects),num_steps);
Elw = zeros(N,length(subjects),num_steps);
Mw = zeros(N,length(subjects),num_steps);
Qw = zeros(length(subjects),num_steps);
Zw = zeros(N,length(subjects),num_steps);
Pw = zeros(N,length(subjects),num_steps);
Lw = zeros(length(subjects),num_steps);
CCw = zeros(N,length(subjects),num_steps);
str = zeros(N,length(subjects),num_steps);
BCw = zeros(N,length(subjects),num_steps);
vw = zeros(N,length(subjects),num_steps);
DCw = zeros(N,length(subjects),num_steps);
Modalcontw = zeros(N,length(subjects),num_steps);
% Binary
Egb = zeros(length(subjects), num_steps);
Elb = zeros(N, length(subjects), num_steps);
Lb = zeros(length(subjects), num_steps);
CCb = zeros(N, length(subjects), num_steps);
BCb = zeros(N, length(subjects), num_steps);
vb = zeros(N, length(subjects), num_steps);
DCb = zeros(N,length(subjects),num_steps);
Mb = zeros(N,length(subjects),num_steps);
Qb = zeros(length(subjects),num_steps);

% Correlation matrices for all subjects
all_subjects_corr = [];
for s = 1:length(subjects)
    clear ind_ts AAL AALt r p W Wt SW TW FTW MST Wmst W_b Wm Wmd ind E tmp; 
    load([ts_path 'ind_ts_' subjects{s} '.mat']);disp([ts_path 'ind_ts_' subjects{s} '.mat'])
    AAL = ind_ts(idx_sub,:);
    AALt = AAL';
    [r,~] = corrcoef(AALt);
    r = round(r * 1000000) / 1000000;
    r(isnan(r)) = 0;
    r_vector = r(:);% Reshape the correlation matrix to a one-dimensional vector
    all_subjects_corr = [all_subjects_corr, r_vector];% Concatenate the reshaped vector as a new column in the all_subjects_corr matrix
end

%% Harmonise corr mat before caluclate graph measures
% Load combat var
%Data inputs for ComBat are:
%A data matrix. The data to harmonize. Rows are features (for instance voxels or brain regions) and columns are participants.
%A batch id vector. A vector (length should be equal to the number of columns in the data matrix) that specifies the id for the batch, site, or scanner to correct for. ComBat only accepts one batch vector. You should provide the smallest unit of the study that you believe introduces unwanted variation. For instance, for a study with 2 sites and 3 scanners (1 site with 1 scanner, 1 site with 2 scanners), the id for scanner should be used.

% Biological variables. Optional design matrix specifying biological covariates that should be protected for during the removal of scanner/site effects, such as disease status, age, gender, etc.
age = data.age_fmri_years;
sex = dummyvar(categorical(data.sex));
% TLE vs Control
% TLE=data.TLE;
% mod = [age sex(:,2) TLE];
% TLE vs ETLE vs Control
Focal=dummyvar(categorical(data.Focal));
mod = [age sex(:,2) Focal(:,2:4)];

% scan var
%A vector (length should be equal to the number of columns in the data matrix) that specifies the id for the batch, site, or scanner to correct for. 
% ComBat only accepts one batch vector. 
% You should provide the smallest unit of the study that you believe introduces unwanted variation. 
% For instance, for a study with 2 sites and 3 scanners (1 site with 1 scanner, 1 site with 2 scanners), the id for scanner should be used.

scanner = string(data.Scanner);%use scanner instead of duration because ComBat harmonization does not exhibit the ability to harmonize the test–retest dataset with 20 subjects or less
scanner_protocol = string(data.Scanner_protocol);
%scan_dur=data.scan_dur;

%batch=scanner;%Used in Initial submission to Epilepsia
batch=scanner_protocol;%Used in revision 1 to Epilepsia

% Perform harmonization after removing constant rows
constant_rows = std(all_subjects_corr, 0, 2) == 0;% Identify rows with constant values across all samples
all_subjects_corr_clean = all_subjects_corr(~constant_rows, :);% Remove rows with constant values
all_subjects_corr_clean = combat(all_subjects_corr_clean, batch, mod, 1);
all_subjects_corr_harmo = zeros(size(all_subjects_corr));
all_subjects_corr_harmo(~constant_rows, :) = all_subjects_corr_clean;% Reinsert the harmonized non-constant rows
all_subjects_corr_harmo(constant_rows, :) = all_subjects_corr(constant_rows, :);% Reinsert the constant rows with their original values

% % Regress out age based on control GLM for both controls and patients % % 23/09/24 Feng: Don't regress at these stage!
% % Controls
% clear age idx
% idx{1} = find(data.TLE == 0); name{1} = '1.5TC';nosub{1} = length(idx{1}); age{1} = data.age_fmri_years(idx{1});
% [n, num_subjects] = size(all_subjects_corr_harmo);
% corr_corrected = zeros(size(all_subjects_corr_harmo));
% for j = 1:n
%     j
%     clear control_b
%     % Fit the age GLM using controls data
%     control_b = fitglm(age{1}, all_subjects_corr_harmo(j, idx{1})'); % Fit GLM with age as predictor
% 
%     % Controls data design matrix
%     control_X = [ones(length(age{1}), 1), age{1}]; 
%     % Get residuals from the fitted model
%     control_residuals = all_subjects_corr_harmo(j, idx{1})' - control_X * control_b.Coefficients.Estimate; 
%     % Store residuals (age-adjusted) in corrected matrix for controls
%     corr_corrected(j, idx{1}) = control_residuals';
% 
%     % Patients data
%     patient_idx = find(~ismember(1:num_subjects, idx{1})); % Select patients' indices (those not in idx{1})
%     age_pat = data.age_fmri_years(patient_idx); % Get patient age data
%     patient_X = [ones(length(age_pat), 1), age_pat]; % Patient design matrix      
%     % Apply control's model coefficients to the patient data
%     patient_residuals = all_subjects_corr_harmo(j, patient_idx)' - patient_X * control_b.Coefficients.Estimate; 
%     % Store residuals (age-adjusted) in corrected matrix for patients
%     corr_corrected(j, patient_idx) = patient_residuals';
% end
% all_subjects_corr_harmo = corr_corrected;

%% NOT harmo
all_subjects_corr_harmo = all_subjects_corr;

%% Generate measurements
for s=1:length(subjects)
	s
	clear ind_ts AAL AALt r p W Wt SW TW FTW MST Wmst W_b Wm Wmd ind E tmp;
	% Reshape the r into a 2D n_rois x n_rois matrix
    r = reshape(all_subjects_corr_harmo(:, s), N, N);

    % Build graph including both posi and negative corr
	W = abs(r);% Abslute value

%     % Build graph solely by positive corr
%     r(r<0) = 0; % set all negative corr=0
%     W = abs(r);
    % % Build graph solely by negative corr
    % r(r>0) = 0; % set all postive corr=0
    % W = abs(r);

	% % Put the diagonal to zeros
	for i=1:size(W,1)
		W(i,i)=0;
        r(i,i)=0;
	end
	corr(:,:,s)=r; % save corr r value(after Fisher Transformation)
	
    %Computes minimal span tree, which includes the edges that connect all nodes with the minimal total weight from Wt.
	Wt=1-W; % Generate an inverse matrix that strong connections (high correlation) are transformed into smaller values, making it suitable for MST algorithms that typically find the minimum
    Wt = (Wt + Wt') / 2;  % Force Wt to be symmetric by averaging with its transpose, if 'Error using graph Adjacency matrix must be symmetric.'   
    SW=sparse(Wt); % Transfer to sparse format   
	%TW=graphminspantree(SW);  

% %     % if graphminspantree doesn't work, use imspantree instead to get TW
    TW_graph = minspantree(graph(SW));% Convert SW to a graph object 
    TW = adjacency(TW_graph);

    FTW=full(TW); % Transfer back to double format
	FTW=FTW+FTW'; % Reconstruct symmetry

	% Generate Minimal Span Tree mask
	MST=FTW; 
	MST(MST>0)=1; %binary matrix generated by MST
	Wmst = W.*MST; %weight matrix generated by MST

	% Generate mask for the rest links
	W_b = W;
	W_b(W_b>0) = 1;
	Wm = W_b - MST; %a binary matrix for the remaining links after removing the MST links
	Wmd = W.*Wm; %a weight matrix for the weights of the remaining links after removing the MST links
    Wmd = triu(Wmd); %ensure symmetry is preserved
    ind=find(Wmd); %Find all links
    E=sortrows([ind Wmd(ind)], -2); %Sort by magnitude

	tmp=1;
	for p=1:num_steps % if num_steps=99, Starting by MST,0.022,then 0.03:0.01:1
		p
        clear Wmdt W_thr WL WL_a Dw;
        Wmdt=Wmd;
        Wmdt(E(en(p)+1:end,1))=0; % Apply threshold
        Wmdt = Wmdt+Wmdt.'; % reconstruct symmetry
		W_thr = Wmdt + Wmst; %the thresholded adjacency matrix at a specific threshold level,including the MST and a varying number of the strongest remaining links based on the current threshold.
		
%       % Measures from weighted graph WITHOUT connection's direction(- or +) 
%       W_ot(:,:,s,tmp) = W_thr; % Original matrix:a 4-dimensional array where the original thresholded networks are stored for each subject (s), and for each thresholding step (tmp).
% 		
%       % % Efficiency
% 		Egw(s,tmp) = efficiency_wei(W_thr); % Global
% 		Elw(:,s,tmp) = efficiency_wei(W_thr,1); % Local	
% 
% 		% % Modularity
% 		[Mw(:,s,tmp),Qw(s,tmp)] = community_louvain(W_thr); % weighted %Qw is the modularity value
% 
% 		% % Within-module degree z-score
% 		Zw(:,s,tmp) = module_degree_zscore(W_thr,Mw(:,s,tmp),0); % weighted		
% 
% 		% % Participation coefficient
% 		Pw(:,s,tmp) = participation_coef(W_thr,Mw(:,s,tmp)); % weighted
% 
%         % % convert to connection-length matrix to be used for some graph measures
% 		WL = weight_conversion(W_thr,'lengths'); 
%         % % Characteristic Path Length
% 		Dw = distance_wei(WL); % weighted
% 		Lw(s,tmp) = charpath(Dw); % weighted
% 
%         % % Clustering Coefficient
% 		CCw(:,s,tmp) = clustering_coef_wu(W_thr); % weighted
% 
% 		% % Strength for weighted
% 		str(:,s,tmp) = strengths_und(W_thr);
% 
%         % % Degree centrality
% 		DCw(:,s,tmp) = degrees_und(W_thr); % weighted
% 
%         % % Eigenvector centrality
% 		vw(:,s,tmp) = eigenvector_centrality_und(W_thr); % weighted
% 
%         % % Betweenness Centrality
% 		BCw(:,s,tmp) = betweenness_wei(WL); % Connection weights are ignored in calculations. DCw should equal to DCb.	

        % Measures from weighted graph with connection's direction(- or +) perserved
        W_thr_d = r.*double(W_thr > 0);% Convert thresholded weighted matrix to weighted matrix with connection's direction        
        W_ot(:,:,s,tmp) = W_thr_d; % Original matrix:a 4-dimensional array where the original thresholded networks are stored for each subject (s), and for each thresholding step (tmp).
        
%         % % Modal Controllability
%         NormA = W_thr_d./(1 + svds(W_thr_d, 1));  % Normalizing the adjacency matrix      
%         [U, T] = schur(NormA, 'real');  % Get Schur decomposition      
%         eigVals = diag(T);  % Get eigenvalues from T
%         N = size(NormA, 1);  % Get the size of the network (number of nodes)
%         phi = zeros(N, 1);  % Initialize modal controllability vector            
%         for i = 1:N % Loop through nodes to compute modal controllability for each node
%             phi(i) = (U(i,:).^2) * (1 - eigVals.^2);  % Modal controllability for each node
%         end     
%         Modalcontw(:,s,tmp) = phi;  % Storing modal controllability in a 4D array (optional)
%        
        % % Efficiency
		Egw(s,tmp) = efficiency_wei(W_thr_d); % Global
		%Elw(:,s,tmp) = efficiency_wei(W_thr_d,1); % Local	

% 		% % Modularity and Participation coefficient may cause error when
%         % negative weights inlcuded, thus not to calculate there
% 
%         % % convert to connection-length matrix to be used for some graph measures
% 		WL = weight_conversion(W_thr_d,'lengths'); 
%         % % Characteristic Path Length
% 		Dw = distance_wei(WL); % weighted
% 		Lw(s,tmp) = charpath(Dw); % weighted
% 
%         % % Clustering Coefficient
% 		%CCw(:,s,tmp) = clustering_coef_wu(W_thr_d); 
%  		CCw(:,s,tmp) = clustering_coef_wu(W_thr); % use W_thr(absolute value) instead of W_thr_d(both posi&negative corr) bc including W_thr_d cause errors

		% % Strength for weighted
		str(:,s,tmp) = strengths_und(W_thr_d);

        % % Degree centrality
		DCw(:,s,tmp) = degrees_und(W_thr_d); % weighted
        %DCw_pos(:,s,tmp) = degrees_und(double(W_thr_d > 0));%store DC based on pos corr
        %DCw_neg(:,s,tmp) = degrees_und(double(W_thr_d < 0));%store DC based on neg corr

%         % % Eigenvector centrality
% 		vw(:,s,tmp) = eigenvector_centrality_und(W_thr_d); % weighted
% 
%         % % Betweenness Centrality
% 		BCw(:,s,tmp) = betweenness_wei(WL); 
        
        % Measures from binary graph
        % % Convert thresholded weighted matrix to binary
        W_thr_bin = double(W_thr > 0);
        W_ot_bin(:,:,s,tmp) = W_thr_bin; % Original matrix:a 4-dimensional array where the original thresholded networks are stored for each subject (s), and for each thresholding step (tmp).
        
        % % Efficiency (binary)
        Egb(s, tmp) = efficiency_bin(W_thr_bin); % Global
%         Elb(:, s, tmp) = efficiency_bin(W_thr_bin, 1); % Local
% 
%         % % Modularity (binary)
% 		[Mb(:,s,tmp),Qb(s,tmp)] = community_louvain(W_thr_bin); %Qb is the modularity value
% 
%         % % Characteristic Path Length (binary)
%         Lb(s, tmp) = charpath(distance_bin(W_thr_bin));
% 
%         % % Clustering Coefficient (binary)
%         CCb(:, s, tmp) = clustering_coef_bu(W_thr_bin); 

		% % Degree centrality
		DCb(:,s,tmp) = degrees_und(W_thr_bin); 

%         % % Eigenvector centrality (binary)
%         vb(:, s, tmp) = eigenvector_centrality_und(W_thr_bin); 
% 
%         % % Betweenness Centrality (binary)
%         BCb(:, s, tmp) = betweenness_bin(W_thr_bin); 

		tmp=tmp+1;
	end
end
%save('Harv-oxthal_all_wei_bin_harmo_TLE_ETLE_Con.mat', 'W_ot_bin', 'W_ot', 'corr', 'Egw', 'Elw', 'str', 'Lw', 'CCw', 'DCw', 'BCw', 'vw', 'Egb', 'Elb', 'DCb', 'Lb', 'CCb', 'BCb', 'vb', 'Mb', 'Qb');
%save('Harv-oxthalnuclei_all_wei_bin_harmo_TLE_ETLE_Con.mat', 'W_ot_bin', 'W_ot', 'corr', 'Egw', 'Elw', 'str', 'Lw', 'CCw', 'DCw', 'BCw', 'vw', 'Egb', 'Elb', 'DCb', 'Lb', 'CCb', 'BCb', 'vb', 'Mb', 'Qb');
save('Harv-oxthalnuclei_all_wei_bin_newharmo_scannerprotocol_TLE_ETLE_Con.mat', 'W_ot_bin', 'W_ot', 'corr', 'Egw', 'str', 'DCw', 'Egb', 'DCb');

%save('Harv-oxthalnuclei_all_wei_bin_harmo_posnegDC.mat', 'W_ot_bin', 'W_ot', 'corr', 'DCw', 'DCw_pos', 'DCw_neg');
%save('Harv-oxthalnuclei_all_wei_bin_harmo.mat', 'Modalcontw', 'W_ot_bin', 'W_ot', 'corr', 'Egw', 'Elw', 'str', 'Lw', 'CCw', 'DCw', 'BCw', 'vw', 'Egb', 'Elb', 'DCb', 'Lb', 'CCb', 'BCb', 'vb', 'Mb', 'Qb');
%save('Harv-oxthal_all_wei_bin.mat', 'W_ot', 'Egw', 'Elw', 'str', 'Lw', 'CCw', 'DCw', 'BCw', 'vw', 'Mw', 'Qw', 'Zw', 'Pw', 'Egb', 'Elb', 'DCb', 'Lb', 'CCb', 'BCb', 'vb');
%save('Harv-oxthal8nuclei_all_wei_bin_harmo.mat', 'CCw', '-append');

%% %% Optional: PLOTS THAl
clear
% Load subjects
infotable = '/Users/fxy/OneDrive - University College London/Data_Code/ICHdatabase/TLE_Con_working_090824.xlsx';
%infotable = 'C:\Users\skgtxfe\OneDrive - University College London\Data_Code\ICHdatabase\TLE_Con_working_090824.xlsx';
datapre = readtable(infotable);% close the excel if unreadable!
data=datapre(datapre.Include_fmri==1,:);
subjects = data.conn_id;%subjects = importdata('\sbj.txt'); % List of your subjects

% 1.5T Controls
idx{1} = find(data.TLE == 0); 
name{1} = '1.5TC';

% %1.5T patients
% idx{2} = find((data.TLE == 1 & ~strcmp(data.Scanner, 'Prisma')));
% name{2} = '1.5TP';

% % 3T patients
% idx{2} = find((data.TLE == 1 & strcmp(data.Scanner, 'Prisma')));
% name{2} = '3TP';

% % Overall patients(1.5T+3T)
% idx{2} = find(data.TLE == 1);
% name{2} = '1.5T+3T P';
% idx{2} = find((data.TLE == 1& strcmp(data.mri_side, 'left')));
% name{2} = '1.5T+3T left TLE';
% idx{2} = find((data.TLE == 1& strcmp(data.mri_side, 'right')));
% name{2} = '1.5T+3T right TLE';

% % Mesial TLE with mts vs. mesial TLE with other pathology, 1.5T
% idx{2} = find((data.mts == 1 & ~strcmp(data.Scanner, 'Prisma')));%find(data.mts == 1);
% name{2} = 'mTLE_mts';
% idx{3} = find((data.mesial == 1 & data.nomts == 1 & ~strcmp(data.Scanner, 'Prisma')));%find(data.mesial == 1 & data.nomts == 1);
% name{3} = 'mTLE_other';

% %%Mesial vs. lateral/neoTLE
% idx{2} = find((data.mesial == 1 ));
% name{2} = 'mesial TLE';
% idx{3} = find((data.lateral == 1 ));
% name{3} = 'lateral TLE';
% 
% % Mts vs. no mts
% idx{2} = find(data.mts == 1);% mts/no-sf all_grp{2}
% name{2} = 'mts';
% idx{3} = find(data.nomts == 1); % no-mts/sf all_grp{3}
% name{3} = 'nomts';

% % Lesion affecting HC vs. not affecting HC, 1.5T
% idx{2} = find((data.HC_lesion == 1 & ~strcmp(data.Scanner, 'Prisma')));
% name{2} = 'HC_lesion';
% idx{3} = find((data.noHC_lesion == 1 & ~strcmp(data.Scanner, 'Prisma')));
% name{3} = 'noHC_lesion';

% %  No sf vs. sf
% idx{2} = find((data.sf_latest == 0 ));
% name{2} = 'nosf';
% idx{3} = find((data.sf_latest == 1 ));
% name{3} = 'sf';
%  No sf vs. sf, 1.5T
% idx{2} = find((data.sf_latest == 0 & ~strcmp(data.Scanner, 'Prisma')));
% name{2} = 'nosf';
% idx{3} = find((data.sf_latest == 1 & ~strcmp(data.Scanner, 'Prisma')));
% name{3} = 'sf';
% % Overall patient group combing NS and SF
% idx{2} = find(((data.sf_latest == 0 | data.sf_latest == 1) & ~strcmp(data.Scanner,'Prisma')));
% name{2} = 'NSF+SF';

% %%No sf vs. sf, 1.5T, ATLR
% idx{2} = find((data.sf_latest == 0 & data.ATLR==1 & ~strcmp(data.Scanner, 'Prisma')));
% name{2} = 'nosf-ATLR';
% idx{3} = find((data.sf_latest == 1 & data.ATLR==1 & ~strcmp(data.Scanner, 'Prisma')));
% name{3} = 'sf-ATLR';

% %fbtcs vs. not-fbtcs
idx{2} = find((strcmp(data.types_focal_to_bilateral_tonic_clonic,'yes')));
name{2} = 'FBCTS';
idx{3} = find((strcmp(data.types_focal_to_bilateral_tonic_clonic,'no')));
name{3} = 'Not-FBCTS';
%fbtcs vs. not-fbtcs,1.5T
% idx{2} = find((strcmp(data.types_focal_to_bilateral_tonic_clonic,'yes')& ~strcmp(data.Scanner, 'Prisma')));
% name{2} = 'FBCTS';
% idx{3} = find((strcmp(data.types_focal_to_bilateral_tonic_clonic,'no')& ~strcmp(data.Scanner, 'Prisma')));
% name{3} = 'Not-FBCTS';

% %motor vs. not-motor
% idx{2} = find((strcmp(data.types_motor,'all_motor') | strcmp(data.types_motor,'mixed')));
% name{2} = 'Motor';
% idx{3} = find((strcmp(data.types_motor,'all_nonmotor')));
% name{3} = 'Nonmotor';

nosub{1} = length(idx{1}); age{1} = data.age_fmri_years(idx{1});epilepsy_duration{1} = data.duration_epilepsy_fmri_years(idx{1});sex{1} = data.sex(idx{1});
nosub{2} = length(idx{2}); age{2} = data.age_fmri_years(idx{2});epilepsy_duration{2} = data.duration_epilepsy_fmri_years(idx{2});sex{2} = data.sex(idx{2});
if length(name) > 2
    nosub{3} = length(idx{3}); age{3} = data.age_fmri_years(idx{3});epilepsy_duration{3} = data.duration_epilepsy_fmri_years(idx{3});sex{3} = data.sex(idx{3});
end
disp([name nosub]);

%% Load graph metrics
clear all_grp*
%cd 'C:\Users\skgtxfe\OneDrive - University College London\Study2-4ThalnucleiFC\Focal_Con\hubness\pre-processed_data_TS\filtered_with_thalamic_subregion\Harv-ox_thal\';
cd '/Users/fxy/OneDrive - University College London/Study2-4ThalnucleiFC/Focal_Con/hubness/pre-processed_data_TS/filtered_with_thalamic_subregion/Harv-ox_thal8nuclei/';
% load Harv-oxthalnuclei_all_wei_bin_harmo_nomst.mat
% nw='Harv-oxthal_nomst';
load Harv-oxthalnuclei_all_wei_bin_harmo_posnegDC.mat
nw='Harv-oxthal_harmo_posneg';

load brainregionFeng.mat
load brainregionFengidx.mat
roi=brainregionFeng;

den_switch=2;
den_thr = 5:50; % num_step=46
%den_thr = 1:99; % num_step=99
%den_thr = 1:100; % num_step=100

for g = 1:numel(name)
    for ii = 1:length(den_thr)
        % Binary graph   
        graph='bin';
%         all_grp_DC{g,ii}(:,:) = DCb(:,idx{g},ii)'; % Degree centrality
        all_grp_DC_neg{g,ii}(:,:) = DCw_neg(:,idx{g},ii)'; 
        all_grp_DC_pos{g,ii}(:,:) = DCw_pos(:,idx{g},ii)'; 
%         all_grp_EC{g,ii}(:,:) = vb(:,idx{g},ii)'; % Eigenvector centrality
%         all_grp_BC{g,ii}(:,:) = BCb(:,idx{g},ii)'; % Betweenness Centrality
%         all_grp_CC{g,ii}(:,:) = CCb(:,idx{g},ii)'; % Clustering Coefficient

%         % Weight graph   
%         graph='wei';
%         all_grp_NS{g,ii}(:,:) = str(:,idx{g},ii)'; % Node strength
%         all_grp_EC{g,ii}(:,:) = vw(:,idx{g},ii)'; % Eigenvector centrality
%         all_grp_BC{g,ii}(:,:) = BCw(:,idx{g},ii)'; % Betweenness Centrality
%         all_grp_CC{g,ii}(:,:) = CCw(:,idx{g},ii)'; % Clustering Coefficient
    end
end

%% at which density, connectivity mat starts to have negative values 
% Loop over each group
for g = 1:numel(name)    
    negative_density{g} = [];% Initialize a list to store densities where negative values are found
    
    for ii = 1:length(den_thr)        
        % Get the matrix for the current group and density threshold
        mat = squeeze(mean(W_ot(107:108,:,idx{g},ii), 3));
        W_ot_grp_mean{g,ii}(:,:)=mat;
               
        if any(mat(:) < 0) % Check if there are any negative values in the matrix
            % Store the density threshold index (or actual threshold) where negative values are found
            negative_density{g} = [negative_density{g}, ii];
        end
    end
end
for g = 1:numel(name)
    if ~isempty(negative_density{g})
        fprintf('Group %s has negative values at density thresholds: %s\n', name{g}, num2str(negative_density{g}));
    else
        fprintf('Group %s has no negative values at any density threshold.\n', name{g});
    end
end

%% Start Plot
% Load ROIs % NOTE:use correct brainregionidx of the atlas!
clear hub_roi titles
% AAL_thal
%sig_node=[75 76]; % overall TLE vs. controls, 1.5T
%sig_node=[37 38]; % NSF vs. SF sig node: Hippocampus_R 38
% Harv-ox 
thalnuclei=length(roi) - 7:length(roi);
ant=[107 108];

hub_roi=ant;
titles = roi(hub_roi,:);

% Select graph measures to plot
clear GM all_grp_GM
% Bin
%all_grp_GM=all_grp_DC;GM='Degree centrality';
%all_grp_GM=all_grp_DC_pos;GM='Degree centrality no. postive correlation';
all_grp_GM=all_grp_DC_neg;GM='Degree centrality no. anticorrelation';
%all_grp_GM=all_grp_EC;GM='Eigenvector centrality';
%all_grp_GM=all_grp_BC;GM='Betweenness centrality';
%all_grp_GM=all_grp_CC;GM='Clustering coefficient';
% Weight
%  all_grp_GM=all_grp_NS;GM='Node strength';
% all_grp_GM=all_grp_EC;GM='Eigenvector centrality';
% all_grp_GM=all_grp_BC;GM='Betweenness centrality';
% all_grp_GM=all_grp_CC;GM='Clustering coefficient';

addpath 'C:\Users\skgtxfe\OneDrive - University College London\Study2-4ThalnucleiFC\Focal_Con\hubness\Feng_code\Piper2021\Rory_matlab_code\';
close all
for kk=1:length(hub_roi)
    disp(roi(hub_roi(kk)))
    options.handle = figure(kk);
    options.error = 'sem';% if 'std',one standard deviation; if 'sem', standard error mean; if 'var',one variance;  if 'c95',95% confidence interval. 

%     % Controls
%     for d = 1:length(den_thr)
%         for q = 1:nosub{1}  
%             if mod(kk, 2) == 1  % Check if kk is odd
%                 tmp_1 = squeeze(all_grp_GM{1,d}(q, hub_roi(kk)));
%                 tmp_2 = squeeze(all_grp_GM{1,d}(q, hub_roi(kk + 1)));
%             else  % kk is even
%                 tmp_1 = squeeze(all_grp_GM{1,d}(q, hub_roi(kk - 1)));
%                 tmp_2 = squeeze(all_grp_GM{1,d}(q, hub_roi(kk)));
%             end
%             tmp_gm = mean([tmp_2 tmp_1]);
%             gm_controls(q,d) = tmp_gm;
%         end
%     end
%     clear tmp1 tmp2
%     plot_areaerrorbar_rory(gm_controls, options,[100 100 255]./300) % blue
%     hold on
     
    % Patients 
    for d=1:length(den_thr)
        for q=1:nosub{2}
            gm_patients(q,d) = all_grp_GM{2,d}(q,hub_roi(kk));
        end
    end
    plot_areaerrorbar_rory(gm_patients, options,[255, 70, 70]./300) % red
    hold on
    
    % G3
    for d=1:length(den_thr)
        for q=1:nosub{3}
            gm_group3(q,d) = all_grp_GM{3,d}(q,hub_roi(kk));
        end
    end
    plot_areaerrorbar_rory(gm_group3, options,[50, 240, 100]./300)% green
    hold on

    %%Adjust asix for neg plot
   set(gca, 'YDir', 'reverse'); % Invert Y axis
   set(gca, 'XAxisLocation', 'top'); % Place x-axis at the top
    %%Move x-axis labels, ticks, and name to the inside of the plot
    set(gca, 'TickDir', 'none'); % Moves tick marks to the inside

   xlabel('% Network Density')
     
     if den_switch ==1
     ylim([0 inf])
     title(titles(kk))
     
     elseif den_switch ==2
    labells = [5, 20, 35, 50];% Define custom labels and their respective positions
    positions = [1, 15, 30, 45]; % Positions to apply the labels on the x-axis
    xticks(positions); 
    xticklabels(num2str(labells')); 
    title(titles(kk));
    xlim([1, 45]);
    ylim([0, inf]); 
     end

     ylabel(GM)
     %legend('Controls SEM','Controls median','Patients SEM','Patients median','Location','southeast')
     %legend('controls SEM','controls median','mts SEM','mts median','no-mts SEM','no-mts median','Location','southeast')
     %legend('Controls SEM','Controls median','NSF SEM','NSF median','SF SEM','SF median','Location','southeast')
     %legend('Controls SEM','Controls median','FBTCS SEM','FBTCS median','NFBTCS SEM','NFBTCS median','Location','southeast')
     %legend('FBTCS SEM','FBTCS median','NFBTCS SEM','NFBTCS median','Location','southeast')
     %legend('Motor SEM','Motor median','Nonmotor SEM','Nonmotor median','Location','southeast')
   
    % Set transparent background, copy to powerpoint to save transparent
     set(gcf, 'Color', 'none'); set(gca, 'Color', 'none');
  
    hold off
    
%     %save results
%     gm_c(:,:,kk) = gm_controls;
%     gm_p(:,:,kk) = gm_patients;
%     %gm_g3(:,:,kk) = gm_group3;
    
    save_roi = convertStringsToChars(titles(kk));
    %saveas(gcf,[GM '_' graph '_' save_roi '_' nw '_' name{2} '.png'])
    %close all    
end

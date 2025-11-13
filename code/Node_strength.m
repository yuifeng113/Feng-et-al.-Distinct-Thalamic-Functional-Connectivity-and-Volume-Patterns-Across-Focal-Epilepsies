
addpath '~\Downloads\2019_03_03_BCT' % add brainconnectivitytoolbox to matlab path

clear

% Load variables
data = readtable('~/subj_data.xlsx');
subjects = data.conn_id

% Load connectivity matrices
cd '~\connectivity-matrices\' 
ts_path='~\connectivity-matrices\';

atlas = load('brainregionFeng.mat').brainregionFeng; % load ROI label
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

%% Harmonise connectivity matrices before compute node strength
%Data inputs for ComBat are:
%A data matrix. The data to harmonize. Rows are features (for instance voxels or brain regions) and columns are participants.
%A batch id vector. A vector (length should be equal to the number of columns in the data matrix) that specifies the id for the batch, site, or scanner to correct for. ComBat only accepts one batch vector. You should provide the smallest unit of the study that you believe introduces unwanted variation. For instance, for a study with 2 sites and 3 scanners (1 site with 1 scanner, 1 site with 2 scanners), the id for scanner should be used.

% Biological variables. 
% Optional design matrix specifying biological covariates that should be protected for during the removal of scanner/site effects, such as disease status, age, gender, etc.
age = data.age_fmri_years;
sex = dummyvar(categorical(data.sex));
Focal=dummyvar(categorical(data.Focal));
mod = [age sex(:,2) Focal(:,2:4)];

% Scan variables
% A vector (length should be equal to the number of columns in the data matrix) that specifies the id for the batch, site, or scanner to correct for. 
% ComBat only accepts one batch vector. 
% You should provide the smallest unit of the study that you believe introduces unwanted variation. 
% For instance, for a study with 2 sites and 3 scanners (1 site with 1 scanner, 1 site with 2 scanners), the id for scanner should be used.

scanner_protocol = string(data.Scanner_protocol);
batch=scanner_protocol;

% Perform harmonization after removing constant rows
constant_rows = std(all_subjects_corr, 0, 2) == 0;% Identify rows with constant values across all samples
all_subjects_corr_clean = all_subjects_corr(~constant_rows, :);% Remove rows with constant values
all_subjects_corr_clean = combat(all_subjects_corr_clean, batch, mod, 1);
all_subjects_corr_harmo = zeros(size(all_subjects_corr));
all_subjects_corr_harmo(~constant_rows, :) = all_subjects_corr_clean;% Reinsert the harmonized non-constant rows
all_subjects_corr_harmo(constant_rows, :) = all_subjects_corr(constant_rows, :);% Reinsert the constant rows with their original values

%% NOT Combat harmonisation
all_subjects_corr_harmo = all_subjects_corr;

%% Generate node strength (and other graph metrics)
for s=1:length(subjects)
	s
	clear ind_ts AAL AALt r p W Wt SW TW FTW MST Wmst W_b Wm Wmd ind E tmp;
	% Reshape the r into a 2D n_rois x n_rois matrix
    r = reshape(all_subjects_corr_harmo(:, s), N, N);

    % Build graph including both posi and negative corr
	W = abs(r);% Abslute value

	% Put the diagonal to zeros
	for i=1:size(W,1)
		W(i,i)=0;
        r(i,i)=0;
	end
	corr(:,:,s)=r; % save corr r value(after Fisher Transformation)
	
    %Computes minimal span tree, which includes the edges that connect all nodes with the minimal total weight from Wt.
	Wt=1-W; % Generate an inverse matrix that strong connections (high correlation) are transformed into smaller values, making it suitable for MST algorithms that typically find the minimum
    Wt = (Wt + Wt') / 2;  % Force Wt to be symmetric by averaging with its transpose, if 'Error using graph Adjacency matrix must be symmetric.'   
    SW=sparse(Wt); % Transfer to sparse format   

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
		
        % Measures from weighted graph with connection's direction(- or +) perserved

        W_thr_d = r.*double(W_thr > 0);% Convert thresholded weighted matrix to weighted matrix with connection's direction        
        W_ot(:,:,s,tmp) = W_thr_d; % Original matrix:a 4-dimensional array where the original thresholded networks are stored for each subject (s), and for each thresholding step (tmp).
        
        % Modal Controllability
        NormA = W_thr_d./(1 + svds(W_thr_d, 1));  % Normalizing the adjacency matrix      
        [U, T] = schur(NormA, 'real');  % Get Schur decomposition      
        eigVals = diag(T);  % Get eigenvalues from T
        N = size(NormA, 1);  % Get the size of the network (number of nodes)
        phi = zeros(N, 1);  % Initialize modal controllability vector            
        for i = 1:N % Loop through nodes to compute modal controllability for each node
            phi(i) = (U(i,:).^2) * (1 - eigVals.^2);  % Modal controllability for each node
        end     
        Modalcontw(:,s,tmp) = phi;  % Storing modal controllability in a 4D array (optional)

        % Efficiency
		Egw(s,tmp) = efficiency_wei(W_thr_d); % Global
		Elw(:,s,tmp) = efficiency_wei(W_thr_d,1); % Local	

        % Convert to connection-length matrix to be used for some graph measures
		WL = weight_conversion(W_thr_d,'lengths'); 

        % Characteristic Path Length
		Dw = distance_wei(WL); % weighted
		Lw(s,tmp) = charpath(Dw); % weighted

        % Clustering Coefficient
		%CCw(:,s,tmp) = clustering_coef_wu(W_thr_d); 
 		CCw(:,s,tmp) = clustering_coef_wu(W_thr); % use W_thr(absolute value) instead of W_thr_d(both posi&negative corr) bc including W_thr_d cause errors

		% Node Strength 
		str(:,s,tmp) = strengths_und(W_thr_d);

        % Degree centrality
		DCw(:,s,tmp) = degrees_und(W_thr_d); % weighted

        % Eigenvector centrality
		vw(:,s,tmp) = eigenvector_centrality_und(W_thr_d); % weighted

        % Betweenness Centrality
		BCw(:,s,tmp) = betweenness_wei(WL); 

		tmp=tmp+1;
	end
end
save('savename.mat', 'W_ot_bin', 'W_ot', 'corr', 'Egw', 'str', 'DCw', 'Egb', 'DCb');

% gac.m
% Growth and Aggregation of Clusters model, attempting to implement via
% Gillespie's algorithm.  not sure if it will produce gains though, each
% clump is effectively like another species.  doens't scale very well.
%

function [V_arr,total_pop_arr,tvec] = gac_gillespie(lr,la,lc,ls,Tmax,n0,nu_A,nu_F)

%% default values for input parameters

% rate of cluster growth
if ~exist('lr','var')||isempty(lr)
    lr = 1;
end

% rate of cluster aggregation
if ~exist('la','var')||isempty(la)
    la = .1;
end

% rate of cluster collapse
if ~exist('lc','var')||isempty(lc)
    lc = .005;
end

% rate of cluster sprouting
if ~exist('ls','var')||isempty(ls)
    ls = .1;
end

% total simulation time
if ~exist('Tmax','var')||isempty(Tmax)
    Tmax = 24;
end

% initial starting population of single cells
if ~exist('n0','var')||isempty(n0)
    n0 = 10;
end

% aggregation exponent
if ~exist('nu_A','var')||isempty(nu_A)
    nu_A = 2/3;
end

% fragmentation exponent
if ~exist('nu_F','var')||isempty(nu_F)
    nu_F = 2/3;
end

% timestep for deterministic growth
% if ~exist('dt','var')||isempty(dt)
%     dt = .00001;
% end

% shuffle random number generator
rng('shuffle');


%% intialize some arrays

% main object of the simulation:  array of cluster volumes
if numel(n0)==1
    V_arr = ones(1,n0);
elseif size(n0,1) > 1
    disp('gac error:  initial cluster size array n0 must have dim 1xM');
    return
else
    V_arr = n0;
end


% array for keeping track of total population size
total_pop_arr = sum(n0);

% pregenerate lots of random numbers for speed.
number_of_random_numbers_to_pre_generate = round(1000);
lots_of_random_numbers = rand(1,number_of_random_numbers_to_pre_generate);

% counter for how many of the pre-generated random numbers have been used.
n = 0;

% time
t = 0;

tvec = t;

% hour marker
disp_time_marker = 1;
disp_time_increment = 1;

% Aggregation matrix
%A_mat = zeros(numel(V_arr));

%% fixed parameters
% carrying capacity of whole system.  regulates cluster growth in 'global' mode.  
max_total_pop = 10^(4);                       

% logical for printing progress to screen
l_print_progress = true;


%% main simulation.  loop over time.
while t < Tmax 
        
    
    % print update to console
    if l_print_progress && t >= disp_time_marker
        disp(['time = ' num2str(t) ' number of clusters = ' num2str(numel(V_arr))])
        disp_time_marker = disp_time_marker + disp_time_increment;
    end   
    
      
    % compute probability of collapse as a function of cluster volume.
    % here, the relationship is linear in linear length dimension.
    lambda_C = lc.*(V_arr).^(1/3);
    
    ids_for_collapse = 1:numel(lambda_C);
    
    label_for_lambda_C = 1.*ones(1,numel(lambda_C));
    

    % compute the probability of aggregation in this timestep.
    % made product kernel 5/1/18
    %lambda_A = la.*dt.*V_arr.^nu_A;
    A_mat = updateA_mat();
    lambda_A = A_mat(A_mat>0); 
    
    if size(lambda_A,1) > 1
        lambda_A = lambda_A';
    end
    
    label_for_lambda_A = 2.*ones(1,numel(lambda_A));
    
    %[rows_of_nonzero_Amat,cols_of_nonzero_Amat] = find(A_mat>0);
    lin_inds_for_Amat = find(A_mat>0);
    
    ids_for_agg = 1:numel(lambda_A);
    
    
    % compute probability of sprouting.  require clusters to be of at least
    % size 2 to sprout.
    lambda_F = ls.*(V_arr >= 2).*(V_arr).^(nu_F);
    
    ids_for_frag = 1:numel(lambda_F);
    
    ids_arr = [ids_for_collapse, ids_for_agg, ids_for_frag];
    
    label_for_lambda_F = 3.*ones(1,numel(lambda_F));
        
    lambda_arr = [lambda_C, lambda_A, lambda_F];
    
    lambda_total = sum(lambda_arr);
    
    all_labels = [label_for_lambda_C, label_for_lambda_A, label_for_lambda_F]; 
    
    % choose reaction time
    delta_t = exprnd(1./lambda_total);
    
    % choose next reaction
   % collect a random number.  If we've run out of pregenerated ones,
    % generate some new ones.
    if n + 1 < number_of_random_numbers_to_pre_generate
        random_number = lots_of_random_numbers(n+1);
        n = n + 1;
    else
        random_number = rand();
    end
    
    reaction_id = find(random_number < cumsum(lambda_arr)./lambda_total,1);
    
    if ~isempty(reaction_id)
        [V_arr,tvec,total_pop_arr] = gac_gillespie_update();
    end
    
    t = t + delta_t;
    
    %tvec = [tvec t];
    
    
    
    %total_pop_arr = [total_pop_arr, sum(V_arr)];
    
    if isempty(V_arr)
        return
    end
    
end
    
    % gillespie update function.  nested to inherit random numbers
    function [V_arr_out,tvec_out,tot_pop_out] = gac_gillespie_update()
        
        this_label = all_labels(reaction_id);
        
        V_arr_out = V_arr;
        
        switch this_label
            case 1
                % collapse

                % remove clusters that have collapsed from the cluster array.
                %V_arr = V_arr(random_numbers >= lambda_arr(reaction_id));
                V_arr_out(reaction_id) = [];
                
            case 2
                % aggregation
                
                
                % find ids of clusters that will aggregate in this timestep
                agg_id = ids_arr(reaction_id);
                
                lin_ind_of_this_agg = lin_inds_for_Amat(agg_id);
                
                [agg_row,agg_col] = ind2sub(size(A_mat),lin_ind_of_this_agg);
                
                V_arr_out(agg_row) = V_arr_out(agg_row) + V_arr_out(agg_col);
                
                V_arr_out(agg_col) = [];
                
                
                
            case 3
                %% sprout//fragment
                
                frag_id = ids_arr(reaction_id);
                
                % reduce sprouted clusters by 1
                V_arr_out(frag_id) = V_arr_out(frag_id) - 1;
                
                % add this collection of new single cells to cluster array
                V_arr_out = [V_arr_out, 1];
                
        end
                
                
        %% growth (deterministic)
        
        dt = delta_t./10;
        tvec_out = [tvec tvec(end):dt:(t+delta_t)];
        
        numsteps = numel(t:dt:(t+delta_t));
        
        tmp_tot_pop = zeros(1,numsteps);
        
        for s = 1:numsteps
        
            % grow all clusters in one deterministic vectorized step
            V_arr_out = V_arr_out + dt.*lr.*V_arr_out.*(1-sum(V_arr_out)./max_total_pop);
            
            % if clusters die out, remove them
            V_arr_out = V_arr_out(V_arr_out>=1);
            
            % update tmp_totpop
            tmp_tot_pop(s) = sum(V_arr_out);
        
        end
        
        tot_pop_out = [total_pop_arr, tmp_tot_pop];
                  
    end

    
    function A_mat = updateA_mat()
        
        A_mat = zeros(numel(V_arr));
        
        %IDs = 1:numel(V_arr);
        
        %for i = 1:numel(IDs)
        for i = 1:(numel(V_arr)-1)
           
            %culled_IDs = IDs(IDs~=i);
            
            %for j = 1:numel(culled_IDs)
            for j = (i+1):numel(V_arr)
                
                %A_mat(j,i) = la.*(V_arr(i).*V_arr(IDs(j))).^nu_A;
                A_mat(j,i) = la.*(V_arr(i).*V_arr(j)).^nu_A;
                
            end
            
        end
        
     
    end
    
    
end










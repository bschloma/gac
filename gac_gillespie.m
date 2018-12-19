% Program:  gac_gillespie.m
% 
% Summary:  Growth and Aggregation of Clusters model, implemented with
%           Gillespie's direct method.  The main object is an array of
%           cluster sizes.  Clusters are subject to 4 rate processes,
%           aggregation, fragmentation, growth, and dispersal, each with a
%           rate parameter that can depend on cluster size.  Since growth 
%           is the fastest of the rates, and we don't care about
%           demographic stochasticity, it is approximated as a
%           deterministic process.  This means that clusters sizes are no
%           longer integers.
%
%           The simulation starts from a specified initial condition and
%           runs for a specified amount of time.  A subset of the
%           simulation results are output as arrays and an option exists to
%           write every cluster size at every time point to a txt file.
%
% Inputs:   lr - (float) growth rate
%           la - (float) aggregation rate
%           lc - (float) collapse/dispersal rate
%           ls - (float) fragmentation/sprout rate
%           Tmax - (float) simulation time to exit
%           n0 - (int or 1x[number of clusters] array of floats/ints) specifies initial
%               conditions. if numel(n0)=1, simulation starts with round(n0)
%               monomers.  if numel(n0) > 1, n0 is taken to be the initial
%               cluster size array.
%           nu_A - (float) exponent determining how the aggregation rate
%                   scales with cluster size via (agg rate) = la*(cluster
%                   size)^nu_A.
%           nu_F - (float) exponent determining how the fragmentation rate
%                   scales with cluster size via (frag rate) = la*(cluster
%                   size)^nu_F.
%           lwritetxt - (logical) 0 for not writing to txt file, 1 for yes.
%           txtdir - (str) full path to dir to save txt files
%           txtname - (str) name of txt file
%
% Outputs:  V_arr - (1x(number of clusters) array of floats) sizes of all
%                   clusters at the final time point
%           total_pop_arr - (1x(number of time steps) array of floats)
%                   total population (sum of all cluster sizes) over time
%           tvec - (1x(number of time steps) array of floats) array of
%                   times at which reactions occured
%           num_clumps_arr - (1x(number of time steps)) number of clusters
%                   over time
%
% Author:   Brandon Schlomann
%
% Date:     Summer 2018 - first written
%
% VCS:      github.com/bschloma/gac
%

function [V_arr,total_pop_arr,tvec,num_clumps_arr] = gac_gillespie(lr,la,lc,ls,Tmax,n0,nu_A,nu_F,max_total_pop,lwritetxt,txtdir,txtname)

%% default values for input parameters

% rate of cluster growth
if ~exist('lr','var')||isempty(lr)
    lr = 1;
end

% rate of cluster aggregation
if ~exist('la','var')||isempty(la)
    la = 1;
end

% rate of cluster collapse
if ~exist('lc','var')||isempty(lc)
    lc = .1;%.005;
end

% rate of cluster sprouting
if ~exist('ls','var')||isempty(ls)
    ls = 20;
end

% total simulation time
if ~exist('Tmax','var')||isempty(Tmax)
    Tmax = 72;
end

% initial starting population of single cells
if ~exist('n0','var')||isempty(n0)
    n0 = 10;
end

% aggregation exponent
if ~exist('nu_A','var')||isempty(nu_A)
    nu_A = 0;
end

% fragmentation exponent
if ~exist('nu_F','var')||isempty(nu_F)
    nu_F = 0;
end

% carrying capacity
if ~exist('max_total_pop','var')||isempty(max_total_pop)
    max_total_pop = 1e5;
end

% logical for saving to txt file
if ~exist('lwritetxt','var')||isempty(lwritetxt)
    lwritetxt = false;
end

% path to dir for saving txt file
if ~exist('txtdir','var')||isempty(txtdir)
    txtdir = pwd;
end

% name for saving txt file
if ~exist('txtname','var')||isempty(txtname)
    txtname = 'gacout';
end

% shuffle random number generator
rng('shuffle');


%% intialize some arrays

% main object of the simulation:  array of cluster volumes
if numel(n0)==1
    V_arr = ones(1,round(n0));
    
    % array for keeping track of number of clumps over time
    num_clumps_arr = n0;

elseif size(n0,1) > 1
    disp('gac error:  initial cluster size array n0 must have dim 1xM');
    return
else
    V_arr = n0;
    
     % array for keeping track of number of clumps over time
    num_clumps_arr = numel(n0);
    
end


% array for keeping track of total population size
total_pop_arr = sum(n0);

% pregenerate lots of random numbers for speed.
number_of_random_numbers_to_pre_generate = round(100000);
lots_of_random_numbers = rand(1,number_of_random_numbers_to_pre_generate);

% counter for how many of the pre-generated random numbers have been used.
n = 0;

% time
t = 0;
tvec = t;

% hour marker
disp_time_marker = 1;
disp_time_increment = 1;

% for printing warning about random number usage
l_first_time_running_out_of_randns = 1;

% exponent for collapse
nu_C = 0;

%% fixed parameters
% carrying capacity, sets maximum total population size.  
%max_total_pop = 10^(4);                       

% growth timestep
fraction_of_delta_t = 1;
baseline_dt = .1;

% logical for printing progress to screen
l_print_progress = false;

% for writing to txt file, set things up, print a header and time zero data
if lwritetxt
    
    % if desired save directory doesn't exist, make it
    if ~exist(txtdir,'dir')
        mkdir(txtdir);
    end
    
    % handle to txt file, open in append mode
    fid = fopen([txtdir filesep txtname],'a');
    
    % write header
    fprintf(fid,'%s %s\n','time','cluster sizes');
    
    % time zero data to write
    outarr = [t,V_arr];
    
    % format string for data
    format_str = [repmat('%f ',1,numel(outarr)-1), '%f\n'];
    
    % write time zero data
    fprintf(fid,format_str,outarr);
    
    % close the file
    fclose(fid);
    
end
    
%% main simulation.  loop over time.
while t < Tmax 
            
    % print update to console
    if l_print_progress && t >= disp_time_marker
        disp(['time = ' num2str(t) ' number of clusters = ' num2str(numel(V_arr))])
        disp_time_marker = disp_time_marker + disp_time_increment;
    end   
          
    % compute probability of collapse as a function of cluster volume.
    % here, the relationship is linear in linear length dimension.
    lambda_C = lc.*(V_arr).^(nu_C);
    
    % create an array of ids labelling every possible collapse reaction
    ids_for_collapse = 1:numel(lambda_C);
    
    % create an accesory array that labels (1) these reactions as "collapse"
    label_for_lambda_C = 1.*ones(1,numel(lambda_C));
    
    % compute the probability of aggregation in this timestep.
    % A_mat is a lower triangular array that keeps track of the
    % proababilities for every possible aggregation reaction.
    A_mat = la.*tril((V_arr'*V_arr).^nu_A,-1);
    
    % extract the non-zero rates into a linear array
    lambda_A = A_mat(A_mat>0); 
    
    % double check dimensions
    if size(lambda_A,1) > 1
        lambda_A = lambda_A';
    end
    
    % create labels (2) that keep track of these possible reactions as
    % "aggregation"
    label_for_lambda_A = 2.*ones(1,numel(lambda_A));
    
    % keep track of the indices of the non-zero reactions in A_mat's
    % coordinates.  this will be used to identify the clusters in V_arr
    % that are participating in the reaction.
    lin_inds_for_Amat = find(A_mat>0);
    
    % assign ids to each of the possible aggregation reactions.
    ids_for_agg = 1:numel(lambda_A);
        
    % compute probability of sprouting.  require clusters to be of at least
    % size 2 to sprout.
    lambda_F = ls.*(V_arr >= 2).*(V_arr).^(nu_F);
    
    % assign ids to each of the possible fragmentation reactions.
    ids_for_frag = 1:numel(lambda_F);
    
    % assemble all of the reaction ids into a single linear array
    ids_arr = [ids_for_collapse, ids_for_agg, ids_for_frag];
    
    % create labels (3) that keep track of these possible reactions as
    % "fragmentation" 
    label_for_lambda_F = 3.*ones(1,numel(lambda_F));
    
    % assemble all reaction probability rates into a single linear array
    lambda_arr = [lambda_C, lambda_A, lambda_F];
    
    % compute total proabability rate of a reaction happening
    lambda_total = sum(lambda_arr);
    
    % assemble all labels into a single linear array
    all_labels = [label_for_lambda_C, label_for_lambda_A, label_for_lambda_F]; 
    
    % choose reaction time based on the total probability rate of a
    % reaction happening.
    
    % collect a random number.  If we've run out of pregenerated ones,
    % generate some new ones.
    if n + 1 < number_of_random_numbers_to_pre_generate
        random_number = lots_of_random_numbers(n+1);
        n = n + 1;
    else
        disp('gac: ran out of pre-generated rns, generating new ones')

        % generate new ones
        lots_of_random_numbers = rand(1,number_of_random_numbers_to_pre_generate);
        
        % collect the one we need
        random_number = lots_of_random_numbers(1);
        
        % reset counter
        n = 1;
        
    end
    
    delta_t = (1/lambda_total).*log(1/random_number);
    
    % choose next reaction using the array of all possible reaction rates
    % BUT if next reaction time is greater than Tmax, don't execute it,
    % only update with growth.  this is denoted with reaction_id = 0.
    if t+delta_t <= Tmax
        
        % collect a random number.  If we've run out of pregenerated ones,
        % generate some new ones.
        if n + 1 < number_of_random_numbers_to_pre_generate
            random_number = lots_of_random_numbers(n+1);
            n = n + 1;
        else
            disp('gac: ran out of pre-generated rns, generating new ones')
            
            % generate new ones
            lots_of_random_numbers = rand(1,number_of_random_numbers_to_pre_generate);
            
            % collect the one we need
            random_number = lots_of_random_numbers(1);
            
            % reset counter
            n = 1;
            
        end
        
        % get the id of the chosen reaction
        reaction_id = find(random_number < cumsum(lambda_arr)./lambda_total,1);
    
    else
        reaction_id = 0;
        
    end

    
    % if there's a reaction to happen, call the update routine and update
    % the output arrays
    if ~isempty(reaction_id)
             
        [V_arr,tvec,total_pop_arr,num_clumps_arr] = gac_gillespie_update();
        
    end
    
    % update time
    t = t + delta_t;
    
    % if the last reaction was dropped because it was scheduled to happen
    % after Tmax, cull output arrays down to Tmax
    if reaction_id==0
        
        total_pop_arr = total_pop_arr(tvec<=Tmax);
        num_clumps_arr = num_clumps_arr(tvec<=Tmax);
        tvec = tvec(tvec<=Tmax);
        
    end
         
    % if desired, write time and cluster size array to txt file
    if lwritetxt
        
        % open the file in append mode
        fid = fopen([txtdir filesep txtname],'a');
        
        % data to be written
        outarr = [t,V_arr];
        
        % format string for data
        format_str = [repmat('%f ',1,numel(outarr)-1), '%f\n'];
        
        % write the data
        fprintf(fid,format_str,outarr);
        
        % close the file
        fclose(fid);
        
    end
    
    % if the population went extinct (no more clusters), stop the
    % simulation
    if isempty(V_arr)
        return
    end
    
end
    
    % gillespie update function.  nested to inherit random numbers
    function [V_arr_out,tvec_out,tot_pop_out,num_clumps_out] = gac_gillespie_update()
        
        % collect the label that denotes which reaction is happening
        if reaction_id > 0
            this_label = all_labels(reaction_id);
        else
            this_label  = 0;
            delta_t = Tmax-t;
        end
        
        % create and output variable.  necessary to avoid referencing issues
        % with the nested function.
        V_arr_out = V_arr;
        
        %% growth (deterministic) - do this first
        
        % construct a timestep for numerical integration.  can change
        % resolution here for better or worse accuracy.  Base timestep is
        % .01, if time between reactions gets small, use a smaller step.
        %dt = min(delta_t./10,.01);
        dt = min(delta_t.*fraction_of_delta_t,baseline_dt);
        
        % number of steps taken between reactions
        %numsteps = round((delta_t - dt)/dt);
        numsteps = ceil(delta_t/dt);

        
        % update time array.  new variable name is needed to avoid
        % referencing issues with the nested function.
        tvec_out = [tvec, linspace(t+dt, t+delta_t, numsteps)];
            
        % temporary arrays to be updated during numerical integration
        tmp_tot_pop = zeros(1,numsteps);
        tmp_num_clumps = zeros(1,numsteps);
        
        % loop over time and update according to growth equation
        for s = 1:numsteps
        
            % grow all clusters in one deterministic vectorized step
            V_arr_out = V_arr_out + dt.*lr.*V_arr_out.*(1-sum(V_arr_out)./max_total_pop);
            
            % if clusters die out, remove them
            %V_arr_out = V_arr_out(V_arr_out>=1);
            
            % update tmp_totpop
            tmp_tot_pop(s) = sum(V_arr_out);
            
            % update tmp_num_clumps
            tmp_num_clumps(s) = numel(V_arr_out);
        
        end
        
        % append total population array
        tot_pop_out = [total_pop_arr, tmp_tot_pop];
        
        % append total number of clumps array
        num_clumps_out = [num_clumps_arr, tmp_num_clumps];
        
        %% update cluster size array based on which reaction is happening
        switch this_label
            case 1
                % collapse

                % remove clusters that have collapsed from the cluster array.
                V_arr_out(reaction_id) = [];
                
            case 2
                % aggregation
                              
                % find ids of clusters that will aggregate in this timestep
                agg_id = ids_arr(reaction_id);
                
                % get index of this aggregation reaction in A_mat
                lin_ind_of_this_agg = lin_inds_for_Amat(agg_id);
                
                % backout the ids of clusters involved in this aggregation
                % reaction (row + column of A_mat)
                [agg_row,agg_col] = ind2sub(size(A_mat),lin_ind_of_this_agg);
                
                % one of the clusters increases in size due to aggregation
                V_arr_out(agg_row) = V_arr_out(agg_row) + V_arr_out(agg_col);
                
                % the other one is removed from the array
                V_arr_out(agg_col) = [];
                                          
            case 3
                % sprout//fragment
                
                % get id of cluster that will fragment
                frag_id = ids_arr(reaction_id);
                
                % reduce sprouted clusters by 1
                V_arr_out(frag_id) = V_arr_out(frag_id) - 1;
                
                % add this collection of new single cells to cluster array
                V_arr_out = [V_arr_out, 1];
            case 0
                
                % do nothing
                
        end
                
                
        %% final output
        % update population array based on reaction
        tot_pop_out(end) =  sum(V_arr_out);
        
        % update total number of clumps array based on reaction
        num_clumps_out(end) =  numel(V_arr_out); 
        
                  
    end

   
        
end










% testing_vigil_model.m

numtrials = 1;

numls = 5;
min_log_ls = -1.5; %-8
max_log_ls = -.5; %-4

ls_arr = logspace(min_log_ls,max_log_ls,numls);
%ls_arr = linspace(.1,1,6);

min_log_la = -1.5;  %-1.5
max_log_la = -.5;     %0
numla = 5;
la_arr = logspace(min_log_la,max_log_la,numla);
%la_arr = .01;
%la_arr = linspace(.1,1,6);

lr = 0;
lc = 0;

nu_A = 1;
nu_F = 1;

K = 1e4;

lwritetxt = 1;
txtdir = '/Users/brandonschlomann/Documents/MATLAB/Simulations/aggregates/testing_fraction';


% generate a typical cluster size distribution at carrying capacity
n0 = 1.5 + randn(1,20);
n0(n0<0) = 0;
n0 =  10.^n0;
n0 = n0(n0 < 500);
n0 = [n0, (1e4-sum(n0))];

%dt = .001.*log10(la_arr)./(log10(la_arr(end)));
dt_arr = .5./(K.*la_arr);
Tmax = 2;
%tvec = 0:dt:Tmax;


all_mean_sizes = zeros(numtrials,1);
all_big_frac = zeros(numtrials,1);

mean_log_size = zeros(numla,numls);

mean_log_big_frac = zeros(numla,numls);

icounter = 0;

%% gac info file
info_file_name = 'gacinfo.txt';

info_file_id = fopen([txtdir filesep info_file_name],'a');

% summary
fprintf(info_file_id,'%s\n','summary of gac simulation');

% clock info
clockinfo = clock;
fprintf(info_file_id,'%s\n\n',num2str(clockinfo));

% aggregation
format_str = [repmat('%f ',1,numel(la_arr)-1), '%f\n'];

fprintf(info_file_id,['%s ', format_str],'A',la_arr);

% fragmentation
format_str = [repmat('%f ',1,numel(ls_arr)-1), '%f\n'];

fprintf(info_file_id,['%s ', format_str],'F',ls_arr);

% growth
format_str = [repmat('%f ',1,numel(lr)-1), '%f\n'];

fprintf(info_file_id,['%s ', format_str],'r',lr);

% dispersal
format_str = [repmat('%f ',1,numel(lc)-1), '%f\n'];

fprintf(info_file_id,['%s ', format_str],'D',lc);

% agg exponent
format_str = [repmat('%f ',1,numel(nu_A)-1), '%f\n'];

fprintf(info_file_id,['%s ', format_str],'nu_A',nu_A);

% frag exponent
format_str = [repmat('%f ',1,numel(nu_F)-1), '%f\n'];

fprintf(info_file_id,['%s ', format_str],'nu_F',nu_F);

% Tmax
format_str = [repmat('%f ',1,numel(Tmax)-1), '%f\n'];

fprintf(info_file_id,['%s ', format_str],'Tmax',Tmax);

% numtrials
format_str = [repmat('%d ',1,numel(numtrials)-1), '%d\n\n'];

fprintf(info_file_id,['%s ', format_str],'numtrials',numtrials);

% header for loop update
fprintf(info_file_id,'%s %s %s %s %s %s %s %s\n','file#','A','F','r','D','nu_A','nu_F','trial#');

% close file
fclose(info_file_id);
    
%% loop
    
for a = 1:numla
    dt = dt_arr(a);
    for s = 1:numls
        disp(['param number' num2str(s)])
        
        for n = 1:numtrials
            
            icounter = icounter + 1;
            
            txtname = ['gacout_' num2str(icounter) '.txt'];
            
            %[V_arr,~] = gac(lr,la_arr(a),lc,ls_arr(s),dt,Tmax,n0,nu_A,nu_F);
            
            [V_arr,~] = gac_gillespie(lr,la_arr(a),lc,ls_arr(s),Tmax,n0,nu_A,nu_F,lwritetxt,txtdir,txtname);
            all_mean_sizes(n) = mean(log10(V_arr));
            
            all_big_frac(n) = sum(V_arr(V_arr==max(V_arr)))./sum(V_arr);
            
            % gacinfo output
            out_arr = [icounter, la_arr(a), ls_arr(s), lr, lc, nu_A, nu_F, n];
            info_file_id = fopen([txtdir filesep info_file_name],'a');
            
            format_str = ['%d %f %f %f %f %f %f %d\n'];

            fprintf(info_file_id,format_str,out_arr);


        end
        
        % abundance
        mean_log_size(a,s) = mean(all_mean_sizes);
        
        % big frac
        mean_log_big_frac(a,s) = mean(all_big_frac);
        
    
    end
end


% %heat map
% figure; hold on;
% contourf(log10(ls_arr),log10(la_arr),mean_log_size);
% %contourf(ls_arr,la_arr,mean_log_size);
% colorbar;
% title(['mean log size, 24hr quench'],'fontsize',24)
% set(gca,'fontsize',24,'linewidth',4)
% xlabel('{\lambda}_s','fontsize',24)
% ylabel('{\lambda}_a','fontsize',24)
% 
% %heat map
% figure; hold on;
% %contourf(log10(ls_arr),log10(la_arr),mean_log_size);
% contourf(log10(ls_arr),log10(la_arr),10.^(mean_log_big_frac));
% colorbar;
% title(['mean gel frac, 24hr quench'],'fontsize',24)
% set(gca,'fontsize',24,'linewidth',4)
% xlabel('{\lambda}_s','fontsize',24)
% ylabel('{\lambda}_a','fontsize',24)


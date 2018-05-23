% testing_vigil_model.m

numtrials = 2;

numls = 4;
min_log_ls = -3; %-8
max_log_ls = -1; %-4

ls_arr = logspace(min_log_ls,max_log_ls,numls);

min_log_la = -2;  %-1.5
max_log_la = 0;     %0
numla = 1;
%la_arr = logspace(min_log_la,max_log_la,numla);
la_arr = .01;

lr = 0;
lc = 0;

nu_A = 1;
nu_F = 1;

K = 1e4;


% generate a typical cluster size distribution at carrying capacity
n0 = 1.5 + randn(1,20);
n0(n0<0) = 0;
n0 =  10.^n0;
n0 = n0(n0 < 500);
n0 = [n0, (1e4-sum(n0))];

%dt = .001.*log10(la_arr)./(log10(la_arr(end)));
dt_arr = .5./(K.*la_arr);
Tmax = 24;
%tvec = 0:dt:Tmax;


all_mean_sizes = zeros(numtrials,1);

mean_log_size = zeros(numla,numls);

for a = 1:numla
    dt = dt_arr(a);
    for s = 1:numls
        disp(['param number' num2str(s)])
        
        for n = 1:numtrials
            
            %[V_arr,~] = gac(lr,la_arr(a),lc,ls_arr(s),dt,Tmax,n0,nu_A,nu_F);
            [V_arr,~] = gac_gillespie(lr,la_arr(a),lc,ls_arr(s),Tmax,n0,nu_A,nu_F);
            all_mean_sizes(n) = mean(log10(V_arr));
        
        end
        
        % abundance
        mean_log_size(a,s) = mean(all_mean_sizes);
        
    
    end
end


% heat map
% figure; hold on;
% contourf(log10(ls_arr),log10(la_arr),mean_log_size);
% colorbar;
% title(['mean log size, 24hr quench'],'fontsize',24)
% set(gca,'fontsize',24,'linewidth',4)
% xlabel('{\lambda}_s','fontsize',24)
% ylabel('{\lambda}_a','fontsize',24)

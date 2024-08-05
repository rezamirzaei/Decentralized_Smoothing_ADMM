% clear;
% close all;
%%
clc;
clear;
  

n                 = 500;   % number of local sample
p                 = 18;    % problem dimension
Num_trials        = 1;
%%%% Initializations %%%%
N=30; %number of nodes

loop=100;
st=1;
 ma=5000;
dist_sub_MCP = zeros(Num_trials,ma);
dist_sub_SCAD = zeros(Num_trials,ma);
dist_sub_L1 = zeros(Num_trials,ma);
dist_siad_MCP = zeros(Num_trials,ma);
dist_siad_SCAD = zeros(Num_trials,ma);

dist_mean_siad_SCAD = zeros(Num_trials,ma);
dist_mean_siad_MCP = zeros(Num_trials,ma);


dist_mean_sub_SCAD = zeros(Num_trials,ma);
dist_mean_sub_MCP = zeros(Num_trials,ma);
dist_mean_sub_L1 = zeros(Num_trials,ma);


acc_sub_MCP = zeros(Num_trials,ma);
acc_sub_SCAD = zeros(Num_trials,ma);
acc_sub_L1 = zeros(Num_trials,ma);
acc_siad_MCP = zeros(Num_trials,ma);
acc_siad_SCAD = zeros(Num_trials,ma);


tau=0.75;
cof=0.2;
s=5;
a_SCAD=3.7;
a_MCP=2.4;
lambda= 0.055;
trial=1;
ex=1;
%%

for trial = 1:3:Num_trials
for ex=1:loop
%%
    snr=trial;
    rng(snr*100+ex,'philox')
    seed = randperm(10000);
    option.seed = seed(1:10000);
    fprintf('trial: %d\n',trial);
    fprintf('expriement: %d\n',ex);
    %[a,b,x_true] = generate_data_phase_sparse_gmm(Num_Nodes,m,n,ceil(n*3/10),10^((-(snr-st)/20)));
    %[a,b,x_true] = generate_data_phase_sparse_gmm(Num_Nodes,m,n,ceil(n*3/10),snr-st);
    %[a,b,x_true] = generate_data_phase_sparse_inside_noise(Num_Nodes,m,n,ceil(n*3/10),snr-st);
    [X,y,x_true,active,maxe] = generate_data_local_dist(n,p,N,3,cof,tau);
 %   [pp,stats]=quantreg(X,y,tau); 

 option.p=0.3;
 option.top='iot';
 Num_Nodes=N;
  if strcmp(option.top,'star')
     net = sign(undirected_star_generator(Num_Nodes))-eye(Num_Nodes);
    lam = svds(net,2);
    disp(['sigma: ', num2str(min(lam(2)))]);
    elseif strcmp(option.top,'ring')
        net = sign(undirected_ring_generator(Num_Nodes))-eye(Num_Nodes);
        lam = svds(net,2);
        disp(['sigma: ', num2str(min(lam(2)))]);
  elseif strcmp(option.top,'iot')
       net = sign(undirected_graph_generator_iot(Num_Nodes));
        lam = svds(net,2);
        disp(['sigma: ', num2str(min(lam(2)))]);
  else
        net = sign(undirected_graph_generator(Num_Nodes, option.p,rand(1)))-eye(Num_Nodes);
        lam = svds(net,2);
        disp(['sigma: ', num2str(min(lam(2)))]);
  end   


 
    
    
    %% SIAD
    %SCAD
     fprintf('DSAD\n');
    option.rho_index=0;%SCAD
      option.net=net;
    option.beta_init=zeros(p,1);
    option.max_iters_outer=100;
    option.max_iters_inner=10000;
    option.max_it=ma;
    option.rho=1;
    option.co=1;
    option.p=0.3;
    w=max(tau,1-tau)*maxe+lambda*n;
    option.w=w;
    option.d=sqrt(20)*w;
    option.c=sqrt(3/2);
    option.beta=1;
    option.type=1/2;
    [beta_siad_SCAD,dist_siad_scad,dist_mean_siad_scad,acc_rec_siad_scad] = smoothing_ADMM(X,y,tau,a_SCAD,lambda,x_true,active,maxe,option);
    dist_siad_SCAD(trial,:) = dist_siad_SCAD(trial,:) + dist_siad_scad'/loop;
    dist_mean_siad_SCAD(trial,:) = dist_mean_siad_SCAD(trial,:) + dist_mean_siad_scad'/loop;
    acc_siad_SCAD(trial,:) = acc_siad_SCAD(trial,:) + acc_rec_siad_scad'/(loop);
    fprintf('dist %d\n',dist_siad_scad(ma))

    option.rho_index=1;%SCAD
      option.net=net;
    option.beta_init=zeros(p,1);
    option.max_iters_outer=100;
    option.max_iters_inner=10000;
    option.max_it=ma;
    option.rho=1;
    option.co=1;
    option.p=0.3;
    w=max(tau,1-tau)*maxe+lambda*n;
    option.w=w;
    option.d=sqrt(20)*w;
    option.c=sqrt(3/2);
    option.beta=1;
    option.type=1/2;
    [beta_siad_MCP,dist_siad_mcp,dist_mean_siad_mcp,acc_rec_siad_mcp] = smoothing_ADMM(X,y,tau,a_MCP,lambda,x_true,active,maxe,option);
    dist_siad_MCP(trial,:) = dist_siad_MCP(trial,:) + dist_siad_mcp'/loop;
     dist_mean_siad_MCP(trial,:) = dist_mean_siad_MCP(trial,:) + dist_mean_siad_mcp'/loop;
    acc_siad_MCP(trial,:) = acc_siad_MCP(trial,:) + acc_rec_siad_mcp'/(loop);
    fprintf('dist %d\n',dist_siad_scad(ma))

    %% SIAD
    %SCAD
     fprintf('DSAD\n');
    option.rho_index=0;%SCAD
    option.beta_init=zeros(p,1);
    option.max_iters_outer=100;
    option.max_iters_inner=10000;
    option.max_it=ma;
    option.rho=1;
    option.co=1;
    [beta_sub_SCAD,dist_sub_scad,dist_mean_sub_scad,acc_rec_sub_scad] = dis_subgradient_quantile(X,y,tau,a_SCAD,lambda,x_true,active,option);
    dist_sub_SCAD(trial,:) = dist_sub_SCAD(trial,:) + dist_sub_scad'/loop;
     dist_mean_sub_SCAD(trial,:) = dist_mean_sub_SCAD(trial,:) + dist_mean_sub_scad'/loop;
    acc_sub_SCAD(trial,:) = acc_sub_SCAD(trial,:) + acc_rec_sub_scad'/(loop);
    fprintf('dist %d\n',dist_sub_scad(ma));

   
    %%
    % option.rho_index=1;%SCAD
    fprintf('DSPQ\n');
     option.rho_index=1;%MCP
    option.beta_init=zeros(p,1);
    option.max_iters_outer=100;
    option.max_iters_inner=10000;
    option.max_it=ma;
    option.rho=1;
    option.co=1;
    [beta_sub_mcp,dist_sub_mcp,dist_mean_sub_mcp,acc_rec_sub_mcp] = dis_subgradient_quantile(X,y,tau,a_MCP,lambda,x_true,active,option);
    dist_sub_MCP(trial,:) = dist_sub_MCP(trial,:) + dist_sub_mcp'/loop;
    dist_mean_sub_MCP(trial,:) = dist_mean_sub_MCP(trial,:) + dist_mean_sub_mcp'/loop;
    acc_sub_MCP(trial,:) = acc_sub_MCP(trial,:) + acc_rec_sub_mcp'/(loop);
    fprintf('dist %d\n',dist_sub_mcp(ma));


     fprintf('DSPQ\n');
     option.rho_index=2;%SCAD
    option.beta_init=zeros(p,1);
    option.max_iters_outer=100;
    option.max_iters_inner=10000;
    option.max_it=ma;
    option.rho=1;
    option.co=1;
    [beta_sub_l1,dist_sub_l1,dist_mean_sub_l1,acc_rec_sub_l1] = dis_subgradient_quantile(X,y,tau,0,lambda,x_true,active,option);
    dist_sub_L1(trial,:) = dist_sub_L1(trial,:) + dist_sub_l1'/loop;
    dist_mean_sub_L1(trial,:) = dist_mean_sub_L1(trial,:) + dist_mean_sub_l1'/loop;
    acc_sub_L1(trial,:) = acc_sub_L1(trial,:) + acc_rec_sub_l1'/(loop);
    fprintf('dist %d\n',dist_sub_l1(ma));
    
%     plot(10*log10(dist_siad_mcp),'linewidth',4);
%    hold on;
%    plot(10*log10(dist_siad_scad),'linewidth',4);

end
end
%%
  figure;
   plot(10*log10(dist_siad_SCAD),'linewidth',4);
   hold on;
   plot(10*log10(dist_siad_MCP),'linewidth',4);
   hold on;
   plot(10*log10(dist_sub_SCAD),'linewidth',4);
   hold on;
   plot(10*log10(dist_sub_MCP),'linewidth',4);
   hold on;
   plot(10*log10(dist_sub_L1),'linewidth',4);
    set(gca,'FontSize',20);
    legend(['DSAD-SCAD'], ...
        ['DSAD-MCP'], ...
        ['DSPQ-SCAD'], ...
        ['DSPQ-MCP'], ...
   ['$l_1$-dQR'],'location','Best','interpreter','latex');
 xlabel('Iterations (k)');  ylabel(['MSE (dB)'],'interpreter','latex');
 %%
  plot((acc_siad_SCAD(2:500)), 'linewidth', 4);
hold on;
plot((acc_siad_MCP(2:500)), 'linewidth', 4);
hold on;
plot((acc_sub_SCAD(2:500)), 'linewidth', 4);
hold on;
plot((acc_sub_MCP(2:500)), 'linewidth', 4);
hold on;
plot((acc_sub_L1(2:500)), 'linewidth', 4);
set(gca, 'FontSize', 20);
   set(gca,'FontSize',20);
       legend(['DSAD-SCAD'], ...
        ['DSAD-MCP'], ...
        ['DSPQ-SCAD'], ...
        ['DSPQ-MCP'], ...
   ['$l_1$-dQR'],'location','Best','interpreter','latex');
 xlabel('Iterations (k)');  ylabel(['Recognition accuracy'],'interpreter','latex');
 %%
   figure;
   plot(10*log10(dist_mean_siad_SCAD(2:ma)),'linewidth',4, 'LineStyle', '-', 'Marker', 'o');
   hold on;
   plot(10*log10(dist_mean_siad_MCP(2:ma)),'linewidth',4,'LineStyle', '--', 'Marker', '+');
   hold on;
   plot(10*log10(dist_mean_sub_SCAD(2:ma)),'linewidth',4, 'LineStyle', ':', 'Marker', 's');
   hold on;
   plot(10*log10(dist_mean_sub_MCP(2:ma)),'linewidth',4,'LineStyle', '-.', 'Marker', 'd');
   hold on;
   plot(10*log10(dist_mean_sub_L1(2:ma)),'linewidth',4,'LineStyle', '-', 'Marker', 'x');
       set(gca,'FontSize',20);
    legend(['DSAD-SCAD'], ...
        ['DSAD-MCP'], ...
        ['DSPQ-SCAD'], ...
        ['DSPQ-MCP'], ...
   ['$l_1$-dQR'],'location','Best','interpreter','latex');
 xlabel('Iterations (k)');  ylabel(['Network MSE (dB)'],'interpreter','latex');

%%



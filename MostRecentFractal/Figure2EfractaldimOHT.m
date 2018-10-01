%% calculate fractal dimension as Rg vs cluster number

% cleanup
close all
clear all
clc

% start timing
tic

 cmapO=cbrewer('seq','Oranges',20);

% % minimum number of cells per cluster
% cluster_definition=2;

% cutoff distance
threshold = 75;

%   load OHTAssayStoreClusters.mat
%   load OHTAssayStoreClusters3exp.mat
  load OHTAssayStoreClustersE6.mat
%  AssayOHT = [1157 54 1210 51 70  1151  1252];
% 
% CellDens = [1000 1000 1000 500  500 500 500 500];
%   color = [12 12 12 5 5 5 5 5];
%   MaxTime= 4*[28 30 40 60 60 60 39 27];
% 
%  AssayOHT = [ 54 67 129 1210 51 70 1251  1252];
% 
% CellDens = [ 1000 1000 1000 1000 500 500 500 500];
%   color = [ 12 12 12 12 5 5 5 5];
%   MaxTime= 4*[ 30 48 39 40 60 60 32 27];
  
%    AssayOHT = [3 22 6 19 9 16 12 13];
% 
% CellDens = [500 500 1000 1000 2000 2000 3000 3000];
%   color = [5 5 12 12 16 16 20 20];
%   MaxTime= 4*[60 60 60 60 36 32 24 27];

%     AssayOHT = [3 19 16 13];
%  CellDens = [500 1000 2000 3000];
%    color = [5 12 16 20];
%    MaxTime= 4*[60 60 27 22];

    AssayOHT = [51 54 64 61];
% 
 CellDens = [500 1000 2000 3000];
   color = [5 12 16 20];
   MaxTime= 4*[60 35 24 20];



% now plot distribution of cluster sizes over this time range
%subplot(2,5,plotCount+5)

KeepIndex = nan(1,numel([storeClusterInfo.Rg]));

storeRg = [storeClusterInfo.Rg];
storeWellNum = [storeClusterInfo.WellNum];
storeTime = [storeClusterInfo.TimeIndex];
storeCellNum = [storeClusterInfo.CellNum];

        
        
 %%plot in reverse order for lowest density last
for i = numel(AssayOHT):(-1):1
        
    fWell = find(storeWellNum == AssayOHT(i));
    
    ft = find(storeTime(fWell) < MaxTime(i));
    % for gray marker [0.8 0.8 0.8]
    plot(storeCellNum(fWell(ft)),...
        storeRg(fWell(ft)),...
        'o',...
        'MarkerEdgeColor',[0.8 0.8 0.8],'LineWidth',0.2,...
        'MarkerFaceColor',cmapO(color(1,i),:),...
        'MarkerSize',7,...
        'DisplayName',num2str(CellDens(i)));
    
    xlabel('Cluster Size')
    ylabel('R_g (\mum)')
    hold on
    
    KeepIndex(fWell(ft)) = 1;
    
end

%legend('show')
% labels = get(legend(), 'String');
% plots = flipud(get(gca, 'children'));
% neworder = [3 4 2 1];
% 
% legend(plots(neworder), labels(neworder),'location','northwest','FontSize',8);
set(gca,'XScale','log','YScale','log')

fKeep = find(~isnan(KeepIndex));

% fit Rg to cluster mass as a power law
f = fit(storeCellNum(fKeep)',...
    storeRg(fKeep)','power1')

% show fit parameters
t=text(100,50,strcat('D_f=',num2str(1/f.b)));

% plot best fit estimate
X = 1:500;
powerFit = (f.a).*X.^(f.b);

plot(X,powerFit,'k-','Linewidth',2);

set(gca,'XScale','log','YScale','log')
ylim([5 inf]);
ax=gca;
box on
grid off
            ax.XColor='black';
            ax.YColor='black';
            set(gca,'fontsize',10);
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperPosition', [0 0 10 7.5])
title('Fractal Dimension: DMSO low EGF 51 54 64 61','FontSize',10);
figname='2d Fractal Dimension: OHT densities';
% saveas(gcf,figname,'epsc');
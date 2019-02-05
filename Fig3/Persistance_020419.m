function Persistance
%Quantify persistance of follower movement

close all; clearvars; clc

LeaderData = csv2struct('Z:\ENG_BBCancer_Shared\group\Swarming Paper Final Data\KymoAnalysis\CombineReps\KymoLeaderSpeedCombined.csv');
FollowerData = csv2struct('Z:\ENG_BBCancer_Shared\group\Swarming Paper Final Data\KymoAnalysis\CombineReps\KymoFollowerSpeedCombined.csv');
ZigData = csv2struct('Z:\ENG_BBCancer_Shared\group\Swarming Paper Final Data\KymoAnalysis\CombineReps\KymoFollowerSpeedCombinedZig.csv');
DMSOwells = [6 19 1257 1258];
OHTwells = [54 67 129 1210];
cmap = cbrewer('qual','Dark2',3);
cmap2 = cbrewer('qual','Set1',3);

%make a beeswarm plot comparing persistance 

%first extract persistance data from the leaders 
ixOHT = ismember(LeaderData.Well,OHTwells);
LeaderDur = LeaderData.Duration_hrs_(ixOHT);

%now extract persistance for followers
ixZ = ismember(ZigData.Well,OHTwells);
FollowerDur = ZigData.Duration_hrs_(ixZ);

%extract general follower zig-zag length
ixY = ismember(FollowerData.Well,OHTwells);
FDY = FollowerData.Duration_hrs_(ixY);

figure;
handles = plotSpread({LeaderDur,FollowerDur},'categoryColors',{'r','b'});
set(handles{1},'MarkerSize',15)
set(handles{1}(1),'color',cmap2(2,:))
set(handles{1}(2),'color',cmap2(3,:))
% set(handles{1}(3),'color',cmap2(3,:))
set(gca,'XTick',1:2,'XTickLabel',{'Leaders','Selected Non-Leaders'})
title('Persistance (time)')
ylabel('Duration (hrs)')
title('Leading Duration')
plot([0 5],[median(LeaderDur),median(LeaderDur)],'--','Color',cmap2(2,:),'linewidth',2);
plot([0 5],[median(FollowerDur),median(FollowerDur)],'--','Color',cmap2(3,:),'linewidth',2);
set(gca,'LineWidth',1)
ax= gca;
ax.XColor='black';
ax.YColor='black';
ax.YGrid = 'on';
ax.LineWidth = 1;
set(gca,'fontsize',8);
set(gca,'YMinorTick','on')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 6.2 10.57/2]);

figname = char(strcat('Z:\ENG_BBCancer_Shared\group\0Zach\Summer 2018\NewKymoBeeswarm (8-09)','\L-Fduration080918'));
% saveas(gcf,figname,'epsc');
% saveas(gcf,figname,'tif');





%compile the speeds for all three groups
ix = ismember(LeaderData.Well,OHTwells); 
LeaderSPD = LeaderData.Speed_microns_per_hr_(ix);
ix = ismember(FollowerData.Well,OHTwells);
FollowerSPD = FollowerData.Speed_microns_per_hr_(ix);
ix = ismember(ZigData.Well,OHTwells);
ZigSPD = ZigData.Speed_microns_per_hr_(ix);

figure;
handles = plotSpread({LeaderSPD,FollowerSPD});
set(handles{1},'MarkerSize',20)
set(handles{1}(1),'color',cmap2(2,:))
set(handles{1}(2),'color',cmap2(3,:))
% set(handles{1}(3),'color',cmap2(3,:))
set(gca,'XTick',1:2,'XTickLabel',{'Leaders','Selected Non-Leaders'})
title('Leading Speed')
set(gca,'LineWidth',1)
ax= gca;
ax.XColor='black';
ax.YColor='black';
ax.YGrid = 'on';
set(gca,'fontsize',8);
set(gca,'LineWidth',1)
figname = char(strcat('Z:\ENG_BBCancer_Shared\group\0Zach\Summer 2018\Updated Graphics (6-30)','\L-Fspeed063018'));
% saveas(gcf,figname,'epsc');
% saveas(gcf,figname,'tif');




end
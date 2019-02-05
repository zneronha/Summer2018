function KYMODATAPLOTTERfollowers063018
close all; clearvars; clc

%path must be in group 0Zach EGFreplicatedata SelectedData

Data = csv2struct('Z:\ENG_BBCancer_Shared\group\Swarming Paper Final Data\KymoAnalysis\CombineReps\KymoLeaderSpeedCombined.csv');
FollowerData = csv2struct('Z:\ENG_BBCancer_Shared\group\Swarming Paper Final Data\KymoAnalysis\CombineReps\KymoFollowerSpeedCombined.csv');
DMSOwells = [6 19 1257 1258];
OHTwells = [54 67 129 1210];
cmap = cbrewer('qual','Dark2',3);
cmap2 = cbrewer('qual','Set1',3);

%DMSO
RepMat = [19	1257;...
19	1257;...
19	1257];

%compile a plot of the E6 E12 Velocity Data the usual way
plotspeed_shaderr(RepMat);

%find all the DMSO wells
ixDMSO = ismember(Data.Well,DMSOwells);
SpeedMAT = Data.Speed_microns_per_hr_; SpeedMAT = SpeedMAT(ixDMSO);
TimeMAT = Data.Frame; TimeMAT = TimeMAT(ixDMSO)./4;
hold on; 
ap = scatter(TimeMAT, SpeedMAT,12,'filled','s');
ap.MarkerFaceColor = cmap2(2,:);
title('DMSO')

%add the follower cells here too!
FixDMSO = ismember(FollowerData.Well,DMSOwells);
FSpeedMAT = FollowerData.Speed_microns_per_hr_; FSpeedMAT = FSpeedMAT(FixDMSO);
FTimeMAT = FollowerData.Frame; FTimeMAT = FTimeMAT(FixDMSO)./4;
hold on; 
bp = scatter(FTimeMAT, FSpeedMAT,10,'filled');
bp.MarkerFaceColor = cmap2(3,:);

ylim([0 60]);


figname = 'Z:\ENG_BBCancer_Shared\group\0Zach\Summer 2018\Kymo Results (6-20)\KymoResultsGraphics\DMSO_overlay071618';
ax= gca;
ax.XColor='black';
ax.YColor='black';
ax.YGrid = 'on';
set(gca,'fontsize',8);
set(gca,'XMinorTick','off','YMinorTick','off')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 6.2575 10.575/2]);

saveas(gcf,figname,'epsc');
print('-dtiff','-r1000',figname)

figure;
%Now for the OHT
RepMat = [54	1210;...
54	1210;...
54	1210];

%compile a plot of the E6 E12 Velocity Data the usual way
plotspeed_shaderr(RepMat);

%find all the OHT wells
ixOHT = ismember(Data.Well,OHTwells);
SpeedMAT = Data.Speed_microns_per_hr_; SpeedMAT = SpeedMAT(ixOHT);
TimeMAT = Data.Frame; TimeMAT = TimeMAT(ixOHT)./4;
hold on; 
ap = scatter(TimeMAT, SpeedMAT,12,'filled','s');
ap.MarkerFaceColor = cmap2(2,:);
title('OHT')

%now for the followers
FixOHT = ismember(FollowerData.Well,OHTwells);
FSpeedMAT = FollowerData.Speed_microns_per_hr_; FSpeedMAT = FSpeedMAT(FixOHT);
FTimeMAT = FollowerData.Frame; FTimeMAT = FTimeMAT(FixOHT)./4;
bp = scatter(FTimeMAT, FSpeedMAT,10,'filled');
bp.MarkerFaceColor = cmap2(3,:);

ylim([0 60]);


figname = '/Users/Zach/Desktop/OHToverlay063018';
ax= gca;
ax.XColor='black';
ax.YColor='black';
ax.YGrid = 'on';
set(gca,'fontsize',8);
set(gca,'XMinorTick','off','YMinorTick','off')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 6.2575 10.575/2]);

% saveas(gcf,figname,'epsc');
% print('-dtiff','-r1000',figname)
% print('-depsc','-r1000',figname)

end

function plotspeed_shaderr(RepMat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot cell proliferation, speeds from all data files
%Updated June 1 to make two seperate figures


% clean up workspace


FramesToHours = 4;
TimeWindow = FramesToHours; % analyze this frequently (i.e. hours)
ImageArea = .900*.900; % mm^2
% FigNameDMSO='DMSO 1000c ProlifPlot E18';
% FigNameOHT='OHT 1000c ProlifPlot E18';
%% read in file names
DirectoryFileList = dir;
FileNames = {DirectoryFileList.name};
FileTypeList = ones(length(FileNames),1);

% pick out .mat files
for i = 1:length(FileNames)
    [filepath,name,ext] = fileparts(FileNames{i});
    
    if ~strcmp(ext,'.mat')
        FileTypeList(i) = 0;
    end
end

f = find(FileTypeList == 0);
FileNames(f) = [];% allocate struct for data
storeCellInfo = struct('FileName',{},...
                       'TimeIndex',{},...
                       'CellCount',{},...
                       'storeX',{},...
                       'storeY',{},...
                       'velX',{},...
                       'velY',{},...
                       'velR',{});


%% now read in data and analyze
for i = 1:length(FileNames)
    
    load(FileNames{i}); % load .mat file
    
    % count cells per frame
    findCellsPresent = ~isnan(storeX); % what is not nan
    tempCellCount = sum(findCellsPresent,1); % sum over each time
    
    maxTime = size(storeX,2);
    
    % calculate speeds
    velX=storeX(:,(1+TimeWindow):1:end)-...
            storeX(:,1:1:end-TimeWindow);
     
    velY=storeY(:,(1+TimeWindow):1:end)-...
            storeY(:,1:1:end-TimeWindow);
    
    velR = sqrt(velX.*velX + velY.*velY);
            
    % fill in struct every hour
    tempCellInfo = struct('FileName',FileNames{i},...
        'TimeIndex',{((1+TimeWindow):1:maxTime)/FramesToHours},...
        'CellCount',{tempCellCount((1+TimeWindow):1:end)},...
        'storeX',{storeX(:,(1+TimeWindow):1:end)},...
        'storeY',{storeY(:,(1+TimeWindow):1:end)},...
        'velX',{velX},...
        'velY',{velY},...
        'velR',{velR});
    
    storeCellInfo = [storeCellInfo; tempCellInfo]; % store struct 
    
    % cleanup
    clear storeX storeY storevelX storevelY 
    clear velX velY velR findCellsPresent tempCellCount tempCellInfo
    
end

storeGMOHTCellCount = [];
storeGMOHTTime = [];
storeAssayOHTCellCount = [];
storeAssayOHTTime = [];
% storeGefOHTCellCount = [];
% storeGefOHTTime = [];

%% consolidate OHT Main Figure 

% clearvars -except E6repsdens

condition77{1,1} = RepMat;





% for i = 1:length(FileNames);
%     
%     % cleanup filename for display
%     tempName = strrep(storeCellInfo(i).FileName,'_','');
%     tempName = strrep(tempName,'.mat','');
%     
%     % normalize cell counts before combining
%     NormCellCount = storeCellInfo(i).CellCount./...
%         storeCellInfo(i).CellCount(1);
    
    
    ConditionNameMat = {'DMSO','OHT'};
    for totalcounter = 1:numel(condition77)
%         close all
    RepMat = condition77{totalcounter,1};
    
    mycounter = 0;
%     totalcounter = 1;
    
    storeVelData3 = cell(size(RepMat,1),1);
    
    loopcounter = 1;
    qualcheck2 = 0;
    for uu = 1:size(RepMat,1)
        mycounter = mycounter + 1;
%         if mycounter == 4
%        figure('Name','-EGF +OHT Counts','NumberTitle','off')
    
        conditionstore = [];
        rawconditionstore = [];
        velstore = [];
     
           for vv = 1:size(RepMat,2)
               mytempwell = RepMat(uu,vv);

               index = find(strcmp({storeCellInfo.FileName}, strcat('EGF(E6)w',num2str(mytempwell),'.mat'))==1);
%convert cell count to density per image area
               cellcount3 = (storeCellInfo(index).CellCount)./ImageArea;
               normcellcount3 = cellcount3./cellcount3(1);
               timeindex = storeCellInfo(index).TimeIndex;
               
               cellvel3 = storeCellInfo(index).velR;
              cellvel3=nanmean(cellvel3,1);
               velstore = [velstore;cellvel3];

              conditionstore = [conditionstore;normcellcount3];
              rawconditionstore = [rawconditionstore;cellcount3];
           end
           
       cellvel4 = nanmean(velstore);
       velstd = nanstd(velstore);
           
       cellcount4 = mean(conditionstore);
       rawcellcount4 = mean(rawconditionstore);
       stdmat = std(conditionstore);
       rawstdmat = std(rawconditionstore);
       

         storeVelData3{uu} = velstore;       
         cmap = cbrewer('qual','Dark2',5); 
%         cmap = cbrewer('qual','Set1',2);

        if mycounter == 1 || mycounter == 3
        else
        boundedline(timeindex,cellvel4,velstd,...
                'cmap',cmap(mycounter,:),'alpha','transparency',0.3); hold on
            box on
%             ylim([0 54])            
            xlim([0 62])
            xlabel('Time (h)','fontsize',9)
            %set(gca,'Ytick',0:2.5:15)
            ylabel('Cell Speed (um/h)','fontsize',9)
%             ylim([0 54]);
            ax= gca;
            ax.XColor='black';
            ax.YColor='black';
            ax.YGrid = 'on';
            box on
%             ax.YMinorGrid = 'on';
%             ax.MinorGridLineStyle = '-';
    %        set(gca,'XMinorTick','off','YMinorTick','on')
           
            set(gca,'fontsize',8);
            
%              title({wellnums;'';'-EGF +OHT Speeds'});
%             title(strcat('Cell Speed over Time',{' '},ConditionNameMat(totalcounter)),'fontsize',11);
% title(strcat('',{' '},ConditionNameMat(totalcounter)),'fontsize',11);
            %legend({'Gefitinib','High EGF','Low EGF'},'Location','northeast','EdgeColor','k','Fontsize',8);

%             title('-EGF +OHT Speeds');
%            set(gcf, 'Position', [0, 0, 300, 400])

%             title(strcat('Cell Speed vs Density',{' '},ConditionNameMat(totalcounter)),'fontsize',11);
% title(strcat('',{' '},ConditionNameMat(totalcounter)),'fontsize',11);
      % legend({'Gefitinib','High EGF','Low EGF'},'Location','northeast','EdgeColor','k','Fontsize',8);
           
       hold on
        end
              j=totalcounter-1;
%          figname=sprintf('%s_%d','speed_E6nE12_shaderr',j);
       
%          (gcf, 'PaperUnits', 'centimeters');
%          set(gcf, 'PaperPosition', [0 0 6.2 10.57])
         
    end
             
    end
end

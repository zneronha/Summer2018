function plotspeed_shaderr_PVALrollingSepCRV_AGM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot cell proliferation, speeds from all data files
%Updated June 1 to make two seperate figures


% clean up workspace
clc; close all;

FramesToHours = 1;
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
%     velX=(storeX(:,(1+TimeWindow):1:end)-...
%             storeX(:,1:1:end-TimeWindow));
%      
%     velY=(storeY(:,(1+TimeWindow):1:end)-...
%             storeY(:,1:1:end-TimeWindow));
    velX=4*(storeX(:,(1+TimeWindow):1:end)-...
            storeX(:,1:1:end-TimeWindow));
     
    velY=4*(storeY(:,(1+TimeWindow):1:end)-...
            storeY(:,1:1:end-TimeWindow));
    if size(velX,2)+1==size(storeX,2)
    else
        a239 = 3;
    end
    

    
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
    disp(FileNames{i})
end

storeGMOHTCellCount = [];
storeGMOHTTime = [];
storeAssayOHTCellCount = [];
storeAssayOHTTime = [];
% storeGefOHTCellCount = [];
% storeGefOHTTime = [];

%% consolidate OHT Main Figure 

% clearvars -except E6repsdens
RepMat = [53	1159	1211;...
54	1157	1210;...
52	1155	128];
condition77{1,1} = RepMat;

RepMat = [5	1111 1259;...
19	1109	1257;...
4	1107	1256];
condition77{2,1} = RepMat;

RepMat = [71	1153	1254;...
51	1151	1251;...
49	1149	1250];
condition77{3,1} = RepMat;

RepMat = [2	1105	125;...
3	1103	124;...
24	1101	121];
condition77{4,1} = RepMat;



% for i = 1:length(FileNames);
%     
%     % cleanup filename for display
%     tempName = strrep(storeCellInfo(i).FileName,'_','');
%     tempName = strrep(tempName,'.mat','');
%     
%     % normalize cell counts before combining
%     NormCellCount = storeCellInfo(i).CellCount./...
%         storeCellInfo(i).CellCount(1);
    
    
    ConditionNameMat = {'OHT1000','DMSO1000','OHT500','DMSO500'};
    for totalcounter = 1:numel(condition77)
        %close all
        figure;
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
        subplot(2,1,1)
        if mycounter == 1
        else
        abc = boundedline(timeindex,cellvel4,velstd,...
                'cmap',cmap(mycounter,:),'alpha','transparency',0.3); hold on
            abc.LineWidth = 1.5;
            box on
            ylim([0 54])            
            xlim([0 62*4/TimeWindow])
            xlabel('Time (h)','fontsize',9)
            %set(gca,'Ytick',0:2.5:15)
            ylabel('Cell Speed (um/h)','fontsize',9)
            ylim([0 54]);
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
title(strcat('',{' '},ConditionNameMat(totalcounter)),'fontsize',11);
            %legend({'Gefitinib','High EGF','Low EGF'},'Location','northeast','EdgeColor','k','Fontsize',8);
            subplot(2,1,2)
            fixMATRIX = [rawcellcount4;cellvel4;velstd];
            fixMATRIX = sortrows(fixMATRIX',1)';
            rawcellcount4 = fixMATRIX(1,:); cellvel4 = fixMATRIX(2,:); velstd = fixMATRIX(3,:);
        cba = boundedline(rawcellcount4,cellvel4,velstd,...
                'cmap',cmap(mycounter,:),'alpha','transparency',0.3); hold on    
            cba.LineWidth = 1.5;
            box on
           % xlim([0 62])
            ylim([0 54]);
            xlabel('Cell Density (cells/mm^{2})','fontsize',9)
            ylabel('Cell Speed (um/h)','fontsize',9)
            ax= gca;
            ax.XColor='black';
            ax.YColor='black';
            ax.YGrid = 'on';
%             ax.YMinorGrid = 'on';
%             ax.MinorGridLineStyle = '-';
            set(gca,'fontsize',8);
            set(gca,'XMinorTick','off','YMinorTick','off')
            
%             title('-EGF +OHT Speeds');
%            set(gcf, 'Position', [0, 0, 300, 400])

%             title(strcat('Cell Speed vs Density',{' '},ConditionNameMat(totalcounter)),'fontsize',11);
title(strcat('',{' '},ConditionNameMat(totalcounter)),'fontsize',11);
      % legend({'Gefitinib','High EGF','Low EGF'},'Location','northeast','EdgeColor','k','Fontsize',8);
           
       hold on
        end
              j=totalcounter-1;
%          figname=sprintf('%s_%d','speed_E6nE12_shaderr',j);
       
         set(gcf, 'PaperUnits', 'centimeters');
         set(gcf, 'PaperPosition', [0 0 6.2 10.57])
         
    end
    %generate and store significance data
    [H1,P1] = ttest(storeVelData3{1},storeVelData3{2});
    [H2,P2] = ttest(storeVelData3{1},storeVelData3{3});
    [H3,P3] = ttest(storeVelData3{2},storeVelData3{3});
    SignificanceStore = [P1;P2;P3];
    subplot(2,1,1);
    ylim([-5 54])
    x1 = double(H1).*timeindex; x1(x1==0) = NaN;
    x2 = double(H2).*timeindex; x2(x2==0) = NaN;
    x3 = double(H3).*timeindex; x3(x3==0) = NaN;
    x1 = x1(1:4:end);
    x2 = x2(1:4:end);
    x3 = x3(1:4:end);
    
    
    hold on
    
% %     t7 = [0 60]; t8 = [2 2];
% %     P=plot(t7,t8);
% %     set(P,'color',[0.5 0.5 0.5],'LineStyle','--','linewidth',1.5);
%     plot(x1,ones(1,size(x1,2)).*0.5*4,'^','MarkerSize',2,'color','k');
    plot(x3,x3.*0-2,'.','MarkerSize',7,'color','k');
%     plot(x3,ones(1,size(x3,2)).*-0.5*4,'o','MarkerSize',2,'color','k');
    
%     xlswrite(char(strcat('Z:\ENG_BBCancer_Shared\group\0Zach\Summer 2018\Updated Graphics (6-30)\Proliferation Plots\PValVel63018GFGM',ConditionNameMat(totalcounter))),SignificanceStore);
    
      figname=char(strcat('Z:\ENG_BBCancer_Shared\group\0Zach\Summer 2018\Replicate Selection (5-11)\Output Graphics 11-10 (difWind)\VelPlotAGM_1frame',ConditionNameMat(totalcounter)));
         set(gcf, 'PaperUnits', 'centimeters');
         set(gcf, 'PaperPosition', [0 0 6.2 10.57]);
         subplot(2,1,1); set(gca,'XMinorTick','on','YMinorTick','on');
         ax = gca;
         ax.LineWidth = 1;
         subplot(2,1,2); set(gca,'XMinorTick','on','YMinorTick','on');
         ax = gca;
         ax.LineWidth = 1;
         set(gcf, 'PaperPosition', [0 0 6.2 10.57]);
         saveas(gcf,figname,'epsc');
         print('-dtiff','-r1000',figname)
             
    end
end

   



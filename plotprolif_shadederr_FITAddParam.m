function plotprolif_shadederr_FITrolling3FitSepCRV_FINALFUN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot cell proliferation, speeds from all data files
%
% 3 May update: include support for fitting the curves to logistic growth
% and exponential models
% 1 June update: split up GFTNB and ASSAY plots

% clean up workspace
clc; clear; close all;

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
    [filepath,name,ext] = fileparts(FileNames{i}); %#ok<ASGLU>
    
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
% storeFitInfo = struct('Well1',{},'Well2',{},...
%                         'ExpA',{},'ExpB',{},...
%                         'ExpR',{},'LogA',{},...
%                         'LogB',{},'LogC',{},...
%                         'LogC',{},'LogR',{});
storeFitInfo = nan(6,12);
loopcounter3 = 1;


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
    disp(FileNames{i})
    
end

storeGMOHTCellCount = [];
storeGMOHTTime = [];
storeAssayOHTCellCount = [];
storeAssayOHTTime = [];
% storeGefOHTCellCount = [];
% storeGefOHTTime = [];


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
        loopcounter3 = 1;
        close all
    RepMat = condition77{totalcounter,1};
        
    

    mycounter = 0;
    storeSIG = cell(size(RepMat,1),1);
   
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
               velstore = [velstore;cellvel3];

              conditionstore = [conditionstore;normcellcount3];
              rawconditionstore = [rawconditionstore;cellcount3];
           end
           
       cellvel4 = nanmean(velstore);
       velstd = std(velstore);
           
       cellcount4 = mean(conditionstore);
       rawcellcount4 = mean(rawconditionstore);
       stdmat = std(conditionstore);
       rawstdmat = std(rawconditionstore);
       
       storeSIG{uu} = conditionstore;
             
         cmap = cbrewer('qual','Dark2',5);  
%         cmap = cbrewer('qual','Set1',2);
        subplot(2,1,1)
        if mycounter == 1
        else
            GCx = cellcount4(2:end)-cellcount4(1:end-1);
            GCx = GCx./cellcount4(2:end);
            blo = plot(timeindex(2:end),GCx,'color',cmap(mycounter,:)); hold on
            
%          blo = boundedline(timeindex,cellcount4,stdmat,...
%                 'cmap',cmap(mycounter,:),'alpha','transparency',0.3); hold on
            blo.LineWidth = 1.5;
        
             box on           
            xlim([0 62])
            xlabel('Time (h)','fontsize',9)
            %set(gca,'Ytick',0:2.5:15)
            ylabel('Normalized Count','fontsize',9)
%             ylim([0 13]);
            ax= gca;
            ax.XColor='black';
            ax.YColor='black';
            ax.YGrid = 'on';


%            ax.YMinorGrid = 'on';
%             ax.MinorGridLineStyle = '-';
    %        set(gca,'XMinorTick','off','YMinorTick','on')
           
            set(gca,'fontsize',8);
            
%             title({wellnums;'';'-EGF +OHT Speeds'});
%             title(strcat('Normalized Cell Count ',{' '},ConditionNameMat(totalcounter)),'fontsize',11);
            title(strcat('',{' '},ConditionNameMat(totalcounter)),'fontsize',11);
%             legend({'High EGF','','Low EGF'},'Location','northwest','EdgeColor','k','Fontsize',8);
            subplot(2,1,2)
            gg(loopcounter3) = boundedline(timeindex,rawcellcount4,rawstdmat,...
                'cmap',cmap(mycounter,:),'alpha','transparency',0.3); hold on   
            a83 = gg(loopcounter3);
            a83.LineWidth = 1.5;
            box on
            xlim([0 62])
            ylim([0 1500]);
            xlabel('Time (hrs)','fontsize',9)
            ylabel('Cell Density (cells/mm^{2})','fontsize',9)
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

%             title(strcat('Cell Dens',{' '},ConditionNameMat(totalcounter)),'fontsize',11);
        title(strcat('',{' '},ConditionNameMat(totalcounter)),'fontsize',11);
%        legend({'High EGF','','Low EGF'},'Location','northwest','EdgeColor','k','Fontsize',8);
           
       hold on
       
        end
        
        loopcounter3 = loopcounter3 + 1;
    end
    %[~,sig12] = ttest(storeSIG{1},storeSIG{2});
    %[~,sig13] = ttest(storeSIG{1},storeSIG{3});
    %[~,sig23] = ttest(storeSIG{2},storeSIG{3});
    %SignificanceStore = [sig12;sig13;sig23];
    %SignificanceStoreH = SignificanceStore<0.05;
    subplot(2,1,1)
    
    %plot significance levels
    %decrease yaxis to allow space for plotting
    %ylim([-1 13])
    %x1 = double(SignificanceStoreH(1,:)).*timeindex; x1(x1==0) = NaN;
    %x2 = double(SignificanceStoreH(2,:)).*timeindex; x2(x2==0) = NaN;
    %x3 = double(SignificanceStoreH(3,:)).*timeindex; x3(x3==0) = NaN;
    %x1 = x1(1:4:end);
    %x2 = x2(1:4:end);
    %x3 = x3(1:4:end);
%     ydu = ones(1,size(SignificanceStoreH,2));
    hold on
%     plot(x1,ones(1,size(x1,2)).*0.5,'^','MarkerSize',2,'color','k');
    %plot(x3,x3.*0-0.4,'.','MarkerSize',7,'color','k');
%     plot(x3,ones(1,size(x3,2)).*-0.5,'o','MarkerSize',2,'color','k');
    
    %set a dashed line to seperate the pvalues from everything else
% %     t7 = [0 60]; t8 = [0.6 0.6];
% %     P=plot(t7,t8);
% %     set(P,'color',[0.5 0.5 0.5],'LineStyle','--','linewidth',1.5);
    % add dashed line to show point of significance  

%             % add dashed line to show point of significance            
%             x=[43,43];
%             y=[0,11.2];
%             P=plot(x,y);
%             set(P,'Color',[cmap(2,:),0.5],'LineStyle','--','linewidth',1);
%             % add dashed line to show point of significance            
%             x=[25,25];
%             y=[0,11.2];
%             P=plot(x,y);
%             set(P,'Color',[cmap(1,:),0.5],'LineStyle','--','linewidth',1);
%         %for windows
%        figname = char(strcat('Z:\ENG_BBCancer_Shared\group\0Zach\Summer 2018\Replicate Selection (5-11)\Output Graphics 6-01 (split)\ProlifPlotAGM',ConditionNameMat(totalcounter),'60118'));
%        xlswrite(char(strcat('Z:\ENG_BBCancer_Shared\group\0Zach\Summer 2018\Replicate Selection (5-11)\Output Graphics 6-01 (split)\FitProlifAGM60118',ConditionNameMat(totalcounter))),storeFitInfo);
%        xlswrite(char(strcat('Z:\ENG_BBCancer_Shared\group\0Zach\Summer 2018\Replicate Selection (5-11)\Output Graphics 6-01 (split)\PValProlifAGM60118',ConditionNameMat(totalcounter))),SignificanceStore);
       %for mac 
       figname = char(strcat('Z:\ENG_BBCancer_Shared\group\0Zach\Summer 2018\Updated Graphics (6-30)\Proliferation Plots (4)\PracticeProlifPlotAGM',ConditionNameMat(totalcounter),'70918'));
%        xlswrite(char(strcat('/Volumes/Research/ENG_BBCancer_Shared/group/0Zach/Summer 2018/Replicate Selection (5-11)/Output Graphics 6-01 (split)\FitProlifAGM60118',ConditionNameMat(totalcounter))),storeFitInfo);
%        xlswrite(char(strcat('Z:\ENG_BBCancer_Shared\group\0Zach\Summer 2018\Updated Graphics (6-30)\Proliferation Plots (2)\PValProlifGFGM63018',ConditionNameMat(totalcounter))),SignificanceStore);

       %        [~,hObj]=legend([gg(1:loopcounter3-1)],{'EGFR Inhibitor','Low EGF','High EGF'},...
%            'EdgeColor','k','Fontsize',8,'location','northwes  t');           % return the handles array
%          hL=findobj(hObj,'type','line');
       
% %          figname=sprintf('%s_%d','SiglineProlifDMSO500_52318',j);
         set(gcf, 'PaperUnits', 'centimeters');
         set(gcf, 'PaperPosition', [0 0 6.2575 10.575]);
         
         
         subplot(2,1,1); set(gca,'XMinorTick','on','YMinorTick','on');
         ax = gca;
         ax.LineWidth = 1;
         subplot(2,1,2); set(gca,'XMinorTick','on','YMinorTick','on');
         ax = gca;
         ax.LineWidth = 1;
%              saveas(gcf,figname,'epsc');
%              print('-depsc','-r1000',figname)
%              print('-dtiff','-r1000',figname)
             
             
    end
clear all; close all; clc;
end
 
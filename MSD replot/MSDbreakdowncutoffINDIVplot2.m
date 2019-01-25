function MSDbreakdowncutoff
%% calculate MSD from only individual cells at early times

% cleanup
close all
clear all
clc

% start timing
tic

% set miscellaneous parameters
FramesToHours = 4;
TimeWindow = FramesToHours; % analyze this frequently (i.e. hours)

PlotInt = 12*FramesToHours;

% bin number of singlets
bins = 0.5:1:21;

%create time array
cmap = cbrewer('qual','Dark2',5);

% minimum number of cells per cluster
cluster_definition=1;

% cutoff distance
threshold = 75;

% set plot parameters
TimeInd = (1:PlotInt:249);
NumTimes = numel(TimeInd);

%% consolidate OHT Assay Figure
AssayOHT = [54];
CellDens = [1000];

plotCount = 1;
figure
% for the selected number of conditions
for condnum = 1:numel(AssayOHT)
    
    % pick out well #
    wellnum = AssayOHT(condnum);
    
    %# generate filename
    tempName = strcat('EGF(E6)w',num2str(wellnum),'.mat');
    
    load(tempName);
    
    % open figure
%     figure('Name',strcat('EGF(E6)w',num2str(wellnum)),'NumberTitle','off')
    
    maxTime = size(storeX,2);
    
    % calculate speeds (again - clean this up later)
    velX=storeX(:,(1+TimeWindow):TimeWindow:end)-...
        storeX(:,1:TimeWindow:end-TimeWindow);
    
    velY=storeY(:,(1+TimeWindow):TimeWindow:end)-...
        storeY(:,1:TimeWindow:end-TimeWindow);
    
    velR = sqrt(velX.*velX + velY.*velY);
    
    %name each single cell into its own cluster (initialize cluster naming)
    storeNeighbors=repmat(1:size(storeX,1),maxTime,1);
    storeNeighbors=storeNeighbors';
    
    %store NaN locations & copy NaN locations to storeNeighbors
    storeNaN=~isnan(storeX);
    storeNeighbors=storeNeighbors.*storeNaN;
    storeNeighbors(storeNeighbors==0)=NaN;
    
    % store cells per neighbor
    storeNeighborNum = zeros(numel(bins)-1,maxTime);
    
    % loop through different distance cutoffts
    for threshold = 75
        
        %calculate and generate storeNeighbors matrix at all times
        for t=1:maxTime
            
            fPresent = find(~isnan(storeX(:,t)));
            
            tempCoords = [storeX(fPresent,t) storeY(fPresent,t)];
            
            [idx,dist] = rangesearch(tempCoords,tempCoords,threshold);
            
            % fill in neighbor numbers
            for i = 1:size(tempCoords,1)
                storeNeighbors(fPresent(i),t) = numel(idx{i})-1;
                
            end
            
        end
    end
    
    % find rows corresponding to individuals at short times
            [row, col] = find(storeNeighbors(:,1:(4*20))==0);
%       [row, col] = find(storeNeighbors(:,1:maxTime) == 0)
    %find(storeNeighbors(:,1:maxTime) == 0)
    
    % pull out unique cell IDs
    IndList = unique(row); %all the rows that are individuals at an early timepoint
    
%     %Kill switch implemenentaiton begins here
%     IndCol = [];
%     
%     %pull out first occurance of single cell behavior:
%     for j = IndList'
%        ix3 = row == j;
%        relCol = col(ix3);
%        sp1 = min(relCol);
%        IndCol = [IndCol; sp1];
%         
%     end
%     
%     Ax3 = [IndList IndCol];
%     %now find where neighbors are found
%     [row1, col1] = find(storeNeighbors(:,1:(10*8))>0);
%         % pull out unique cell IDs
%     IndList1 = unique(row1); %all the rows that are individuals at an early timepoint
%     
%     IndCol1 = [];
%     
%     %pull out first occurance of single cell behavior:
%     for j = IndList1'
%        ix31 = row1 == j;
%        relCol1 = col1(ix31);
%        sp11 = min(relCol1);
%        IndCol1 = [IndCol1; sp11];
%     end
%     
%     Ax4 = [IndList1 IndCol1];
%     uniquCol = unique(Ax4(:,1));
%     ismem = ismember(uniquCol,Ax3(:,1));
%     Ax4 = Ax4(ismem,:); 
%     
%     for q = 1:size(Ax3,1)
%        numCell = Ax3(q,1);
%        det1 = numCell == Ax4(:,1);
%        if sum(det1) == 0
%        else
%            if Ax3(q,2)>Ax4(det1,2)
%               Ax3(q,:) = NaN; 
%            end
%            
%        end
%        
%         
%     end
%     
%     %Kill switch implementaiton ends here
    
    sumcount61 = 0;
    MyList = []; 
    for k6 = 1:numel(IndList)
       jsum = sum(row==IndList(k6)); 
       if jsum>4*5
           MyList = [MyList,IndList(k6)];
           sumcount61 = sumcount61+1;
       end
        
    end
    fprintf('There are %d cells tracked\n',sumcount61)
    IndList = MyList;
    % preallocate cell for MSD analysis
    tracks = cell(length(IndList),1);
    
    % add each cell ID as a unique cell to tracks
    for ID = 1:length(IndList) %for each cellID
        
        % make sure cell is present
        fPresent = find(~isnan(storeX(IndList(ID),:)));
        
        if ~isempty(fPresent)
            % add as time, X, Y
            tracks{ID} = [1/FramesToHours*fPresent'...
                storeX(IndList(ID),fPresent)'...
                storeY(IndList(ID),fPresent)'];
        end
    end
    
    % cleanup and remove empty cells
    tracks = tracks(~cellfun(@isempty,tracks));
    
    StoreMSDdata = cell(numel(tracks),1);
    StoreFitParam = [];
    
    for jj = 1:numel(tracks)
        practicetracks = {tracks{jj}};
        % set msd analyzer to microns and hours
        ma = msdanalyzer(2,'um','h');

        % run msdanalyzer on tracks
        ma = ma.addAll(practicetracks);

        % calculate MSD info
        allMSD = ma.computeMSD;
        calcMSD = allMSD.msd; % pull out MSD info for analysis

        % calculate ensemble averaged MSD
        calcMeanMSD = ma.getMeanMSD(1:length(practicetracks));

        % start sampling at 1 h, due to nuclear deformation
        f = find(calcMeanMSD(:,1) >= 1);

        % pull out MSD values below 10 h
        fearly = find(calcMeanMSD(f,1) <= 10);
%        condnum2=[1 1 2 2 3 3]; 
        
        plot(calcMeanMSD(f(fearly),1),calcMeanMSD(f(fearly),2),'color',cmap(condnum,:))
        hold on
        
            % fit to power law
    fitD = fit(calcMeanMSD(f(fearly),1),...
        calcMeanMSD(f(fearly),2),'power1')
    

    % plot power law for comparison, offset for clarity
    plot(calcMeanMSD(f(fearly),1),...
       1*fitD.a.*calcMeanMSD(f(fearly),1).^fitD.b,'--','color',cmap(condnum,:),...
        'Linewidth',4); hold on
      StoreFitParam = [StoreFitParam; fitD.a fitD.b];
      StoreMSDdata{jj,1} = [calcMeanMSD(f(fearly),1)';calcMeanMSD(f(fearly),2)'];
      
    
    
    
    % show fit parameters
    t=text(11,fitD(10),strcat('\alpha=',num2str(fitD.b)),'color',cmap(condnum,:),'FontSize',10);
        
    end
    
    % set msd analyzer to microns and hours
    ma = msdanalyzer(2,'um','h');
    
    % run msdanalyzer on tracks
    ma = ma.addAll(tracks);
    
    % calculate MSD info
    allMSD = ma.computeMSD;
    calcMSD = allMSD.msd; % pull out MSD info for analysis
    
    % calculate ensemble averaged MSD
    calcMeanMSD = ma.getMeanMSD(1:length(tracks));
    
    % start sampling at 1 h, due to nuclear deformation
    f = find(calcMeanMSD(:,1) >= 1);
    
    
    % pull out MSD values below 10 h
    fearly = find(calcMeanMSD(f,1) <= 10);
      
    % note this is standard error of the mean
%     errorbar(calcMeanMSD(f(fearly),1),calcMeanMSD(f(fearly),2),...
%         calcMeanMSD(f(fearly),3)./sqrt(calcMeanMSD(f(fearly),4)),'o-',...
%         'Color',cmap(plotCount+1,:),...
%         'Linewidth',2); hold on
  
        
    % fit to power law
    fitD = fit(calcMeanMSD(f(fearly),1),...
        calcMeanMSD(f(fearly),2),'power1')
    

    % plot power law for comparison, offset for clarity
    plot(calcMeanMSD(f(fearly),1),...
       1*fitD.a.*calcMeanMSD(f(fearly),1).^fitD.b,'--','color',cmap(condnum,:),...
        'Linewidth',4); hold on
    
    
    
    % show fit parameters
    t=text(11,fitD(10),strcat('\alpha=',num2str(fitD.b)),'color',cmap(condnum,:),'FontSize',11);
   
    % set to log log plot
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
        xlabel('\tau (h)')
    ylabel('MSD (\mum^2)')
    grid off
    xlim([0.8 60])
    hold on
%     
%     figname = 'Z:\ENG_BBCancer_Shared\group\0Zach\Summer 2018\MSD replot\GraphicTest\Testall10cutoff_w54';
%     saveas(gcf,figname,'epsc');
%     print('-dtiff','-r1000',figname)
end

   x=[1,10];
   y1= 0.8*fitD.a.*x;
   y2= 0.8*fitD.a.*x.^2;
   plot(x,y1,'k','Linewidth',4); hold on;
   plot(x,y2,'k','Linewidth',4); hold on;
    
    
   t=text(11,2196,strcat('\alpha=',num2str(1)),'color','k','FontSize',11);
   t=text(11,21960,strcat('\alpha=',num2str(2)),'color','k','FontSize',11);
toc
end

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

% for the selected number of conditions
for condnum = 1:numel(AssayOHT)
    
    % pick out well #
    wellnum = AssayOHT(condnum);
    
    %# generate filename
    tempName = strcat('EGF(E6)w',num2str(wellnum),'.mat');
    
    load(tempName);
    
    % open figure
    figure('Name',strcat('EGF(E6)w',num2str(wellnum)),'NumberTitle','off')
    
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
    [row, col] = find(storeNeighbors(:,1:(10*4)) > -1);
    
    % pull out unique cell IDs
    IndList = unique(row); %all the rows that are individuals at an early timepoint
    
    sumcount61 = 0;
    MyList = []; 
    for k6 = 1:numel(IndList)
       jsum = sum(row==IndList(k6)); 
       if jsum>-1
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
    
%     for jj = 1:numel(tracks)
%         practicetracks = {tracks{jj}};
%         % set msd analyzer to microns and hours
%         ma = msdanalyzer(2,'um','h');
% 
%         % run msdanalyzer on tracks
%         ma = ma.addAll(practicetracks);
% 
%         % calculate MSD info
%         allMSD = ma.computeMSD;
%         calcMSD = allMSD.msd; % pull out MSD info for analysis
% 
%         % calculate ensemble averaged MSD
%         calcMeanMSD = ma.getMeanMSD(1:length(practicetracks));
% 
%         % start sampling at 1 h, due to nuclear deformation
%         f = find(calcMeanMSD(:,1) >= 1);
% 
%         % pull out MSD values below 10 h
%         fearly = find(calcMeanMSD(f,1) <= 10);
%         
%         plot(calcMeanMSD(f(fearly),1),calcMeanMSD(f(fearly),2),'color',cmap(2,:))
%         hold on
%         
%     end
    
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
    errorbar(calcMeanMSD(f(fearly),1),calcMeanMSD(f(fearly),2),...
        calcMeanMSD(f(fearly),3)./sqrt(calcMeanMSD(f(fearly),4)),'o-',...
        'Color',cmap(plotCount+1,:),...
        'Linewidth',2); hold on
  
        
    % fit to power law
    fitD = fit(calcMeanMSD(f(fearly),1),...
        calcMeanMSD(f(fearly),2),'power1')
    
%     %fit to biased diffusion model...
%     F3 = @(z,tdat)z(1)*((z(2)*tdat)-(z(2)^2)*(1-exp(-tdat/z(2))));
%     z0 = [1000,0.001];
%     options=optimset('disp','iter','LargeScale','off','TolFun',.001,'MaxIter',100000000000,'MaxFunEvals',100000000000)
%     %options=optimset('disp','iter','LargeScale','off','TolFun',.001,'MaxIter',100000,'MaxFunEvals',100000);
%     [z,resnorm,~,exitflag,output] = nlinfit(calcMeanMSD(f(fearly),1),...
%         calcMeanMSD(f(fearly),2),F3,z0);
%     disp(exitflag)
%     disp(output)
%     z

%     problem = createOptimProblem('lsqcurvefit','x0',z0,'objective',F3,...
%             'xdata',calcMeanMSD(f(fearly),1),'ydata',calcMeanMSD(f(fearly),2));
%         
%         ms = MultiStart('PlotFcns',@gsplotbestf);
%         [xmulti,errormulti] = run(ms,problem,2000);
%         meanres = sum((mean(calcMeanMSD(f(fearly),2))-calcMeanMSD(f(fearly),2)).^2);
%         rsquare = 1-(errormulti/meanres);
%         x = xmulti;
    
    % plot power law for comparison, offset for clarity
    plot(calcMeanMSD(f(fearly),1),...
       2*fitD.a.*calcMeanMSD(f(fearly),1).^fitD.b,'k-',...
        'Linewidth',2); hold on
    
    
    % show fit parameters
    t=text(1,3000,strcat('\alpha=',num2str(fitD.b)));
    
    % set to log log plot
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
        xlabel('\tau (h)')
    ylabel('MSD (\mum^2)')
    grid on
    xlim([0.8 20])
%     
    figname = 'Z:\ENG_BBCancer_Shared\group\0Zach\Summer 2018\MSD replot\GraphicTest\Testall10cutoff_w54';
    saveas(gcf,figname,'epsc');
    print('-dtiff','-r1000',figname)
end

toc
end

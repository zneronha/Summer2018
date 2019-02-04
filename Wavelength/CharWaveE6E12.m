function CharWaveE6E12
%31 May 2018 <Zachary_Neronha@brown.edu> 
%this code creates a scatter plot of clusters and number of leaders
%Updated on 6 June for some improvements and curve fitting
%Updated on 14 June for log log plotting equal 
%Updated on 2 July to comptue the radius of Gyration too

%to change between E6 and E12 change: wells, load path, ws modifier

%cleanup
clearvars; close all; clc

%define replicates horizontally and conditions vertically
%wells = [57 58; 9 10];
wells = [6 19 1257 1258; 54 67 129 1210];

%allocate data storage
storedata4 = cell(2,2);
%create a color matrix
cpal = cbrewer('qual','Set1',4,'cubic');
colorcounter = 1; %initiate a color counter

storePTS2 = cell(size(wells,1),size(wells,2));

%loop through the conditions
for u2 = 1:size(wells,1)
    
    %loop through the replicates
    for u3 = 1:size(wells,2)
        storePTS = nan(400,4);
        indexer7 = 1;
        %extract the well of interest
        well = wells(u2,u3);
        storedata3 = []; %#ok<NASGU>
        
        if length(num2str(well)) > 2
        %load the leader data
        load('Z:\ENG_BBCancer_Shared\group\0Zach\Leader Data\TriLeaderData31418\DataProcessing060418\EGF(E12)_Thr1.compiledleaderdata060418.mat',...
            'CompiledDataStore');
            ww = num2str(well);
            wellc = ww(3:end);
            wellc = str2double(wellc);
        else
        load('Z:\ENG_BBCancer_Shared\group\0Zach\Leader Data\TriLeaderData31418\DataProcessing031518\EGF(E6)_Thr1.compiledleaderdata031518.mat',...
            'CompiledDataStore');
            wellc = well;
        end
        
        %get to the appropriate well's leader data
        leaderdata = CompiledDataStore{wellc}; %#ok<IDISVAR,USENS>
        
        %extract the number of frames and convert to hours
        framenum = size(leaderdata,2);
        
        %load the appropriate well
        %ws = strcat('12',num2str(well));
        ws = strcat('',num2str(well));
        
        %load the position data for this well
        load(strcat('Z:\ENG_BBCancer_Shared\group\0Zach\EGFReplicateData\EGF(E6)w',...
            ws,'.mat'),'storeX','storeY');

        %at each frame compile a list of all unique leaders
        for timept = 1:framenum
            %make a list of the unique leaders
            leaderlist2 = cell(3,1);
            
            %cycle through each user
            for user = [1 3 5]
                %extract and store their leaders
                leaderlist1 = leaderdata{user,timept};
                leaderlist2{(user + 1)/2,1} = leaderlist1;
            end
            %compile those leaders
            leaderlist2 = [leaderlist2{1,1};leaderlist2{2,1};leaderlist2{3,1}];
            
            %convert time point to frame 
            framenum = (timept * 5 * 4)-3;
            
            %now pull out the unique leaders and their positions
            leaderlist3 = unique(leaderlist2,'rows');
            %cleave off the position data
            leaderPOScheck = leaderlist3(:,[2 3]);
            leaderlist3 = leaderlist3(:,1);
            %define our positions from the storeX and storeY matricies
            txx = storeX(leaderlist3,framenum);
            tyy = storeY(leaderlist3,framenum);
            
            %check to make sure book keeping is working correctly here
            if sum(leaderPOScheck(:,1)-txx)==0 && sum(leaderPOScheck(:,2)-tyy)==0
            else
                error('Leader positions are not compatable. Check data sources')
            end

            %compute what the clustered cells are using our webplot fun
            [~,~,~,clusterlist] = WebPlotProxSearch2(storeX,storeY,framenum,2);
            clusterlist = [clusterlist zeros(size(clusterlist,1),1)]; %#ok<AGROW> false error message
            
            %loop through the leaders we identified
            for jj = leaderlist3'
                %extract that particular leader
                vv = leaderlist3 == jj;
                
                %see if it matches any of our clustered cells 
                u = clusterlist(:,4)==leaderlist3(vv);
                if sum(u) == 1 %if it does
                    %first verify that the position data lines up with what
                    %weve recorded elsewhere
                    
                    %determine just what cell that is in our clustered
                    %matrix and record it!
                    vv = clusterlist(:,4) == jj;
                    clusterlist(vv,5) = 1;
                else
                    %tie up loose ends, display a warning message if needed
                    if sum(u) == 0
                    else
                        error('Crash: leader appears more than once!')
                    end
                end

            end
            
            %store and extract the unique clusters 
            clTRUN = clusterlist(:,[1,5]);
            clTRUN2 = unique(clTRUN(:,1));
            clTRUN2 = [clTRUN2 nan(size(clTRUN2,1),1)];  %#ok<AGROW> another false error
            
            %compute radius of gyration for each cluster
            
            
            %for every cluster
            for jj = 1:size(clTRUN2,1)
                %which cells are in that cluster
                uu = clTRUN(:,1) == clTRUN2(jj,1);
                %compute and record how many of them are leaders
                clTRUN2(jj,2) =  sum(clTRUN(uu,2));
            end
            
            clTRUN2 = [clTRUN2 nan(size(clTRUN2,1),1)]; %#ok<AGROW>
            %now assign the clustercell count to the matricies
            counter3 = 1;
            for yz = clTRUN2(:,1)'
                ix3 = clTRUN(:,1) == yz;
                clTRUN2(counter3,3) = sum(ix3);
                counter3 = counter3+1; 
            end
            
            %compute radius of gyration for each cluster
            clTRUN2 = [clTRUN2 nan(size(clTRUN2,1),1)]; %#ok<AGROW>
            mcc = 1;
            for rad3 = clTRUN2(:,1)'
                ww = clusterlist(:,1) == rad3;
                ptPOS = clusterlist(ww,2:3);
                [clTRUN2(mcc,4)] = ROG(ptPOS);
               
                mcc = mcc+1;
            end
            
            

            storePTS((indexer7:indexer7+(size(clTRUN2,1))-1),:) = clTRUN2;
            indexer7 = indexer7+size(clTRUN2,1);

        end
            storePTS2{u2,u3} = storePTS;
           colorcounter = colorcounter + 1;
           

    end
end
%combine data and plot
AssayDMSO = [storePTS2{1,1};storePTS2{1,2};storePTS2{1,3};storePTS2{1,4}];
AssayOHT = [storePTS2{2,1};storePTS2{2,2};storePTS2{2,3};storePTS2{2,4}];

%get rid of clusters without leaders
AssayDMSO = AssayDMSO(AssayDMSO(:,2)>0,:);
AssayOHT = AssayOHT(AssayOHT(:,2)>0,:);


% AssayOHT(:,4) = AssayOHT(:,4)./50;
% 
% %make sure offsetting is consistant across log axes 
% AssayDMSO(:,2) = AssayDMSO(:,2)+


% plot(AssayDMSO(:,3),AssayDMSO(:,2)+0.3,'.','color',cpal(2,:),'MarkerSize',25,'linewidth',1.3)
% plot(AssayOHT(:,3),AssayOHT(:,2),'.','color',cpal(3,:),'MarkerSize',25,'linewidth',.85)
% f2 = scatter(AssayDMSO(:,3),AssayDMSO(:,2)+0.05.*AssayDMSO(:,2),'MarkerFaceColor',cpal(2,:),'MarkerEdgeColor',[1 1 1]);
% f2.SizeData = 60;
% f2.MarkerFaceAlpha = 0.8;

% f3 = scatter(AssayOHT(:,3),AssayOHT(:,2)-0.05.*AssayOHT(:,2),'MarkerFaceColor',cpal(3,:),'MarkerEdgeColor',[1 1 1]);
% f3.SizeData = 60;
% f3.MarkerFaceAlpha = 0.8;

f3 = scatter(AssayOHT(:,4),AssayOHT(:,2),'MarkerFaceColor',cpal(2,:),'MarkerEdgeColor',[1 1 1]);
f3.SizeData = 60;
% f3.MarkerFaceAlpha = 0.8;

xa = AssayOHT(AssayOHT(:,2)>0,4); ya = AssayOHT(AssayOHT(:,2)>0,2);

% y = y(1:13); x = x(1:13);
st = [1 1];
L1 = [0 1];
U1 = [inf 1];
fo1 = fitoptions('method','NonlinearLeastSquares','Lower',L1,'Upper',U1);
set(fo1,'Startpoint',st);
ft1 = fittype('A*(x^B)',...
     'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'A', 'B'});
 [cf,G] = fit(xa,ya,ft1,fo1)
 hello = plot(cf);
 hello.LineWidth = 2;
 hello.Color = [0 0 0];

    
set(gca, 'XScale', 'lin')
set(gca, 'YScale', 'lin')
%ylim([0.8 20]); 
%xlim([2 500]);
legend({'Assay OHT'},'location','nw')
xlabel('Cluster Radius of Gyration')
ylabel('Number of Leaders Identified')
% legend('Assay DMSO','Assay DMSO','Assay OHT','Assay OHT','location','nw')
% xlabel('Time (hours)')
% ylabel('Maximum Leaders in a Cluster')


ax=gca;
box on
grid off
ax.XColor='black';
ax.YColor='black';
set(gca,'fontsize',8);
set(gca,'XMinorTick','off','YMinorTick','off')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 10 7.5])

title('Fractal Dimension: OHT low EGF','FontSize',10);
figname=strcat('C:\Users\zjner\GitHub\Summer2018\Wavelength\ClusterSizeNumLeader122519');
saveas(gcf,strcat(figname,'.pdf'))
print('-dtiff','-r1000',figname)
end

function [RadGyr] = ROG(Data)

%remove nan
di = ~isnan(Data(:,1));
Data = Data(di,:);

Xcom = mean(Data(:,1)); Ycom = mean(Data(:,2));
MS = [(Data(:,1)-Xcom).^2 (Data(:,2)-Ycom).^2];
MSS = MS(:,1) + MS(:,2);
SumVect = sum(MSS);
xx = SumVect./size(Data(:,1),1);
RadGyr = sqrt(xx);

end

function [clustercellcount, singlecellcount, numCLU,clusterlist] =...
    WebPlotProxSearch2(storeX,storeY,frame,cluster_definition)
%This code is intended to be a training code to determine the the ideal
%buffer region to quantify clusters.

%24 APRIL 2018: NEW TO THIS VERSION IS CELL COUNT CLUSTER EXCLUSION AND
%FINAL CLUSTER COUNTS

% clearvars; close all; clc

% frame = 5;
range = 75;
cluthresh = cluster_definition;



tempX = storeX(:,frame);
tempY = storeY(:,frame);
cluID = [1:size(tempX,1)]';
lx = isnan(tempX) == 0;
ly = isnan(tempY) == 0;
tempX = tempX(lx);
tempY = tempY(ly);
cluID = cluID(lx);
hold on

[ix] = rangesearch([tempX tempY],[tempX tempY],range);

%make a list which states which cluster a cell is in
clusterlist = nan(numel(tempX),1);
clusterlist = [clusterlist tempX tempY cluID];
% %add a redundant indexing column
% zz = 1:numel(tempX);
% clusterlist = [zz' clusterlist ];
%initiate the cluster counter
clustercounter = 1;

%loop through all the cells
for ii = 1:numel(tempX)
    %determine which cells have been assigned to a clsuter
    uu = isnan(clusterlist);
    
    %if none of the cells this cell is in a cluster with are clustered
    %already
    if sum(uu(ix{ii}))== numel(tempX(ix{ii}))
        %assign all of them to this cluster
        clusterlist(ix{ii},1) = clustercounter;
        %and add one to the cluster counter
        clustercounter = clustercounter + 1;
    else
        %otherwise we must determine which cells are clustered already
        %first locate the ID's of those cells already in clusters
        q3 = find(uu(:,1)==0);
        %then locate those cells which are in our current cluster
        q4 = ix{ii};
        %find all cells that are members of both groups
        PCC = q4(ismember(q4',q3));
        %now determine the id's of those clusters
        pccIDs = unique(clusterlist(PCC,1));
        %take the lowest number
        minID = min(pccIDs);
        
        %now label everbody with those ids
        clusterlist(ix{ii},1) = minID;
        for uv = pccIDs'
           v3 = find(clusterlist(:,1)==uv);
           clusterlist(v3,1)=minID;
        end
        
    end
    
end
    
    %calculate total number of cells
    totalcell = sum(isnan(tempX)==0);
    %calculate number of single cells
    clustercellcount = 0;
    singlecellcount = 0;
    
    %first get rid of any cluster with too few cells
    for oo = unique(clusterlist(:,1))'
       clusize6 = sum(clusterlist(:,1)==oo);
       if clusize6<cluthresh
           killmarker = clusterlist(:,1) == oo;
           clusterlist(killmarker,:) = NaN;
           tempX(killmarker) = NaN;
           tempY(killmarker) = NaN;
           singlecellcount = singlecellcount + clusize6;
       else
           clustercellcount = clustercellcount + clusize6;
       end

    end
    
    %consolidate positions
    ix8 = isnan(tempX)==0;
    clusterlist = clusterlist(ix8,:);
%     tempX = tempX(ix8);
%     tempY = tempY(ix8);


numCLU = numel(unique(clusterlist));
end

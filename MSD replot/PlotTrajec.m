function PlotTrajec
%Purpose: plot trajectories for cells to verify "randomness"

cmap = cbrewer('qual','Set1',8);
AssayOHT = [51 54 57 60];
wellnum = 54;

%load the data
tempName = strcat('EGF(E6)w',num2str(wellnum),'.mat');
load(tempName);

for ii = 1:size(storeX,1)
    
    %find first recorded position 
    for jj = 1:size(storeX,2)
       qq = isnan(storeX(ii,jj));
       if qq == 0
          testpos = jj; 
       end
        
    end
    
    
    gx = storeX(ii,:)-storeX(ii,jj);
    gy = storeY(ii,:)-storeY(ii,jj);
    
    plot(gx,gy)
    hold on
    
end




end
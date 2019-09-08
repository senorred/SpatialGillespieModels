%% Bursting video

c=jet(400);

fig = figure;
% if sum(CrowdSpace(:))>0
%     [crow,ccol,v]=find(CrowdSpace);
%     scatter(crow,ccol,250,'s','c')
% end
% hold on
% The final plot.
 
%% Set up the movie.
writerObj = VideoWriter('out.avi'); % Name it.
writerObj.FrameRate = 2; % How many frames per second.
open(writerObj); 

for i=1:tmax      
    
    pause(0.1);
    figure(fig); % Makes sure you use your desired frame.
    
    mRNASpace=mRNASpaceTrack(:,:,:,i);
    RibosomeSpace=RiboSpaceTrack(:,:,:,i);
    
    [row,col,dep]=ind2sub([VoxLength VoxWidth VoxHeight],find(mRNASpace));
    vm=mRNASpace(mRNASpace~=0);
    [rowR,colR,depR]=ind2sub([VoxLength VoxWidth VoxHeight],find(RibosomeSpace));
    vR=RibosomeSpace(RibosomeSpace~=0);

    if sum(CrowdSpace(:))>0
    [crow,ccol,cdep]=ind2sub([VoxLength VoxWidth VoxHeight],find(CrowdSpace));
    scatter3(crow,ccol,cdep,250,'*','c') %'MarkerFaceColor','c')
    hold on
    end
    scatter3(row,col,dep,vm*10,'x','k')
    hold on
    scatter3(rowR,colR,depR,vR*10)
    hold off
    
    ax=gca;
%     ax.XLim=[0,10];
%     ax.YLim=[0,10];
%     ax.ZLim=[0,10];
    
    %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
        
    %end
 
end
%hold off
close(writerObj); % Saves the movie.
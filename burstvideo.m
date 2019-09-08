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
writerObj.FrameRate = 5; % How many frames per second.
open(writerObj); 

for i=1:tmax      
    
    pause(0.1);
    figure(fig); % Makes sure you use your desired frame.
    
    indexM= find(mRNATrack(:,1,i));
    m=mRNATrack(indexM,:,i);
    mRNASpace=zeros(20,20);
    for j=1:length(indexM)
        mRNASpace(m(j,1),m(j,2))=mRNASpace(m(j,1),m(j,2))+1;
    end
    
    indexR= find(RibosomeTrack(:,1,i));
    r=RibosomeTrack(indexR,:,i);
    RibosomeSpace=zeros(20,20);
    for j=1:length(indexR)
        RibosomeSpace(r(j,1),r(j,2))=RibosomeSpace(r(j,1),r(j,2))+1;
    end
    
    [row,col,v]=find(mRNASpace);
    [rowR,colR,vR]=find(RibosomeSpace);
    %[rowB,colB,vB]=find(RibosomeSpace & mRNASpace);

%     mRNANum =length(row);
%     destNum =length(rowR);
%     bNum =length(rowB);    
    %hold on
    if sum(CrowdSpace(:))>0
    [crow,ccol]=find(CrowdSpace);
    scatter(crow,ccol,250,'s','MarkerFaceColor','c')
    hold on
    end
    scatter(row,col,v*100,'x','k')
    hold on
    scatter(rowR,colR,vR*10)
    hold off
    
    ax=gca;
    ax.XLim=[0,20];
    ax.YLim=[0,20];
    
    
    %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
        
    %end
 
end
%hold off
close(writerObj); % Saves the movie.
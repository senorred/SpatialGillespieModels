%% Hippoplotamus SXPM 2.3
% Judges winners by occupancy over time

tInt=100;

mRNASpace=mRNASpaceTrack(:,:,tInt);
RibosomeSpace=RiboSpaceTrack(:,:,tInt);

ccc=prism(10000);
cc=ones(10000,3);
k=find(RibosomeSpace & mRNASpace);
c=jet(length(k));
cb=brighten(jet(length(k)),0.8);
cd=brighten(jet(length(k)),-0.6);

% make every 2nd value light and every 3rd value dark
for n=1:length(k)
    if rem(n,2)==1
        c(n,:)=cb(n,:);
%     elseif rem(10*n,2)==2
%         c(n,:)=cb(n,:);
    else
        c(n,:)=cd(n,:);
    end
end

% A=randperm(length(k));
% %Randomize the hsv colormap
% for n=1:length(k)
%     c(n,:)=c(A(n),:);
% end

for n=1:length(k)
    cc(k(n),:)=c(n,:);
end

i=1;
fNumStore=[];
figure('Position',[0 0 800 800])
if sum(CrowdSpace(:))>0
    [crow,ccol]=find(CrowdSpace);
    scatter(crow,ccol,800,'s','MarkerFaceColor','k','MarkerEdgeColor','k')
end
hold on

indWinnerStore=[];


for x=1:VoxLength
    for y=1:VoxWidth
        
        if CrowdSpace(x,y)<1
        % Find all the ribosomes that existed at some point in time in the
        % chosen voxel
        [riboNum, tOcc]= find(squeeze(RibosomeTrack(:,1,1:tInt))==x & squeeze(RibosomeTrack(:,2,1:tInt))==y);
            if isempty(riboNum)
                scatter(x,y,800,'s','MarkerEdgeColor','k','MarkerFaceColor','w')
                continue
            end
        riboNum=unique(riboNum);
        finalDest=zeros(length(riboNum),2);
        
        % Find the final destination of all found Ribos at end time
        for r=1:length(riboNum)
            finalDest(r,:) = squeeze(RibosomeTrack(riboNum(r),:,tInt));
        end
        
        % What is the most common final destination?
        indDest= sub2ind([VoxLength,VoxWidth],finalDest(:,1),finalDest(:,2));
        indWinner= mode(indDest);
        
        % Color the chosen voxel the winning destination's index color
        scatter(x,y,1600,'s','MarkerEdgeColor','k','MarkerFaceColor',cc(indWinner,:))
        indWinnerStore=[indWinnerStore;indWinner];
        end
    end
end

indWinners=unique(indWinnerStore);
watershedNum=length(indWinners);

areaWinners= zeros(watershedNum,1);
for var=1:watershedNum
    A=indWinners(var);
    areaWinners(var)= length(indWinnerStore(indWinnerStore==A));
end

[winX,winY]=ind2sub([VoxLength,VoxWidth],indWinners);

%scatter(winX,winY,150,'*','MarkerEdgeColor','k')


hold off
    

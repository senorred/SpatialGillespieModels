c=jet(400);
% cc=[];
% for n=1:bNum
% k=find(rowR== rowB(n) & colR== colB(n));
% cc(n,:)= c(k,:);
% end

%fig = figure;
%hold on
% The final plot.
figure
for i=1:tmax      
    %figure(fig); % Makes sure you use your desired frame.
    
    indexM= find(mRNATrack(:,1,i));
    m=mRNATrack(indexM,:,1);
    mRNASpace=zeros(20,20);
    for j=1:length(indexM)
        mRNASpace(m(j,:))=mRNASpace(m(j,:))+1;
    end
    
    indexR= find(RibosomeTrack(:,1,i));
    r=RibosomeTrack(indexR,:,1);
    RibosomeSpace=zeros(20,20);
    for j=1:length(indexR)
        RibosomeSpace(r(j,:))=RibosomeSpace(r(j,:))+1;
    end
    
    [row,col,v]=find(mRNASpace);
    [rowR,colR,vR]=find(RibosomeSpace);
    [rowB,colB,vB]=find(RibosomeSpace & mRNASpace);

    mRNANum =length(row);
    destNum =length(rowR);
    bNum =length(rowB);
    
    
    %hold on
    scatter(row,col,v*10,'x','k')
    scatter(rowR,colR,vR*10,c)
    %drawnow
    %hold off
    
 
end

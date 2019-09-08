%% Master Plot for SXPM 2.0 trials
% Description: Wrapper for plotting multiple simulation runs
% 1000 second runs, 50 iterations per condition
% 20x20 flat XY area
% 400 Ribosomes, seeded randomly
% alpha= kCM*0.001
% dRibo= 1.0
% dPoly= 0.0001
% Reflective Boundaries

%% MODEL -- SXPM 2.0

alphaConditions= [0.001];
crowdConditions= [0, 10, 20, 30, 40, 50];
iterations= 100;
tInt=50
voccIntRStore=zeros(tInt,iterations,length(crowdConditions));
figure %% for all the traces
hold on

for I=1:length(alphaConditions)
    II= alphaConditions(I)
    for J= 1:length(crowdConditions)
        JJ= crowdConditions(J)
%         figure  % for each CF
%         hold on
        for K= 1:iterations
            KK= K;
            
            c=hsv(length(crowdConditions));
            
            
%% clear variables, seed, and name the run
clearvars -except alphaConditions crowdConditions iterations II JJ KK I J K c tInt voccIntRStore

Name = sprintf('alpha%gCrowd%giter%g.mat',II,JJ,KK);
load(Name)
voccR=[];
vsumR=[];
for i=1:tInt

    mRNASpace=mRNASpaceTrack(:,:,i);
    RibosomeSpace=RiboSpaceTrack(:,:,i);
    
    [row,col,v]=find(mRNASpace);
    [rowR,colR,vR]=find(RibosomeSpace);
    [rowB,colB]=find(RibosomeSpace & mRNASpace);
    mRNANum =length(row);
    destNum =length(rowR);
    bNum =length(rowB); % Number of Polysomes, burst frequency

    %c=jet(bNum);
    occupiedribos=zeros(bNum,1);
    for n=1:bNum
        k=find(rowR== rowB(n) & colR== colB(n));
        occupiedribos(n)= vR(k);
    end
    
    sumOR= sum(occupiedribos);
    freeR= Ribosomes-sumOR;
    
    voccR(i)= sumOR;
    vfreeR(i)= freeR;
end

voccIntR= cumsum(voccR); % integral of ribo occupancy (analogous to protein)
voccIntRStore(:,K,J)=voccIntR;


plot(tspan(1:tInt), voccIntR(1:tInt), 'color', c(J,:))
title('Occupied Ribosomes over time','FontSize',15)
set(gca,'fontsize',10);
%set(gca,'YScale','log');
%set(gca,'XScale','log');
xlabel('Time','FontSize',15)
%axis([200 1800 0 0.01])
%ax = gca;
%ax.XTick = 2:2:70;
%grid on
ylabel('Occupied Ribosomes','FontSize',15)


        end
    end
%     hold off % for each CF
end
hold off %% for all the traces

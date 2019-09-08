%% Description 

tInt=200;

alphaConditions= [0.001];
crowdConditions= [0, 10, 20, 30, 40, 50];
iterations= 100;

burstSizeSum=zeros(200,6);

for I=1:length(alphaConditions)
    II= alphaConditions(I)
    for J= 1:length(crowdConditions)
        JJ= crowdConditions(J)
        for K= 1:iterations
            KK= K;
            
%% clear variables, seed, and name the run
clearvars -except alphaConditions crowdConditions iterations II JJ KK I J K c tInt burstSizeSum

Name = sprintf('alpha%gCrowd%giter%g.mat',II,JJ,KK);
load(Name)

%% Script


PolyIntSpace=zeros(VoxLength,VoxWidth);
PolySpaceTrack=RiboSpaceTrack;
PolySpaceTrack(~(PolySpaceTrack&mRNASpaceTrack))=0; % PolySpaceTrack: Ribosome intensity in locations where mRNA and Ribosomes exist over time

i=1;
for i=1:tInt %integrated over this many seconds
    PolySpace=PolySpaceTrack(:,:,i);
    PolyIntSpace= PolyIntSpace+PolySpace;   
end

nummRNA=length(find(mRNATrack(:,1,tInt)));
birthTime=NaN(nummRNA,1);
burstSize=NaN(nummRNA,1);
proteinPower=NaN(nummRNA,1);

for j=2:nummRNA+1
    birthTime(j-1)= find(mRNATrack(j,1,1:tInt),1,'first');  % Time when mRNA appeared
    burstSize(j-1)= PolySpaceTrack(mRNATrack(j,1,tInt),mRNATrack(j,2,tInt),tInt); % Number of ribos at end of time
    proteinPower(j-1)= PolyIntSpace(mRNATrack(j,1,tInt),mRNATrack(j,2,tInt)); % Protein production estimate over time
end

rankMat=[birthTime, burstSize, proteinPower];

burstSizeSum(1:nummRNA,J)=burstSizeSum(1:nummRNA,J)+burstSize;
%hist(F)
%title('Histogram of Ribosome Occupancy per voxel over time, per Crowding Fraction')


        end
    end
end

burstSizeMean=burstSizeSum/iterations;
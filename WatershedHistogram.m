%% Analysis wrapper SXPM 2.3
% Description: Find the number and areas of the watersheds
% 100 second runs, 50 iterations per condition
% Var x Var flat XY area
% Var x Var Ribosomes, seeded randomly
% alpha= 0.001
% dRibo= 1.0
% dPoly= 0.0001
% Reflective Boundaries

%% MODEL -- SXPM 2.3

areaConditions= [10 20 30 40 50];
iterations= 50;
areaWinnersStore=nan(150,50,5); 
tInt=100;

for I=1:length(areaConditions)
    II= areaConditions(I)
    for K= 1:iterations
        KK= K

%% clear variables and name
clearvars -except areaConditions iterations II KK I K tInt areaWinnersStore

Name = sprintf('area%gx%giter%g.mat',II,II,KK);
load(Name)

indWinnerStore=[];


for x=1:VoxLength
    for y=1:VoxWidth       
        if CrowdSpace(x,y)<1
        % Find all the ribosomes that existed at some point in time in the
        % chosen voxel
        [riboNum,~]= find(squeeze(RibosomeTrack(:,1,1:tInt))==x & squeeze(RibosomeTrack(:,2,1:tInt))==y);
            if isempty(riboNum)
                %scatter(x,y,800,'s','MarkerEdgeColor','k','MarkerFaceColor','w')
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
        %scatter(x,y,800,'s','MarkerEdgeColor','k','MarkerFaceColor',cc(indWinner,:))
        % Store the index of the winning destination
        indWinnerStore=[indWinnerStore;indWinner];
        end
    end
end

% Find how many winners there are, this is how many watersheds there are
indWinners=unique(indWinnerStore);
watershedNum=length(indWinners);

% Count the number of times the winner appears in storage, this is the area
% of the corresponding watershed
areaWinners= zeros(watershedNum,1);
for var=1:watershedNum
    A=indWinners(var);
    areaWinners(var)= length(indWinnerStore(indWinnerStore==A));
end

[winX,winY]=ind2sub([VoxLength,VoxWidth],indWinners);
            
areaWinnersStore(1:watershedNum,K,I)= areaWinners;
%%Wrapper
        
    end
end

areaHist=reshape(areaWinnersStore, [7500,5]);
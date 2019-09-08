%% Parameter sweep wrapper SXPM 3.0
% Description: Trial wrapper for making multiple simulation runs
% 200 second runs, 50 iterations per condition
% 10x10x10 volume
% 1000 Ribosomes, seeded randomly
% alpha=0.001
% dRibo= 1.0
% dPoly= 0.0001
% Reflective Boundaries

%% MODEL -- SXPM2_0

alphaConditions= [0.001];
crowdConditions= [0, 30, 50];
iterations= 50;

for I=1:length(alphaConditions)
    II= alphaConditions(I)
    for J= 1:length(crowdConditions)
        JJ= crowdConditions(J)
        for K= 1:iterations
            KK= K
%% PASTE MODEL HERE
%%
%%
%% clear variables, seed, and name the run
clearvars -except alphaConditions crowdConditions iterations II JJ KK I J K alpha
rand('state',sum(100*clock)); %#ok<RAND>

Name = sprintf('alpha%gCrowd%giter%g.mat',II,JJ,KK);
disp(Name)

%% Specify Parameters
% runtime
tmax= 200;  %100
dt = 1; %step to record state
tspan = 0:dt:tmax;

% grid size
VoxLength=10;
VoxWidth=10;
VoxHeight=10;

% Ribo/Polysome diffusion
dRibo = 1.0;    % 1.0
dPoly = 0.0001; % 0.0001
% Translation rate
kP=1.0;

% Species
Ribosomes= 1000;
InitialmRNA=0;
mRNACount= 0;
mRNATrack2=[];

% Crowders
CrowdVol= JJ;
NumCrowders = floor(VoxLength*VoxWidth*VoxHeight*(CrowdVol/100));
rReject=0;

% mRNA birth rate
%kCM= 1/((100-CrowdVol)/100);
alpha= II;

%% Arrays for storage
% population of ribos per voxel
% coords of mRNAs, ribos
    mRNASpace = zeros(VoxLength,VoxWidth,VoxHeight);
    RibosomeSpace = zeros(VoxLength,VoxWidth,VoxHeight);
    CrowdSpace= zeros(VoxLength,VoxWidth,VoxHeight);

%% Initialize Model

    %track spatial species on an individual basis
    %rows,for each species, xyz location, for all time.
    RibosomeTrack = zeros(Ribosomes,3,length(tspan));
    mRNATrack = zeros(1,3,length(tspan));
    %CrowdTrack= zeros(NumCrowders, 3,length(tspan));
    
    mRNASpaceTrack=zeros(VoxLength, VoxWidth, VoxHeight, length(tspan));
    RiboSpaceTrack=zeros(VoxLength, VoxWidth, VoxHeight, length(tspan));
    
            %Randomly distribute Crowders into compartments, only one
            %crowder per compartment
    i=1;        
    while i <= NumCrowders
        comp = ceil([rand*VoxLength, rand*VoxWidth, rand*VoxHeight]);
        if CrowdSpace(comp(1),comp(2),comp(3))== 0
            CrowdSpace(comp(1),comp(2),comp(3)) = CrowdSpace(comp(1),comp(2),comp(3)) + 1;
            %CrowdTrack(i,:,1) = comp;
            i=i+1;
        else
        end
    end
    
        %Randomly distribute Ribos into compartments NOT OCCUPIED BY
        %CROWDERS
    i=1;    
    while i <= Ribosomes
        comp = ceil([rand*VoxLength, rand*VoxWidth, rand*VoxHeight]);
        if CrowdSpace(comp(1),comp(2),comp(3))==0
            RibosomeSpace(comp(1),comp(2),comp(3)) = RibosomeSpace(comp(1),comp(2),comp(3)) + 1;
            RibosomeTrack(i,:,1) = comp;
            i=i+1;
        else 
        end
    end

RiboTrackInit=RibosomeTrack(:,:,1);
%CrowdTrackInit=CrowdTrack(:,:,1);

%SpatialRxnMatrix [mRNA X  Y  Z]
SpatialRxnMatrix = [ 1  0  0  0; %mRNA born
                     0  1  0  0; %move ribo up 1 in x
                     0 -1  0  0; %move ribo down 1 in x
                     0  0  1  0; %move ribo up 1 in y
                     0  0 -1  0; %move ribo down 1 in y
                     0  0  0  1; %move ribo up 1 in z
                     0  0  0 -1]; %move ribo down 1 in z

%Track propensities across each compartment space
aSpace = zeros(3,1);

mRNACurrent = mRNATrack(:,:,1);
%mRNACurrent=[];
RibosomeCurrent = RibosomeTrack(:,:,1);


TTrack = zeros(length(tspan), 1); %Time tracking
T = 0; %Start Time
RxnCount = 1;  %Reaction counter
tempBurstCount = 0;
RecordTime = dt; %Recording time
RecordCount = 1;

a = zeros(VoxLength,VoxWidth,VoxHeight,7);
amn = zeros(VoxLength,VoxWidth,VoxHeight);

%% Model
    
    Timer = 0;
    %disp('Gillespie Start')
    while T < tmax %Run gillespie until time is up

%         %Timer for console display
%         if floor(T) > Timer
%             Timer = Timer + 1;
%             disp(Timer)
%         end

        %calculate Propensities
        for x = 1:VoxLength
            for y = 1:VoxWidth
                for z = 1:VoxHeight
                
                    a(x,y,z,1)= alpha;
                    
                    if mRNASpace(x,y,z) >= 1
                        a(x,y,z,2:7) = dPoly*RibosomeSpace(x,y,z);
                    else
                        a(x,y,z,2:7) = dRibo*RibosomeSpace(x,y,z);
                    end
                    
                    % Since crowders don't move, I assume this is okay...
                    if CrowdSpace(x,y,z)==1
                        a(x,y,z,:)=0;
                    end
                    
                end
            end
        end
        a0 = sum(a(:));

        rtime = rand();
        ra = rand();

        tau = (1/(a0))*log(1/ rtime);  %time-step increment
        T = T + tau;   %increment current time by tau
        flag = 0;
        
        % Choose a reaction, get the X and Y and Z coordinates and the reaction
        % type, 'q'
        q=[];
        avec=reshape(a,[1,(VoxLength*VoxWidth*VoxHeight*7)]);
        ia = find((cumsum(avec) >= ra*a0),1,'first');
        [xVox, yVox, zVox,q]=ind2sub(size(a),ia);
        
        % "Does a crowder exist here" gut check-- shouldn't 
        if CrowdSpace(xVox,yVox,zVox)>0
            disp('whoa partner, your code sucks!')
        else
        % According to which reaction type is chosen, execute reaction
        if q == 1   %makes mRNA in chosen voxel
            mRNASpace(xVox,yVox,zVox) = mRNASpace(xVox,yVox,zVox) + 1;
            mRNACount= mRNACount+1;
            %disp('Made mRNA!')
            mRNACurrent=vertcat(mRNACurrent, [xVox, yVox, zVox]);
            mRNATrack= vertcat(mRNATrack, zeros(1,3,length(tspan)));

        elseif q>=2  %moves Ribosome from chosen voxel
            % calculate destination coordinates
            dest= [xVox,yVox,zVox]+ SpatialRxnMatrix(q,2:4);
            
                % Reflective Boundary Conditions
                if dest(1)> VoxLength
                    dest(1)=1;
                elseif dest(1)==0
                    dest(1)=VoxLength;
                end

                if dest(2)> VoxWidth
                    dest(2)=1;
                elseif dest(2)==0
                    dest(2)=VoxWidth;
                end
                
                if dest(3)> VoxHeight
                    dest(3)=1;
                elseif dest(3)==0
                    dest(3)=VoxHeight;
                end
                
                if CrowdSpace(dest(1),dest(2),dest(3))==0
    
                    RibosomeSpace(xVox,yVox,zVox) = RibosomeSpace(xVox,yVox,zVox) - 1;  %decrement conc of Rib in source voxel
                    RibosomeSpace(dest(1),dest(2),dest(3)) = RibosomeSpace(dest(1),dest(2),dest(3)) + 1; %increment conc of Rib in destination voxel

                    % Take the first ribosome found at source voxel and change address to
                    % destination, and update RibosomeCurrent
                    chooseR = find(RibosomeCurrent(:,1)==xVox & RibosomeCurrent(:,2)==yVox & RibosomeCurrent(:,3)==zVox, 1,'first'); 
                    RibosomeCurrent(chooseR,:)=dest;
                else
                    %reject! no movement
                    rReject=rReject+1;
                end
        else
            % something has gone wrong with the number of reactions
            disp('something odd happened with q, your code sucks')
        end
        end
        
        while (T >= RecordCount*(dt))   %sampled time
            mRNATrack2(RecordCount) = sum(mRNASpace(:));
            mRNATrack(:,:,RecordCount)= mRNACurrent;
            mRNASpaceTrack(:,:,:,RecordCount)=mRNASpace;
            RibosomeTrack(:,:,RecordCount)= RibosomeCurrent;
            RiboSpaceTrack(:,:,:,RecordCount)= RibosomeSpace;
            RecordCount = RecordCount + 1;
        end
    RxnCount=RxnCount+1;
    end

    
save(Name);
%%
%%            
%% MODEL ENDS HERE            
        end
    end
end
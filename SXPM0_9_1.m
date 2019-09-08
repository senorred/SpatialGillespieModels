%% Spatial Gillespie Model for XMAS Crowding, 0.9.1
% Purpose: The purpose of this model is to demonstrate the spatial
% accumulation of ribosomes around early mRNAs, and the resultant impact on
% protein production in these regions. The model will later be manipulated
% by introducing crowders to the space. 
%
% Description: A Spatial 2D Gillespie model. Based on Chuck's
% SpatialResourceCluster.m file.
%
% Version Description: second try!

%% Generate Seed
clear all
close all

rand('state',sum(100*clock));

jobNum=1

%% Model Description

%Turn on commenting display? Mainly for debugging
Comments = 0;
%Turn on gene circuit?
Expression = 1;
%reflective[0] or periodic[1] boundaries
Boundary = 0;
%Sticky effects on diffusion?
StickyEffects = 1;
%Steady State[0] or Transient[1]
Experiment = 1;

%% Parameters

%Running variables
tmax = 100;
dt = 1; %step to record state
tspan = 0:dt:tmax;
Runs = 1;

%Rate Constants for expression
kONStart = 0.1;
kOFFStart = 10;

%changing kON
kOFFArray = [.1,10];%ones(1,2)*kOFFStart;
kONArray = [.1,10];

StdArray = [0.5,1,2,5];
%changing kOFF
%kOFFArray = [10,.1];
%kONArray = ones(1,2)*kONStart;

kM = 10; %k*V
kP = 0;
gammam = .1;
gammap = 0;
lambda = 0;

varRun = ceil(jobNum/20);
kOFF = kOFFArray(varRun);
kON = kONArray(varRun);

ONfrac = kON / (kON + kOFF);

%Species
TXResources = 1;
Ribosomes = 0; 
MaxGenes = 2;
InitialGenes = 0;
VaryMean = 0;
InitialmRNA = 0;
Protein = 0;


%Space
VoxLength = 15; %number of voxels along x direction
VoxWidth = 15; %number of voxles along y direction
VoxHeight = 1;
h = 1; %Side length in um

%Diffusion Constants
DRes = 60; %Resource Bulk Diffusion Constant
DRib = 60; %Ribosome Bulk Diffusion Constant
DGene = 0; %Gene Bulk Diffusion Constant
DRibA = .60; %Ribosome Active Diffusion Constant
DResA = .60; %Resource Active Diffusion Constant

PullFactor = .5;

%Crowders
Crowding = 0; %0 no crowding, 1 crowding
DynamicCrowding = 0; %0 static crowders, 1 moving crowders
CrowdVol = 0; % percent volume crowded
CrowdD = .0006; %Crowder Diffusion Constant
CrowdWallD = .0006; %Crowder Diffusion Constant on Wall

Name = sprintf('%gx%gx%gkON%gkOFF%gS%gSSResults%g.mat',VoxLength,VoxWidth,VoxHeight,...
    kON,kOFF,StickyEffects,jobNum);
%% Initial Conditions
tic
%Batch running%


%Inital calculations
ActiveGenes = zeros(1,MaxGenes);
ActiveGenes(1:InitialGenes) = 1;
Volume = (VoxLength*VoxWidth*VoxHeight)*h^3;
%kM = kM*(Volume/h^3);
if CrowdVol == 0
    NumCrowders = 0;
else
    NumCrowders = floor(VoxLength*VoxWidth*VoxHeight*(CrowdVol/100));
end
MaxPassage = 5000; %expected maximum number of passages per particle

%kM = kM/h^3; %effective kM second order reaction

% %Initialize seperate species tracking variables
%disp('initializing')
MaxGeneBurst = 100000;
AllGene1Track = zeros(MaxGeneBurst,2,Runs);
Gene1ONTimesTrack = zeros(MaxGeneBurst,1,Runs);
AllGene2Track = zeros(MaxGeneBurst,2,Runs);
Gene2ONTimesTrack = zeros(MaxGeneBurst,1,Runs);
Gene1RhitTrack = zeros(MaxGeneBurst,1,Runs);
Gene2RhitTrack = zeros(MaxGeneBurst,1,Runs);
Gene1TTRTrack =  zeros(MaxGeneBurst,1,Runs);
Gene2TTRTrack =  zeros(MaxGeneBurst,1,Runs);

%% Initialize Model

for ii = 1:Runs
    disp(ii)
    
    %Initialize Gillespie%


    %Tracking Arrays for storing data
    %Protein is not tracked spatially, can just be placed in single array
    ProteinTrack = zeros(1,length(tspan));
    GeneONTrack = ProteinTrack;
    Gene1Track = zeros(MaxGeneBurst,2);

    Gene2Track = Gene1Track;
    Gene1Track(1,:) = [0,0];
    Gene2Track(1,:) = [0,0];
    Gene1TrackCount = 1;
    Gene2TrackCount = 1;
    Gene1RhitCount = 1;
    Gene2RhitCount = 1;
    Gene1TTRCount = 1;
    Gene2TTRCount = 1;
    Gene1ONFlag = 0;
    Gene2ONFlag = 0;
    Gene1CurrentFlag = 0;
    Gene2CurrentFlag = 0;
    
    ProteinTrack(1) = Protein;%round(Protein + randn*sqrt(Protein)); %First point is initial protein count


    mRNATrack = zeros(1,length(tspan));
    mRNATrack(1) = InitialmRNA;


    ActiveRibosomeTrack = zeros(1,length(tspan));
    OldEfficiencyTrack = zeros(1,length(tspan));

    %SpatialRxnMatrix
    SpatialRxnMatrix = [ 1  0  0; %move up 1 in x
                        -1  0  0; %move down 1 in x
                         0  1  0; %move up 1 in y
                         0 -1  0; %move down 1 in y
                         0  0  1; %move up 1 in z
                         0  0 -1]; %move down 1 in z
    %NonSpatial RxnMatrix
    %mRNA, Protein, Gene ON, Gene OFF
    %RxnMatrix = 1;%transcription
     RxnMatrix = [ 1  0  0  0; %transcription
                   0  1  0  0; %translation
                  -1  0  0  0; %mRNA decay
                   0 -1  0  0; %Protein decay
                   0  0 -1  1; %burst OFF
                   0  0  1 -1];%burst ON

    %track spatial species on an individual basis
    %rows,for each species, xyz location, for all time.
    RibosomeTrack = zeros(Ribosomes,3,length(tspan));
    ResourceTrack = zeros(TXResources,3,length(tspan));
    GeneTrack = zeros(MaxGenes,3,length(tspan));
    ActiveGeneTrack = zeros(length(ActiveGenes),length(tspan));
    CrowdTrack = zeros(NumCrowders,3,length(tspan));

    PassageDurationTrack = zeros(Ribosomes,MaxGenes,MaxPassage);
    PassageDurationLocTrack = zeros(Ribosomes,MaxGenes,MaxPassage);
    PassageTimesTrack = zeros(Ribosomes,MaxGenes,MaxPassage);
    PassageCounterMatrix = ones(Ribosomes,MaxGenes);
    PassageLocationsMatrix = zeros(VoxLength, VoxWidth, VoxHeight);
    FirstPassageMatrix = NaN(Ribosomes,MaxGenes);
    MeetMatrix = zeros(Ribosomes,MaxGenes);
    RibosomeWallhit = zeros(Ribosomes,1);
    RibosomeVoxDuration = zeros(VoxLength,VoxWidth,VoxHeight);
    RibosomeVoxCurrentTime = zeros(Ribosomes,1);
    TranscriptWallhit = zeros(MaxGenes,1);
    RibosomeCrowderhit = zeros(Ribosomes,1);
    TranscriptCrowderhit = zeros(MaxGenes,1);
    %Species tracked on a compartment basis. 2D 1to1 map xy, z species, k over
    %time
    SpeciesSpace = zeros(VoxLength,VoxWidth,VoxHeight,3);
    ResourceSpace = zeros(VoxLength,VoxWidth,VoxHeight);
    RibosomeSpace = zeros(VoxLength,VoxWidth,VoxHeight);
    GeneActiveSpace = zeros(VoxLength,VoxWidth,VoxHeight);
    GeneONSpace = zeros(VoxLength,VoxWidth,VoxHeight);
    GeneOFFSpace = zeros(VoxLength,VoxWidth,VoxHeight);
    mRNASpace = zeros(VoxLength,VoxWidth,VoxHeight);
    ProteinSpace = zeros(VoxLength,VoxWidth,VoxHeight);

    %Track propensities across each compartment space
    aSpace = zeros(4,1);

    %Distribute crowders in the system
    CrowdSpace = zeros(VoxLength,VoxWidth,VoxHeight);
    MaxTrials = 500;

    %Randomly distribute Species into compartments
    for i = 1:TXResources
        comp = ceil([rand*VoxLength, rand*VoxWidth, rand*VoxHeight]);
        ResourceSpace(comp(1),comp(2),comp(3)) = ResourceSpace(comp(1),comp(2),comp(3)) + 1;
        ResourceTrack(i,:,1) = comp;
    end

    for i = 1:Ribosomes
        comp = ceil([rand*VoxLength, rand*VoxWidth, rand*VoxHeight]);
        RibosomeSpace(comp(1),comp(2),comp(3)) = RibosomeSpace(comp(1),comp(2),comp(3)) + 1;
        RibosomeTrack(i,:,1) = comp;
    end

    for i = 1:MaxGenes
        comp = ceil([rand*VoxLength, rand*VoxWidth, rand*VoxHeight]);
        GeneOFFSpace(comp(1),comp(2),comp(3)) = GeneOFFSpace(comp(1),comp(2),comp(3)) + 1;
        GeneTrack(i,:,1) = comp;
        GeneIdx(i,1) = comp(1);
        GeneIdx(i,2) = comp(2);
    end
    if MaxGenes > 1
        GeneDist = sqrt((GeneIdx(2,2) - GeneIdx(1,2))^2 + (GeneIdx(2,1) - GeneIdx(1,1))^2);
    end

    SpaceIdx = RibosomeSpace + GeneOFFSpace;
    SpaceIdx(SpaceIdx > 0) = 1;

    mRNACurrent = mRNATrack(1);
    ProteinCurrent = ProteinTrack(1);
    ResourceCurrent = ResourceTrack(:,:,1);
    RibosomeCurrent = RibosomeTrack(:,:,1);
    GeneCurrent = GeneTrack(:,:,1);
    CrowdCurrent = CrowdTrack(:,:,1);

    x = [InitialmRNA, Protein];
    TTrack = zeros(length(tspan), 1); %Time tracking
    T = 0; %Start Time
    RxnCount = 1;  %Reaction counter
    tempBurstCount = 0;
    RecordTime = dt; %Recording time
    RecordCount = 1;

    a = zeros(VoxLength,VoxWidth,VoxHeight,6);
    amn = zeros(VoxLength,VoxWidth,VoxHeight);
    kPo = kP;

    %% Spatial Model Run
    
    Timer = 0;
    disp('Gillespie Start')
    while T < tmax %Run gillespie until time is up

        %Timer for console display
        if floor(T) > Timer
            Timer = Timer + 1;
            disp(Timer)
        end

        %calculate Propensities
        kP = kPo*exp(-T*lambda);
        for x = 1:VoxLength
            for y = 1:VoxWidth
                for z = 1:VoxHeight
                    if ResourceSpace(x,y,z) >= 1
                        tempIdx = 1;
                    else
                        tempIdx = 0;
                    end
                    a(x,y,z,1) = kM*GeneONSpace(x,y,z)*tempIdx; %compute propensities in voxel (x,y,z)

                    if mRNASpace(x,y,z) >= 1 %if mRNA present, produce protein at rate kP
                        a(x,y,z,2) = kP*RibosomeSpace(x,y,z);
                    else 
                        a(x,y,z,2) = 0;
                    end
                    a(x,y,z,3) = gammam*mRNASpace(x,y,z); %compute propensities in voxel (x,y,z)
                    a(x,y,z,4) = gammap*ProteinSpace(x,y,z);
                    a(x,y,z,5) = kOFF*GeneONSpace(x,y,z); %burst OFF
clear                    a(x,y,z,6) = kON*GeneOFFSpace(x,y,z); %Burst ON
                    if (mRNASpace(x,y,z) >= 1) %if mRNA present, diffusion much harder
                        a(x,y,z,7) = DRibA*RibosomeSpace(x,y,z);
                    else
                        a(x,y,z,7) = DRib*RibosomeSpace(x,y,z);
                    end
                    if (GeneONSpace(x,y,z) >= 1) %if gene on present, diffusion much harder                       
                        a(x,y,z,8) = DResA*ResourceSpace(x,y,z);
                    else
                        a(x,y,z,8) = DRes*ResourceSpace(x,y,z);
                    end

                    amn(x,y,z) = sum(a(x,y,z,:));

                end
            end
        end
        a0 = sum(amn(:));

        rtime = rand();
        rvox = rand();
        ra = rand();
        rmove = rand();

        tau = (1/(a0))*log(1/ rtime);  %time-step increment
        T = T + tau;   %increment current time by tau
        ajk = 0;
        flag = 0;
        for j = 1:VoxLength   %select voxel (j,k,l)
            for k = 1:VoxWidth
                for l = 1:VoxHeight
                    ajk = ajk + amn(j,k,l);
                    if ajk > rvox*(a0)
                        flag = 1;
                        break;
                    end
                end
                if (flag == 1) 
                    break;
                end
            end
            if (flag == 1) 
                break;
            end
        end
        q = 1;
        ajkl = a(j,k,l,q);
        while ajkl <= ra*(amn(j,k,l))
            q = q + 1;
            ajkl = ajkl + a(j,k,l,q);
        end
        if (q <= 6)   %reaction chosen
            mRNASpace(j,k,l) = mRNASpace(j,k,l) + RxnMatrix(q,1);
            ProteinSpace(j,k,l) = ProteinSpace(j,k,l) + RxnMatrix(q,2);
            GeneONSpace(j,k,l) = GeneONSpace(j,k,l) + RxnMatrix(q,3);
            GeneOFFSpace(j,k,l) = GeneOFFSpace(j,k,l) + RxnMatrix(q,4);

            %Track bursting
            if q == 1
                currentSpace = [j,k,l];
                if GeneIdx(i,1) == j && GeneIdx(i,2) == k
                    Gene1Track(Gene1TrackCount,1) = Gene1Track(Gene1TrackCount,1)+1;
                else
                    Gene2Track(Gene2TrackCount,1) = Gene2Track(Gene2TrackCount,1)+1;
                end
            end
            if q == 5
                currentSpace = [j,k,l];
                if GeneIdx(i,1) == j && GeneIdx(i,2) == k
                    Gene1ONFlag = 0;
                    Gene1CurrentFlag = 0;
                else
                    Gene2ONFlag = 0;
                    Gene2CurrentFlag = 0;
                end
            end
            if q == 6
                currentSpace = [j,k,l];
                if GeneIdx(i,1) == j && GeneIdx(i,2) == k
                    Gene1Track(Gene1TrackCount,2) = T;
                    Gene1ONFlag = 1;
                    Gene1CurrentFlag = 1;
                    Gene1TrackCount = Gene1TrackCount + 1;
                else
                    Gene2Track(Gene2TrackCount,2) = T;
                    Gene2ONFlag = 1;
                    Gene2CurrentFlag = 1;
                    Gene2TrackCount = Gene2TrackCount + 1;
                end
            end

        elseif (q == 7)  %particle move chosen
            reject = 0;
            dir = ceil(6 * rmove);
            m = j + SpatialRxnMatrix(dir,1);  %x-index of target voxel
            n = k + SpatialRxnMatrix(dir,2);  %y-index of target voxel
            p = l + SpatialRxnMatrix(dir,3);  %z-index of target voxel


            if (m > VoxLength)  
                reject = 1;
            end %hard wall boundary conditions
            if (m < 1) 
                reject = 1;
            end
            if (n > VoxWidth) 
                reject = 1;
            end
            if (n < 1) 
                reject = 1;
            end
            if (p > VoxHeight) 
                reject = 1;
            end
            if (p < 1) 
                reject = 1;
            end

            if (reject == 0)

                if (CrowdSpace(m,n,p) == 1)  % if crowder at target voxel
                else

                RibosomeSpace(j,k,l) = RibosomeSpace(j,k,l) - 1;  %decrement conc of Rib in vox (j,k)
                RibosomeSpace(m,n,p) = RibosomeSpace(m,n,p) + 1; %increment conc of Rib in vox (m,n)
                end
            end
        else
            reject = 0;
            dir = ceil(6 * rmove);
            m = j + SpatialRxnMatrix(dir,1);  %x-index of target voxel
            n = k + SpatialRxnMatrix(dir,2);  %y-index of target voxel
            p = l + SpatialRxnMatrix(dir,3);  %z-index of target voxel


            if (m > VoxLength)  
                reject = 1;
            end %hard wall boundary conditions
            if (m < 1) 
                reject = 1;
            end
            if (n > VoxWidth) 
                reject = 1;
            end
            if (n < 1) 
                reject = 1;
            end
            if (p > VoxHeight) 
                reject = 1;
            end
            if (p < 1) 
                reject = 1;
            end

            if (reject == 0)

                if (CrowdSpace(m,n,p) == 1)  % if crowder at target voxel
                else

                ResourceSpace(j,k,l) = ResourceSpace(j,k,l) - 1;  %decrement conc of Rib in vox (j,k)
                ResourceSpace(m,n,p) = ResourceSpace(m,n,p) + 1; %increment conc of Rib in vox (m,n)
                ResourceCurrent(1,:) = [m;n;p];
                currentSpace = [m,n,p];
                if GeneIdx(1,1) == m && GeneIdx(1,2) == n
                    Gene1RhitTrack(Gene1RhitCount,1) = T;
                    Gene1RhitCount = Gene1RhitCount + 1;
                    if Gene1ONFlag == 1 && Gene1CurrentFlag == 1
                        Gene1TTRTrack(Gene1TTRCount) = T - Gene1Track(Gene1TrackCount-1,2);
                        Gene1CurrentFlag = 0;
                        Gene1TTRCount = Gene1TTRCount + 1;
                    end
                elseif GeneIdx(2,1) == m && GeneIdx(2,2) == n
                    Gene2RhitTrack(Gene2RhitCount,1) = T;
                    Gene2RhitCount = Gene2RhitCount + 1;
                    if Gene2ONFlag == 1 && Gene2CurrentFlag == 1
                        Gene2TTRTrack(Gene2TTRCount) = T - Gene2Track(Gene2TrackCount-1,2);
                        Gene2CurrentFlag = 0;
                        Gene2TTRCount = Gene2TTRCount + 1;
                    end
                end
                end
            end
        end



        while (T >= RecordCount*(dt))   %sampled time
            mRNATrack(RecordCount) = sum(mRNASpace(:));
            ProteinTrack(RecordCount) = sum(ProteinSpace(:));
            GeneONTrack(RecordCount) = sum(GeneONSpace(:));
            Gene1OFFTrack(RecordCount) = GeneOFFSpace(GeneIdx(1,1),GeneIdx(1,2),1);
            Gene2OFFTrack(RecordCount) = GeneOFFSpace(GeneIdx(2,1),GeneIdx(2,2),1);
            Gene1ONTrack(RecordCount) = GeneONSpace(GeneIdx(1,1),GeneIdx(1,2),1);
            Gene2ONTrack(RecordCount) = GeneONSpace(GeneIdx(2,1),GeneIdx(2,2),1);
            ResourceTrack(:,:,RecordCount) = ResourceCurrent;
            RecordCount = RecordCount + 1;
        end

    end %end gillespie
    
    if Gene1TrackCount-1 > 0
        AllGene1Track(:,:,ii) = Gene1Track;
    end
    if Gene2TrackCount-1 > 0
        AllGene2Track(:,:,ii) = Gene2Track;
    end
    AllGeneDist(ii) = GeneDist;
    AllmRNATrack(:,ii) = mRNATrack;

end %end all runs
% mRNATrack(end) = [];
% ProteinTrack(end) = [];
% GeneONTrack(end) = [];
ElapsedTime = toc;
%save(Name);

save(Name,'mRNATrack','Gene1Track','Gene2Track','ONfrac','GeneIdx','ResourceTrack',...
    'Gene1RhitTrack','Gene2RhitTrack','Gene1TTRTrack','Gene2TTRTrack','GeneDist');
quit();


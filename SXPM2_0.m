%% Spatial XMAS Pocket Model 2.0
% Purpose: To describe the burst size and burst frequency over time for
% mRNAs acquiring ribosomes
%
% Description: Spatial gillespie model that tracks mRNA birth times and
% ribosome positioning over time. Can do either uniform distribution of
% ribosomes or randomly seed them.
%
% Version Description: "More accurate" version which reads through all
% possible reactions to find the reaction voxel, instead of reading the
% most likely voxel and running through the possible reactions

%% Provide seed
% seed
clear all
close all
rand('state',sum(100*clock)); %#ok<RAND>

%% Specify Parameters
% runtime
tmax= 1000;  %100
dt = 1; %step to record state
tspan = 0:dt:tmax;

% grid size
VoxLength=20;
VoxWidth=20;

% Ribo/Polysome diffusion
dRibo = 1.0;    % 1.0
dPoly = 0.0001; % 0.0001
% Translation rate
kP=1.0;

% Species
Ribosomes= 400;
InitialmRNA=0;
mRNACount= 0;
mRNATrack2=[];

% Crowders
CrowdVol= 50;
NumCrowders = floor(VoxLength*VoxWidth*(CrowdVol/100));
rReject=0;

% mRNA birth rate
kCM= 1/((100-CrowdVol)/100);
alpha= kCM*0.001;   %0.001

%% Arrays for storage
% population of ribos per voxel
% coords of mRNAs, ribos
    mRNASpace = zeros(VoxLength,VoxWidth);
    RibosomeSpace = zeros(VoxLength,VoxWidth);
    CrowdSpace= zeros(VoxLength,VoxWidth);

%% Initialize Model

    %track spatial species on an individual basis
    %rows,for each species, xyz location, for all time.
    RibosomeTrack = zeros(Ribosomes,2,length(tspan));
    mRNATrack = zeros(1,2,length(tspan));
    CrowdTrack= zeros(NumCrowders, 2,length(tspan));
    
    mRNASpaceTrack=zeros(VoxLength, VoxWidth, length(tspan));
    RiboSpaceTrack=zeros(VoxLength, VoxWidth, length(tspan));
    
            %Randomly distribute Crowders into compartments, only one
            %crowder per compartment
    i=1;        
    while i <= NumCrowders
        comp = ceil([rand*VoxLength, rand*VoxWidth]);
        if CrowdSpace(comp(1),comp(2))== 0
            CrowdSpace(comp(1),comp(2)) = CrowdSpace(comp(1),comp(2)) + 1;
            CrowdTrack(i,:,1) = comp;
            i=i+1;
        else
        end
    end
    
        %Randomly distribute Ribos into compartments NOT OCCUPIED BY
        %CROWDERS
    i=1;    
    while i <= Ribosomes
        comp = ceil([rand*VoxLength, rand*VoxWidth]);
        if CrowdSpace(comp(1),comp(2))==0
            RibosomeSpace(comp(1),comp(2)) = RibosomeSpace(comp(1),comp(2)) + 1;
            RibosomeTrack(i,:,1) = comp;
            i=i+1;
        else 
        end
    end

% % uniform dist of ribos (MUST ADD CROWDER EXCEPTION)
%     RibosomeSpace=ones(20,20);
%     A=1:400;
%     [I,J]=ind2sub([20,20],A');
%     RibosomeTrack(:,1,1)=I;
%     RibosomeTrack(:,2,1)=J;
%     RiboInit=RibosomeSpace;

RiboTrackInit=RibosomeTrack(:,:,1);
CrowdTrackInit=CrowdTrack(:,:,1);

%SpatialRxnMatrix [ mRNA X Y]
SpatialRxnMatrix = [ 1  0  0; %mRNA born
                     0  1  0; %move ribo up 1 in x
                     0 -1  0; %move ribo down 1 in x
                     0  0  1; %move ribo up 1 in y
                     0  0 -1]; %move ribo down 1 in y

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

a = zeros(VoxLength,VoxWidth,5);
amn = zeros(VoxLength,VoxWidth);

%% Model
    
    Timer = 0;
    %disp('Gillespie Start')
    while T < tmax %Run gillespie until time is up

        %Timer for console display
        if floor(T) > Timer
            Timer = Timer + 1;
            disp(Timer)
        end

        %calculate Propensities
        for x = 1:VoxLength
            for y = 1:VoxWidth
                
                    a(x,y,1)= alpha;
                    
                    if mRNASpace(x,y) >= 1
                        a(x,y,2:5) = dPoly*RibosomeSpace(x,y);
                    else
                        a(x,y,2:5) = dRibo*RibosomeSpace(x,y);
                    end
                    
                    % Since crowders don't move, I assume this is okay...
                    if CrowdSpace(x,y)==1
                        a(x,y,:)=0;
                    end
            end
        end
        a0 = sum(a(:));

        rtime = rand();
        ra = rand();

        tau = (1/(a0))*log(1/ rtime);  %time-step increment
        T = T + tau;   %increment current time by tau
        flag = 0;
        
        % Choose a reaction, get the X and Y coordinates and the reaction
        % type, 'q'
        q=[];
        avec=reshape(a,[1,(VoxLength*VoxWidth*5)]);
        ia = find((cumsum(avec) >= ra*a0),1,'first');
        [xVox, yVox, q]=ind2sub(size(a),ia);
        
        % "Does a crowder exist here" gut check-- shouldn't 
        if CrowdSpace(xVox,yVox)>0
            disp('whoa partner, your code sucks!')
        else
        % According to which reaction type is chosen, execute reaction
        if q == 1   %makes mRNA in chosen voxel
            mRNASpace(xVox,yVox) = mRNASpace(xVox,yVox) + 1;
            mRNACount= mRNACount+1;
            disp('Made mRNA!')
            mRNACurrent=vertcat(mRNACurrent, [xVox, yVox]);
            mRNATrack= vertcat(mRNATrack, zeros(1,2,length(tspan)));

        elseif q>=2  %moves Ribosome from chosen voxel
            % calculate destination coordinates
            dest= [xVox,yVox]+ SpatialRxnMatrix(q,2:3);
            
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
                
                if CrowdSpace(dest(1),dest(2))==0
    
                    RibosomeSpace(xVox,yVox) = RibosomeSpace(xVox,yVox) - 1;  %decrement conc of Rib in source voxel
                    RibosomeSpace(dest(1),dest(2)) = RibosomeSpace(dest(1),dest(2)) + 1; %increment conc of Rib in destination voxel

                    % Take the first ribosome found at source voxel and change address to
                    % destination, and update RibosomeCurrent
                    chooseR = find(RibosomeCurrent(:,1)==xVox & RibosomeCurrent(:,2)==yVox, 1,'first'); 
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
            mRNASpaceTrack(:,:,RecordCount)=mRNASpace;
            RibosomeTrack(:,:,RecordCount)= RibosomeCurrent;
            RiboSpaceTrack(:,:,RecordCount)= RibosomeSpace;
            RecordCount = RecordCount + 1;
        end
    RxnCount=RxnCount+1;
    end

function [gainOverNoisedB,R,pilotIndex,APpositions,UEpositions,distance_bw_UEs_APs] = generateSimulationSetup3(L,K,N,paramValues,chDist,shadowFlag)
%Generate realizations of the simulation setup described in:
%
%Zakir Hussain Shaik, Emil Bjornson, and Erik G. Larsson,
%"MMSE-Optimal Sequential Processing for Cell-Free Massive MIMO With Radio
%Stripes," IEEE Transactions on Communications, To appear.
%
%Download article: https://arxiv.org/pdf/2012.13928.pdf
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.
%
%INPUT:
%L                  = Number of APs for the Cell-free system
%K                  = Number of UEs in the network
%N                  = Number of antennas per AP
%
%OUTPUT:
%R                  = Matrix with dimension N x N x K x L
%                     where (:,:,k,l) is the spatial correlation matrix
%                     between AP l and UE k, normalized by noise
%                     power
%pilotIndex         = Matrix with dimension K x 1 containing the
%                     pilot assigned to the UEs for a given setup
%gainOverNoisedB    = Matrix with dimension K x L where
%                     element (k,l) is the channel gain (normalized by
%                     the noise power) between AP l and UE k
%UEpositions        = Vector of length K with UE positions, where the real
%                     part is the horizontal position and the imaginary
%                     part is the vertical position
%APpositions        = Vector of length L with the AP locations, measured in
%                     the same way as UEpositions


%% Define simulation setup
% For reference, Parameters indexing is as below
%paramValues(1)  : Bandwidth
%paramValues(2)  : Noise Figure
%paramValues(3)  : Coherence block length/Time
%paramValues(4)  : Pilot Time or number of pilots per coherence block
%paramValues(5)  : Power allocated to UEs
%paramValues(6)  : Alpha - path loss constant
%paramValues(7)  : PL_constantTerm
%paramValues(8)  : sigma_sf
%paramValues(9)  : decorr
%paramValues(10) : distanceVertical
%paramValues(11) : antennaSpacing
%paramValues(12) : ASDdeg
%paramValues(13) : squareLength

%{B,noiseFigure,tau_c,tau_p,allocatedPowUEs,alpha,PL_constantTerm,sigma_sf,decorr,distanceVertical,antennaSpacing,ASDdeg,radioStripeLength};
%Communication bandwidth
B = paramValues{1};

%Noise figure (in dB)
noiseFigure = paramValues{2};

%Compute noise power:
%noiseVariancedB = -204 + 10*log10(B) + noiseFigure
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Length of the coherence block
%tau_c = paramValues{3};

%Number of pilots per coherence block
tau_p = paramValues{4};

%Uplink transmit power per UE (in mW)
%pow_UE = paramValues{5};

%Pathloss parameters
alpha = paramValues{6}; % Pathloss exponent
constantTerm = paramValues{7}; % Median channel gain at a reference distance of 1Km in (dB)

%Standard deviation of the shadow fading
sigma_sf = paramValues{8};

%Decorrelation distance of the shadow fading
decorr = paramValues{9}; % Distance beyond which UEs are decorrelated

%Height difference between an AP and a UE
distanceVertical =  paramValues{10};

%Define the antenna spacing (in number of wavelengths)
antennaSpacing =  paramValues{11};%Half wavelength distance

%Angular standard deviation around the nominal angle (measured in degrees)
ASDdeg =  paramValues{12};

%Size of the coverage area (as a square with wrap-around)
radioStripeLength = paramValues{13}; %meter

%Per side length of the setup, it can be greater than radio stripe length/4
perSideLength = radioStripeLength/4;

pilotIndex = zeros(K,1);
% if sum(abs(diff(APpositionss))<antennaSpacing)>0
%
%     error('Radio stripe cannot accommodate APs with half wavelength spacing for given L and radio stripe length');
%
% end

%Number of APs per dimension on the grid
%nbrAPsPerDim = L/4;

%Distance between BSs in vertical/horizontal direction
%interAPDistance = radioStripeLength/nbrAPsPerDim;

%% Deploy APs on the grid
% singleAxisAPpositions = linspace(0,perSideLength,L/4+1);
% APpositions = zeros(L,1);
% APpositions(1:L/4+1,1) = singleAxisAPpositions.';
% APpositions(1+L/4:2*L/4+1,1) = perSideLength + 1i*singleAxisAPpositions.';
% APpositions(1+2*L/4:3*L/4+1,1) = flipud(singleAxisAPpositions.') + 1i*perSideLength;
% APpositions(1+3*L/4:L,1) = 1i*flipud(singleAxisAPpositions(2:end).');

APLocations = linspace(0,radioStripeLength-1,L+1); APpositions=APLocations;

APpositions(APLocations<=perSideLength) = 1i*APpositions(APLocations<=perSideLength);
APpositions(perSideLength<APLocations & APLocations<=2*perSideLength) = APpositions(perSideLength<APLocations & APLocations<=2*perSideLength)-(perSideLength) + 1i*(perSideLength);
APpositions(2*perSideLength<APLocations & APLocations<=3*perSideLength) = (perSideLength) + 1i*((3*perSideLength) - APpositions(2*perSideLength<APLocations & APLocations<=3*perSideLength));
APpositions(3*perSideLength<APLocations) = ((4*perSideLength) - APpositions(3*perSideLength<APLocations));
APpositions(end) = [];

APpositions = APpositions.';


%% Deploy UEs on the grid randomly
%Random UE locations with uniform distribution
UEpositions = (0.1 + (0.9-0.1).*rand(K,1))*perSideLength + 1i*((0.1 + (0.9-0.1).*rand(K,1))* perSideLength) ;

% Uncomment below code if one wants display generated simulation setup
% tmp = 0:perSideLength;
% APtmp = [1i*tmp,tmp + 1i*perSideLength, perSideLength + 1i*flip(tmp),tmp];
% plot(APtmp,'b--');
% hold on;
% plot(APpositions,'ro');
% hold on;
% plot(UEpositions,'g*');
% axis([-2 perSideLength+1 -2 perSideLength+1]);
% grid on;
% hold off;

%Prepare to save results
R = zeros(N,N,K,L);

%Compute the correlation matrices for the shadow fading

%Distance between UEs
distance_bw_UEs = abs(repmat(UEpositions,[1,K]) - repmat(transpose(UEpositions),[K,1])); %No Wrap around
shadowCorrMatrix_UEs = 2.^(-distance_bw_UEs/decorr);

%Generate shadow fading realizations
f = sigma_sf*sqrtm(shadowCorrMatrix_UEs)*randn(K,L);

%Simulation Models
% No wrap around
distance_bw_UEs_APs_inPlane0 = (repmat(UEpositions,[1,L]) - repmat(transpose(APpositions),[K,1])); % Vector position difference useful to compute angle in line 175 of the code
distance_bw_UEs_APs_inPlane = abs(repmat(UEpositions,[1,L]) - repmat(transpose(APpositions),[K,1])); % Distance
distance_bw_UEs_APs = sqrt(distanceVertical^2+distance_bw_UEs_APs_inPlane.^2);

%Terminal powers is in dBm
gainOverNoisedB = constantTerm - alpha*10*log10(distance_bw_UEs_APs) - noiseVariancedBm + 30;

if shadowFlag==1
    % Add shadow fading to all channels from APs to UE k that have a distance larger than 50 meters
    gainOverNoisedB(distance_bw_UEs_APs>50) = gainOverNoisedB(distance_bw_UEs_APs>50) + f(distance_bw_UEs_APs>50);
    
end


% Nominal Angle between UEs to APs
nAngle_bw_UEs_APs = angle(distance_bw_UEs_APs_inPlane0);

for k = 1:K
    for l = 1:L
        
        if strcmp(chDist,'uncorrelated')
            
            correlationMatrix = eye(N);
            
        else
            correlationMatrix = functionRlocalscattering(N,nAngle_bw_UEs_APs(k,l),ASDdeg,antennaSpacing,chDist);
        end
        
        R(:,:,k,l) = db2pow(gainOverNoisedB(k,l))*correlationMatrix;
        
    end
    
    %Generate pilot indexing
    %Determine the best AP for UE k by looking for AP with best
    %channel condition
    [~,bestAP] = max(gainOverNoisedB(k,:));
    
    if k <= tau_p
        
        pilotIndex(k) = k;
        
    else %Assign pilot for remaining UEs
        
        %Compute received power from to the master AP from each pilot
        pilotinterference = zeros(tau_p,1);
        
        for t = 1:tau_p
            
            pilotinterference(t) = sum(db2pow(gainOverNoisedB(pilotIndex(1:k-1,1)==t,bestAP)));
            
        end
        
        %Find the pilot with the least receiver power
        [~,bestpilot] = min(pilotinterference);
        pilotIndex(k,1) = bestpilot;
        
    end
    
    
end



end

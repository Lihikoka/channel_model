% function [Channel, Fading] = MIMO_channel(FadingMatrixType, ... 
% 			     FadingOversamplingFactor, Speed_ms, ... 
% 			     CarrierFrequency_Hz, ... 
% 			     RunningTimeInstant, ChipOversamplingFactor, ... 
% 			     ChipRate_Hz, NumberOfChipsPerIteration, PDP_linear,... 
%                  NumberOfTxAntennas, NumberOfTransmittedSamples, ... 
%                  NumberOfRxAntennas, NumberOfReceivedSamples) 
              
function [Channel, Fading] = MIMO_channel(FadingMatrixType, ... 
			      FadingOversamplingFactor, Speed_ms, ... 
		 	      CarrierFrequency_Hz, ... 
	 		      RunningTimeInstant, ChipOversamplingFactor, ... 
			      ChipRate_Hz, NumberOfChipsPerIteration, PDP_linear,... 
                NumberOfTxAntennas, NumberOfTransmittedSamples, ... 
                NumberOfRxAntennas, NumberOfReceivedSamples)
% 
% Generates the correlated tap coefficients of the MIMO tapped delay line 
% model to be used during one iteration of the main loop. The function 
% performs a double interpolation, first in the fading vector domain, 
% to collect the fading samples corresponding to the (sub)chip-spaced time 
% samples, then in the tap domain, to match the delays of the PDP to 
% the sample instants of the simulation. This second interpolation complies 
% with the description in TR 25.869 v 0.1.1, p.23. 
% 
% The correlated tap coefficients are delivered through two different 
% matrices Channel and Fading which embed the same information but 
% exhibit two different structures. 
% 
% Inputs 
% 
% * 2-D matrix FadingMatrixType of size (NumberOfPaths * NumberOfTxAntennas 
%   * NumberOfRxAntennas) x FadingNumberOfIterations containing the spatially 
%   correlated fading process of the NumberOfPaths * NumberOfTxAntennas 
%   * NumberOfRxAntennas tap coefficients of the MIMO model, including 
%   an optional Rice component. 
% * Variable FadingOversamplingFactor, oversampling factor of the fading 
%   processes, used by the first interpolation to match fading instants 
%   with simulation instants. 
% * Variable Speed_ms, speed of the UE in m/s 
% * Variable CarrierFrequency_Hz, carrier frequency in Hz 
% * Variable RunningTimeInstant, time reference 
% * Variable ChipOversamplingFactor, oversampling factor of the main loop 
% * Variable ChipRate_Hz, chip rate 
% * Variable NumberOfChipsPerIteration, span of an iteration in chips 
% * Variable PDP_linear, PDP of the channel 
% * Variable NumberOfTxAntennas, number of antenna elements at Tx 
% * Variable NumberOfTransmittedSamples, time span of the main loop at Tx 
% * Variable NumberOfRxAntennas, number of antenna elements at Rx 
% * Variable NumberOfReceivedSamples, time span of the main loop at Rx, 
%   equals to NumberOfReceivedSamples plus the maximum delay expressed 
%   in samples (refer to IST METRA D3.2 for more details). 
% 
% Outputs 
% 
% * 2-D matrix Fading of size (NumberOfTxAntennas * NumberOfRxAntennas * 
%   NumberOfPaths) x (NumberOfChipsPerIteration * ChipOversamplingFactor) 
%   containing the fading samples of an iteration for all taps of the PDP. 
% * 3-D matrix Channel of size NumberOfTxAntennas x (NumberOfRxAntennas * 
%   NumberOfReceivedSamples) x NumberOfTransmittedSamples to be used to 
%   perform channel filtering according to the formalism described in 
%   IST METRA Deliverable D3.2, pp. 22-28. 
% 
% 
% STANDARD DISCLAIMER 
% 
% CSys is furnishing this item "as is". CSys does not provide any 
% warranty of the item whatsoever, whether express, implied, or 
% statutory, including, but not limited to, any warranty of 
% merchantability or fitness for a particular purpose or any 
% warranty that the contents of the item will be error-free. 
% 
% In no respect shall CSys incur any liability for any damages, 
% including, but limited to, direct, indirect, special, or 
% consequential damages arising out of, resulting from, or any way 
% connected to the use of the item, whether or not based upon 
% warranty, contract, tort, or otherwise; whether or not injury was 
% sustained by persons or property or otherwise; and whether or not 
% loss was sustained from, or arose out of, the results of, the 
% item, or any services that may be provided by CSys. 
% 
% (c) Laurent Schumacher, AAU-TKN/IES/KOM/CPK/CSys - February 2002 
 
% First interpolation 
% Interpolation in the fading vector domain 
 
TimeStep_s   = 1/(ChipOversamplingFactor*ChipRate_Hz); 
FadingStep_s = 3e8/(2*Speed_ms*CarrierFrequency_Hz*FadingOversamplingFactor); 
TimeInstants = RunningTimeInstant:TimeStep_s:RunningTimeInstant+ ... 
               ((NumberOfChipsPerIteration*ChipOversamplingFactor)-1)*TimeStep_s; 
 
% The mod(ulo) operator enables to loop on the fading vector along simulation. 
% This requires that the length of this fading vector is at least 40 wavelengths, 
% to be able to claim decorrelation between successive samples. 
FadingInstants         = mod(TimeInstants./FadingStep_s, size(FadingMatrixType,2)-1); 
FloorFadingInstants    = floor(FadingInstants); 
RemainerFadingInstants = FadingInstants - FloorFadingInstants; 
% Actual interpolation 
for ii = 1:(NumberOfChipsPerIteration*ChipOversamplingFactor) 
    Fading(:,ii) = (1-RemainerFadingInstants(ii))... 
                   .*FadingMatrixType(:,FloorFadingInstants(ii)+1)+... 
                   RemainerFadingInstants(ii)... 
                   .*FadingMatrixType(:,FloorFadingInstants(ii)+2); 
end; 
 
% Second interpolation 
% Interpolation in the tap domain 
% The delay values of the PDP do not necessarily match with the sampling instants 
% of the simulation. This second interpolation step derives the values of the taps 
% by linearly interpolating the correlated fadings. Note that the coefficients 
% of the linear interpolation are given as square root values of the distance 
% between samples, since the interpolation is linear in powers to maintain 
% the property that sum(PDP) = 1 
% 
% Note: this step could be optimised. Indeed, the interpolation in the tap domain 
% Could be performed out of the loop, since it always delivers the same result. 
% Tap delays and sampling instants do not change from iteration to the other. 
PDPInstants         = PDP_linear(2,:)./TimeStep_s; 
FloorPDPInstants    = floor(PDPInstants); 
RemainerPDPInstants = PDPInstants - FloorPDPInstants; 
Channel             = zeros(NumberOfTxAntennas,... 
	                  NumberOfRxAntennas*NumberOfReceivedSamples,... 
	                  NumberOfTransmittedSamples); 
% Actual interpolation 
for ii = 1:NumberOfTransmittedSamples 
    for jj = 1:size(PDP_linear,2) 
        shift = ((ii-1)*ChipOversamplingFactor) + FloorPDPInstants(jj) + 1; 
        for kk = 1:NumberOfTxAntennas 
            for ll = 1:NumberOfRxAntennas 
                Channel(kk, shift + (ll - 1) * NumberOfReceivedSamples, ii) = ... 
                    Channel(kk, shift + (ll - 1) * NumberOfReceivedSamples, ii) + ... 
                    sqrt(1 - RemainerPDPInstants(jj)) * Fading((kk - 1) * NumberOfRxAntennas + ll, ii); 
                Channel(kk, shift + (ll - 1) * NumberOfReceivedSamples + 1, ii) = ... 
                    Channel(kk, shift + (ll - 1) * NumberOfReceivedSamples + 1, ii) + ... 
                    sqrt(RemainerPDPInstants(jj)) * Fading((kk - 1) * NumberOfRxAntennas + ll, ii); 
            end; 
        end; 
    end; 
end; 
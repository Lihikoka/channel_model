% Plots the impulse responses, the PDPs, the correlation properties 
% and the Doppler spectra of the MIMO channel stored in the global 
% variable H produced by script example_MIMO.m. plot_MIMO.m also 
% relies on other global variables of example_MIMO.m 
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
 
close all; 
 
figure(1) 
colour = ['b','r','k','m','g','c']; 
abscissa = (1:size(H,4)).*NumberOfChipsPerIteration.*ChipOversamplingFactor./ChipRate_Hz; 
for ii = 1:NumberOfTxAntennas 
    for jj = 1:NumberOfRxAntennas 
        for kk = 1:NumberOfPaths 
            subplot(NumberOfTxAntennas,NumberOfRxAntennas,((ii-1)*NumberOfRxAntennas)+jj),... 
                semilogy(abscissa,abs(reshape(H(ii,jj,kk,:),1,NumberOfIterations)),colour(kk)), ... 
                hold on; 
        end; 
        grid; 
        xlabel('Time [s]'); 
        title(['Tx#',num2str(ii),' - Rx#',num2str(jj)]); 
    end; 
end; 
 
% Check PDP 
 
figure(2); 
for ii = 1:NumberOfTxAntennas 
    for jj = 1:NumberOfRxAntennas 
        pdp = zeros(1,NumberOfPaths); 
        for kk = 1:NumberOfPaths 
            pdp(kk) = sum((abs(reshape(H(ii,jj,kk,:),1,NumberOfIterations))).^2)./NumberOfIterations; 
        end; 
        if ((sum(pdp) < .9) | (sum(pdp) > 1.1)) 
            disp('Warning! sum(pdp) <> 1'); 
        end; 
        subplot(NumberOfTxAntennas,NumberOfRxAntennas,(ii-1)*NumberOfRxAntennas+jj),... 
            stem(1:NumberOfPaths,10*log10(PDP_linear(1,:)./PDP_linear(1,1)),'r--'), hold on, ... 
            stem(1:NumberOfPaths,10*log10(pdp./pdp(1)),'b'), grid,... 
            title(['Tx#',num2str(ii),' - Rx#',num2str(jj),' - Sum PDP = ',num2str(sum(pdp))]); 
        axe    = axis; 
        axe(1) = 1; 
        axis(axe); 
        if (jj==1) 
            subplot(NumberOfTxAntennas,NumberOfRxAntennas,(ii-1)*NumberOfRxAntennas+jj),... 
                ylabel('Power [dB]'); 
        end; 
        if (ii==NumberOfTxAntennas) 
            subplot(NumberOfTxAntennas,NumberOfRxAntennas,(ii-1)*NumberOfRxAntennas+jj),... 
                xlabel('Tap index'); 
        end; 
    end; 
end; 
 
% Check correlation 
 
figure(3) 
Rcomplex = zeros(NumberOfPaths,NumberOfTxAntennas*NumberOfRxAntennas,NumberOfTxAntennas*NumberOfRxAntennas); 
for ii=1:NumberOfPaths 
    for jj=1:NumberOfTxAntennas 
        for kk=1:NumberOfRxAntennas 
            ll = (jj-1)*NumberOfRxAntennas + kk; 
            for mm=1:NumberOfTxAntennas 
                for nn=1:NumberOfRxAntennas 
                    oo                 = (mm-1)*NumberOfRxAntennas + nn; 
                    ha                 = reshape(H(jj,kk,ii,:),NumberOfIterations,1); 
                    hb                 = reshape(H(mm,nn,ii,:),NumberOfIterations,1); 
                    temp               = corrcoef(ha,hb); 
                    Rcomplex(ii,ll,oo) = temp(1,2); 
                end; 
            end; 
        end; 
    end; 
end; 
for ii = 1:NumberOfPaths 
    index_low  = (ii-1)*NumberOfTxAntennas*NumberOfRxAntennas+1; 
    index_high = ii*NumberOfTxAntennas*NumberOfRxAntennas; 
    for jj = 1:NumberOfTxAntennas*NumberOfRxAntennas 
        for kk = 1:NumberOfTxAntennas*NumberOfRxAntennas 
            temp(jj,kk) = abs(Rcomplex(ii,jj,kk)).^2; 
        end 
    end 
    subplot(2,NumberOfPaths,ii),... 
        mesh(1:NumberOfTxAntennas*NumberOfRxAntennas,1:NumberOfTxAntennas*NumberOfRxAntennas,temp),... 
        title(['\langle h_{ij}^',num2str(ii),', h_{kl}^',num2str(ii),'\rangle']); 
    subplot(2,NumberOfPaths,NumberOfPaths+ii),... 
        plot(1:NumberOfTxAntennas*NumberOfRxAntennas, temp(1,:), 'b', ... 
             1:NumberOfTxAntennas*NumberOfRxAntennas, ... 
             abs(reshape(R(index_low:index_high,index_low), 1, ... 
             NumberOfTxAntennas*NumberOfRxAntennas)).^2, 'r--'),... 
             title(['\langle h_{11}^',num2str(ii),', h_{kl}^',num2str(ii),'\rangle']),... 
             xlabel(['(k-1)*',num2str(NumberOfRxAntennas),' + l']); 
    axe    = axis; 
    axe(1) = 1; 
    axe(3) = 0; 
    axis(axe); 
    grid; 
    if (ii == 1) 
        subplot(2,NumberOfPaths,ii), zlabel('Correlation coefficient'); 
        subplot(2,NumberOfPaths,NumberOfPaths+ii), ylabel('Correlation coefficient'); 
    end; 
end; 
 
% Check Doppler spectrum 
 
figure(4) 
normalised_frequency = (-NumberOfIterations/2+1:1:NumberOfIterations/2)... 
                       * (ChipRate_Hz * ChipOversamplingFactor) ... 
                       / (NumberOfChipsPerIteration * NumberOfIterations * Max_Doppler_shift); 
for ii = 1:NumberOfTxAntennas 
    for jj = 1:NumberOfRxAntennas 
        for kk = 1:NumberOfPaths 
            spectrum = (abs(fftshift(fft(reshape(H(ii, jj, kk, :), 1, NumberOfIterations))))./NumberOfIterations); 
            subplot(NumberOfTxAntennas*NumberOfRxAntennas,NumberOfPaths,... 
                (((ii-1)*NumberOfRxAntennas)+(jj-1))*NumberOfPaths+kk),plot(normalised_frequency,spectrum); 
            axe    = axis; 
            axe(1) = max(-2,axe(1)); 
            axe(2) = min(2,axe(2)); 
            axis(axe); 
            lower_bound = line(ones(1,2),[axe(3) axe(4)]); 
            set(lower_bound,'Color',[1 0 0],'LineStyle',':'); 
            upper_bound = line(-1*ones(1,2),[axe(3) axe(4)]); 
            set(upper_bound,'Color',[1 0 0],'LineStyle',':'); 
            title(['Tap h_{',num2str(ii),num2str(jj),'}^',num2str(kk)]); 
        end; 
    end; 
end;

%close all
clear all
clc

%% Initialization.
M = 2; % BPSK
N = 256; % number of carriers
CP = 32; % length of cyclic prefix
fs = 3.84e6; % Hz

n = 8; % number of subcarriers per group
K = [3 5]; % number of active subcarriers per group
g = N/n;
num_possibilities = zeros(1,length(K));
Zp_start = ones(1,length(K)+1);
for i = 1:length(K)
    num_possibilities(i) = M^K(i)*nchoosek(n,K(i));
    if i ~= length(K)
        Zp_start(i+1) =  Zp_start(i) + num_possibilities(i);
    end
end

SE = floor(log2(sum(num_possibilities)))/n;
SE_wCP = SE*N/(N+CP);

Zp_start(end) = 2^(SE*n)+1;

% Look-up tables
p1_range = zeros(length(K),2);
for i = 1:length(K)
    Zp_min_bin = de2bi(Zp_start(i),SE*n,'left-msb');
    Zp_min_forSAP = Zp_min_bin(1:end-K(i)*log2(M));
    p1_range(i,1) = bi2de(Zp_min_forSAP,'left-msb');
    
    Zp_max_bin = de2bi(Zp_start(i+1)-2,SE*n,'left-msb');
    Zp_max_forSAP = Zp_max_bin(1:end-K(i)*log2(M));
    p1_range(i,2) = bi2de(Zp_max_forSAP,'left-msb');
end
num_SAPs = p1_range(:,2) - p1_range(:,1) + 1;

LUT3 = nchoosek(1:8,K(1));
LUT5 = nchoosek(1:8,K(2));
LUT5 = LUT5(1:num_SAPs(end),:);

framen = 140; % 10ms frame -> framen is actually number of symbols for 1 frame and not the number of frames
repetitions = 100; % just here so i dont have to delete the indentations :)
total = (SE*N)*framen;

randomgenerationchannel = 1;
addnoise = 1; % change this to 1 before the simulation >:(
repetctr = 0;
%
interleaved = 1;
snrval = 0:1:5;
%snrval = 20;

berml = zeros(1,length(snrval));
%aveberml = berml+10e-10;
aveberml = berml;
aveserml = berml;
aveierml = berml;
total_indices = 0;

pathDelays_ETU = [0 50 120 200 230 500 1600 2300 5000]*1e-9; % ETU/LTE
avgPathGains_ETU = [-1 -1 -1 0 0 0 -3 -5 -7]; % dB
pathDelays_EPA = [0 30 70 90 110 190 410]*1e-9; % EPA/LTE
avgPathGains_EPA = [0 -1 -2 -3 -8 -17.2 -20.8]; % dB
pathDelays_EVA = [0 30 150 310 370 710 1090 1730 2510]*1e-9; % EVA/LTE
avgPathGains_EVA = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7 -12 -16.9]; % dB

PDP = 'ETU'; % change this based on factor levels
fD = 5;

tic
%% loop
for snrctr = 1:length(snrval)
    for rpt=1:repetitions
        %% Generating and coding data
        frame_size = SE*N;
        totalbits = frame_size*framen;

        %data=randi([0,1],totalbits,1)'; %bit stream
        data = zeros(1,totalbits);
            
        grouped = reshape(data,SE*n,g*framen); % will have g rows per frame
        grouped = transpose(grouped);
        dec_grouped=bi2de(grouped,'left-msb'); % will convert each bitstream
                                               % for each group into a
                                               % decimal value
        
        %% Mapping
        p1 = zeros(1,g);
        p2 = zeros(1,g);
        databins = zeros(framen,N);
        origactive = [];
        origactive_marker = [];
        num_K_perframe = zeros(1,framen);
        
        % the following loop maps each row (bitstream) to an SAP and k
        % M-ary constellation symbols
        for framectr = 1:framen
            for i = 1:g
                for j = 1:length(Zp_start)-1
                    if dec_grouped((framectr-1)*g+i) >= Zp_start(j)-1
                        if dec_grouped((framectr-1)*g+i) < Zp_start(j+1)-1
                            % grouping the bits
                            num_forBPSK = K(j)*log2(M);
                            num_forSAP = SE*n - num_forBPSK;
                            bits_forBPSK = grouped((framectr-1)*g+i,end-num_forBPSK+1:end);
                            bits_forSAP = grouped((framectr-1)*g+i,1:end-num_forBPSK);
                            dec_p1 = bi2de(bits_forSAP,'left-msb');

                            % LUT
                            if K(j) == 3
                                LUT = LUT3;
                            elseif K(j) == 5
                                LUT = LUT5;
                            end

                            % BPSK modulation
                            y_bpsk = pskmod(bits_forBPSK,M);
                            % should change this to a decimal value kung
                            % hindi bpsk! but since constrained to bpsk oks
                            % lang to

                            % Distributing into IFFT bins
                            norm_dec_p1 = dec_p1 - p1_range(j,1) + 1;
                            for i2 = 1:K(j)
                                ind = (i-1)*n+LUT(norm_dec_p1,i2);
                                origactive = [origactive ind];
                                
                                % the following if else statement marks the
                                % start of a new frame
                                if length(origactive)==1
                                    origactive_marker = [origactive_marker framectr];
                                elseif origactive(end)<origactive(end-1)
                                    origactive_marker = [origactive_marker framectr];
                                else
                                    origactive_marker = [origactive_marker 0];
                                end
                                
                                databins(framectr,ind) = y_bpsk(i2);
                            end
                            
                            num_K_perframe(framectr) = num_K_perframe(framectr)+K(j);
                        end
                    end
                end
            end
        end
        databins_beforeintlv = databins;
        
        
        %% INTERLEAVING 
        if interleaved == 1
            trmatind = 1:N;
            trmatindmat = reshape(trmatind,n,g);
            trmatind = reshape(transpose(trmatindmat),1,N);
            
            trmateye = eye(N);
            Imat = zeros(N);
            for trmatctr = 1:N
                Imat(trmatctr,:) = trmateye(trmatind(trmatctr),:);
            end
            
            for i = 1:framen
                databins_prime(:,i) = Imat*transpose(databins(i,:));
            end
            databins = transpose(databins_prime); % results in a framenxN matrix
        end
            
        
        %% IFFT
        %ifft_sig=ifft(databins_precoded,N);

        ifftmat = conj(dftmtx(N))/N;
        ifft_sig = zeros(N,framen);
        for i = 1:framen
        	ifft_sig(:,i) = ifftmat*transpose(databins(i,:));
        end
        ifft_sig = transpose(ifft_sig);
        
        %% CYCLIC PREFIX
        x_wCP = zeros(framen,N+CP);
        x_wCP(:,1:CP) = ifft_sig(:,end-CP+1:end);
        x_wCP(:,CP+1:end) = ifft_sig;
            
        xt_wCP = transpose(reshape(transpose(x_wCP),1,[]));
        
        %% CHANNEL
        
        if PDP == 'ETU'
            pathDelays = pathDelays_ETU;
            avgPathGains = avgPathGains_ETU;
        elseif PDP == 'EPA'
            pathDelays = pathDelays_EPA;
            avgPathGains = avgPathGains_EPA;
        elseif PDP == 'EVA'
            pathDelays = pathDelays_EVA;
            avgPathGains = avgPathGains_EVA;
        end
        
        rayleighchan = comm.RayleighChannel('SampleRate',fs, ...
        'PathDelays',pathDelays, ...
        'AveragePathGains',avgPathGains, ...
        'MaximumDopplerShift',fD,'PathGainsOutputPort',1,'NormalizePathGains',0);
    
        channel_filtcoeff = info(rayleighchan).ChannelFilterCoefficients;
        channel_filtdelay = info(rayleighchan).ChannelFilterDelay;


        [yt,pathgains] = rayleighchan(xt_wCP);
        % outputs are tall vectors of height framen*(N+CP)

        h = zeros(framen,size(channel_filtcoeff,2));
        frame_start = [0:framen-1].*(N+CP)+1;
        for framectr = 1:framen
            h_gains = mean(pathgains(frame_start(framectr):frame_start(framectr)+N+CP-1,:),1);
            for i = 1:length(h_gains)
                h(framectr,:) = h(framectr,:) + h_gains(i)*channel_filtcoeff(i,:);
            end
        end
        %h_wodelay = h(channel_filtdelay+1:end);

            
        ht_ext = zeros(framen,N);
        %ht_ext(:,1:size(h_gains,2)) = h_gains;
        ht_ext(:,1:size(h,2)) = h;
        % this is a framen by N matrix with h for each frame
            
        % fft applied to each column if X is a matrix
        htf = transpose(fft(transpose(ht_ext)));
        %htf = ones(1,N);
        
        %yt = Hmatrix*ifft_sig.'; % Hmatrix is an NxN matrix, ifft_sig is an Nx1 matrix
        %yt = transpose(xt_wCP); % no channel effect
        %Hmatrix = eye(N);

        %figure;plot(f,abs(fftshift(fft(ifft_sig))));
        %figure;plot(f,abs(fftshift(fft(yt))));
        
    
        %% NOISY CHANNEL
        if addnoise == 1
            ypluswt = awgn(yt,snrval(snrctr),'measured');
        else
            ypluswt = yt;
        end
        noise = ypluswt - yt;
        noisevar = var(noise);
        noisevarf = (g/N)*noisevar;
       

        %figure;plot(f,abs(fftshift(fft(ypluswt))));


        %% REMOVING CP

        %rcvdsig = ypluswt(CP+1:length(ypluswt)); % removing cyclic prefix
        rcvd_byframe = transpose(reshape(ypluswt,N+CP,framen));
        rcvd_woCP = rcvd_byframe(:,CP+1:end);
       	%yt_woCP = transpose(reshape(transpose(rcvd_woCP),1,[]));
        
        %% ML Detector
        %tic
        % deinterleaving matrix
        revtrmatind = 1:N;
        revtrmatindmat = reshape(revtrmatind,g,n);
        revtrmatind = reshape(transpose(revtrmatindmat),1,N);
            
        revtrmateye = eye(N);
        revImat = zeros(N);
        for revtrmatctr = 1:N
            revImat(revtrmatctr,:) = revtrmateye(revtrmatind(revtrmatctr),:);
        end      
        
        mln = zeros(1,N);
        grouplen = N*M/g;
        s = pskmod(0:M-1,M);
        sf = s;
        %yf = fft(ypluswt);
        yf = fft(transpose(rcvd_woCP));
        
        if interleaved == 1
            for i = 1:framen
                yf(:,i) = revImat*yf(:,i);
            end
        end
        

        % choose one posible set of indices and generate a M^2-length vector and
        % determine which is the minimum

        % generate a vector containing length(indeq) elements and determine which
        % is the minumum

        %possibleind = zeros(1,k);
        rcvdactive = [];
        rcvdactive_marker = [];
        rcvdactive_frame = [];
        xrcvdml_frame = zeros(framen,N);
        %possibleindcomb = permn(1:M,k); % all possible combinations of symbols of
                                        % QAM symbols for k active subcarriers
        num_combs = 2^(SE*n);
        bitsml =[];
        for framectr = 1:framen
            rcvdactive_frame = [];
            
            % h changes for each frame so the ff part gets Hmatrix for
            % the current frame
            Hmatrix = zeros(N,N);
            ht_frame = ht_ext(framectr,:);
            Hmatrix(:,1) = ht_frame.';

            for i=2:N
               	Hmatrix(:,i) = [ht_frame(N-i+2:N) ht_frame(1:N-i+1)].';
            end
            
            %Hmatrix = eye(N);
            if interleaved == 1
                A = revImat*dftmtx(N)*Hmatrix*ifftmat*Imat;
            else
                A = dftmtx(N)*Hmatrix*ifftmat;
            end
            % A*transpose(databins_beforeintlv) is comparable to yf so it
            % must be correct up to here
            
            for i=1:g
                minresid = zeros(1,sum(num_SAPs)); % contains all min residuals for each SAP
                minresidpos = zeros(sum(num_SAPs),max(K)); % contains the symbols corresponding to the min residual
                for i5=1:length(K) %for every k
                    if i5 == 1
                        LUT = LUT3;
                        k=K(1);
                    elseif i5 == 2
                        LUT = LUT5;
                        k=K(2);
                    end
                    minresidposdec = zeros(1,k); % contains the corresponding decimal values of the symbols
                    possibleindcomb = permn(1:M,k); % contains all possible combinations of symbols

                    for i3=1:size(LUT,1)
                        possiblesymb = zeros(1,M^k); % contains residual norm for each combination
                        residml = yf((i-1)*n+1:i*n,framectr); % gets current subblock from the current frame
                        Assample = zeros(n,k,M); %a matrix composed of column vectors
                        possibleind = LUT(i3,:) + (i-1)*n;
                        for i4=1:k %for every active index
                            for i1=1:M %for every possible symbol per subcarrier
                                Assample(:,i4,i1) = A((i-1)*n+1:i*n,possibleind(i4))*sf(i1); %Nx1 vector
                            end
                        end

                        % generate all possible combinations
                        for i2=1:M^k
                            residml = yf((i-1)*n+1:i*n,framectr);
                            for i6=1:k
                                residml = residml - Assample(:,i6,possibleindcomb(i2,i6));
                            end
                            possiblesymb(i2) = norm(residml);
                        end

                        % now that we've generated all possible combinations, get the minimum
                        % for this set of indices
                        [minposssymbval,minposssymbpos] = min(possiblesymb);
                        minresid(sum(num_SAPs(1:i5-1))+i3) = minposssymbval;
                        minresidposdec = possibleindcomb(minposssymbpos,:);
                        minresidpos(sum(num_SAPs(1:i5-1))+i3,1:k) = sf(minresidposdec);
                    end
                    %minresid % delete later
                end

                [val post] = min(minresid);
                if post <= sum(num_SAPs(1))
                    LUT = LUT3;
                    k=K(1);
                    post_norm = post;
                    start = p1_range(1,1);
                else
                    LUT = LUT5;
                    k=K(2);
                    post_norm = post - sum(num_SAPs(1));
                    start = p1_range(2,1);
                end

                rcvdactive = [rcvdactive LUT(post_norm,:)+(i-1)*n];
                rcvdactive_frame = [rcvdactive_frame LUT(post_norm,:)+(i-1)*n];
                rcvdactive_marker = [rcvdactive_marker framectr zeros(1,length(LUT(post_norm,:))-1)];
                xrcvdml_frame(framectr,LUT(post_norm,:) + (i-1)*n) = minresidpos(post,1:k);
                subblock = zeros(1,n);
                subblock(LUT(post_norm,:)) = minresidpos(post,1:k);
                num_bits = SE*n - k*log2(M);
                [bits_subblock,bitsind,bitsqam] = subblockrec(subblock,post_norm,LUT(post_norm,:),n,k,M,num_bits,start);

                bitsml = [bitsml bits_subblock];
            end
            % error counting for indices
            indctr = 1;
            origactive_frame_ind = [];
            while length(origactive_frame_ind) < 2 
                if origactive_marker(indctr) == framectr
                    origactive_frame_ind = [origactive_frame_ind indctr];
                    
                    if framectr==framen
                        origactive_frame_ind = [origactive_frame_ind length(origactive)];
                    end
                elseif framectr<framen && origactive_marker(indctr) == framectr+1
                    origactive_frame_ind = [origactive_frame_ind indctr-1];
                end
                indctr = indctr+1;
            end
            origactive_frame = origactive(origactive_frame_ind(1):origactive_frame_ind(2));
            mlinderror = numerror(origactive_frame,rcvdactive_frame,'IER')
            aveierml(snrctr) = aveierml(snrctr)+mlinderror
        end
        
        xrcvdml = reshape(transpose(xrcvdml_frame),1,[]);
        
        %rcvdactive = sort(rcvdactive);
        %timerML = toc;
        
        %%
        origactive;

        %% Error counting
        % QAM error counting
        databins_forcomp = reshape(transpose(databins_beforeintlv),1,[]);
        mlsyme = numerror(databins_forcomp,xrcvdml,'SER')

        % bitstream error counting
        mlbite = numerror(data,bitsml,'BER')
       
        aveserml(snrctr) = aveserml(snrctr)+mlsyme
        aveberml(snrctr) = aveberml(snrctr)+mlbite
        
        total_indices = total_indices+length(origactive);
        if rpt == repetitions
            aveierml(snrctr) = aveierml(snrctr)/total_indices;
            total_indices = 0;
        end
        
        repetctr = repetctr + 1
    end
end

aveberml = aveberml/total
aveserml = aveserml/(N*repetitions*framen)
aveierml
timesim = toc

save('aveber_snr0to5_wrayleigh.mat','aveberml')
save('aveser_snr0to5_wrayleigh.mat','aveserml')
save('aveier_snr0to5_wrayleigh.mat','aveierml')


%{
figure;
semilogy(snrval,aveberml,'ko-.')
grid on
title('BER vs SNR')
xlabel('SNR (dB)')
ylabel('Bit Error Rate (BER)')

figure;
semilogy(snrval,aveserml,'ko-.')
grid on
title('BER vs SNR')
xlabel('SNR (dB)')
ylabel('Symbol Error Rate (BER)')

figure;
semilogy(snrval,aveierml,'ko-.')
grid on
title('BER vs SNR')
xlabel('SNR (dB)')
ylabel('Index Error Rate (BER)')
%}
%save('run0.mat','aveberml');

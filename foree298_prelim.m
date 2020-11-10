%close all
clear all
clc

%% Initialization.
M = 2; % BPSK
N = 256; % number of carriers
CP = 16; % length of cyclic prefix

n = 8; % number of subcarriers per group
K = [1 2 3 4 5 6]; % number of active subcarriers per group
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

LUT1 = nchoosek(1:8,K(1));
LUT2 = nchoosek(1:8,K(2));
LUT3 = nchoosek(1:8,K(3));
LUT4 = nchoosek(1:8,K(4));
LUT5 = nchoosek(1:8,K(5));
LUT6 = nchoosek(1:8,K(6));
LUT6 = LUT6(1:num_SAPs(end),:);

repetitions = 3000;
total = (SE*N)*repetitions;

randomgenerationchannel = 1;
addnoise = 1;
repetctr = 0;
%
interleaved = 1;
snrval = 0:1:15;
%snrval = 0;

berml = zeros(1,length(snrval));
%aveberml = berml+10e-10;
aveberml = berml;
aveserml = berml;
aveierml = berml;
total_indices = 0;

tic
%% loop
for snrctr = 1:length(snrval)
    for rpt=1:repetitions
        %% Generating and coding data
        frame_size = SE*N;
        framen = 1;
        totalbits = frame_size*framen;

        data=randi([0,1],totalbits,1)'; %bit stream
        %data = ones(1,totalbits);

        s2=length(data);
        j_intlv=s2/(SE*n);
        intlvddata=reshape(data,j_intlv,SE*n);
        intlvddatastream = reshape(transpose(intlvddata),1,[]);
        
        %{
        data_2 = randi([0,1],SE*n,1)';
        intlvddatastream = [];
        for i=1:g
           intlvddatastream = [intlvddatastream data_2]; 
        end
        %}
            
        grouped = reshape(intlvddatastream,SE*n,g);
        grouped = transpose(grouped);
        dec_grouped=bi2de(grouped,'left-msb');
        
        %% Mapping
        p1 = zeros(1,g);
        p2 = zeros(1,g);
        databins = zeros(1,N);
        origactive = [];
        for i = 1:g
            for j = 1:length(Zp_start)-1
                if dec_grouped(i) >= Zp_start(j)-1
                    if dec_grouped(i) < Zp_start(j+1)-1
                        % grouping the bits
                        num_forBPSK = K(j)*log2(M);
                        num_forSAP = SE*n - num_forBPSK;
                        bits_forBPSK = grouped(i,end-num_forBPSK+1:end);
                        bits_forSAP = grouped(i,1:end-num_forBPSK);
                        dec_p1 = bi2de(bits_forSAP,'left-msb');
                        
                        % LUT
                        if K(j) == 1
                            LUT = LUT1;
                        elseif K(j) == 2
                            LUT = LUT2;
                        elseif K(j) == 3
                            LUT = LUT3;
                        elseif K(j) == 4
                            LUT = LUT4;
                        elseif K(j) == 5
                            LUT = LUT5;
                        else
                            LUT = LUT6;
                        end
                        
                        % BPSK modulation
                        y_bpsk = pskmod(bits_forBPSK,M);
                        
                        % Distributing into IFFT bins
                        norm_dec_p1 = dec_p1 - p1_range(K(j),1) + 1;
                        for i2 = 1:K(j)
                            ind = (i-1)*n+LUT(norm_dec_p1,i2);
                            origactive = [origactive ind];
                            databins(ind) = y_bpsk(i2);
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
            
            databins_prime = Imat*transpose(databins);
            databins = transpose(databins_prime);
        end
            
        
        %% IFFT
        %ifft_sig=ifft(databins_precoded,N);

        ifftmat = conj(dftmtx(N))/N;
        ifft_sig = ifftmat*transpose(databins);
        ifft_sig = transpose(ifft_sig);

        fs = 10e6;
        t = 0:1/fs:(N-1)*(1/fs);
        f = linspace(-fs/2,fs/2,N);

        %{
        %% CHANNEL
        % circular symmetric complex gaussian random variables CN(0,1/v)
        % distribution -- 10 taps
        taps = 10;
        variance = 1/taps;

        if randomgenerationchannel == 1
            ht = sqrt(variance/2)*(randn(1,taps)+1i*randn(1,taps));
        else
            htS = load('ht.mat');
            ht = htS.ht;
        end
            %ht = sqrt(variance/2)*(randi([5 10],1,taps)+1i*randi([5 10],1,taps))/10;
            %ht = 1:taps;


        ht_ext = zeros(1,N);
        ht_ext(1:taps) = ht;

        htf = fft(ht_ext);


        Hmatrix = zeros(N,N);
        Hmatrix(:,1) = ht_ext.';

        for i=2:N
            Hmatrix(:,i) = [ht_ext(N-i+2:N) ht_ext(1:N-i+1)].';
        end
        %}
        %yt = Hmatrix*ifft_sig.'; % Hmatrix is an NxN matrix, ifft_sig is an Nx1 matrix
        yt = ifft_sig; % no channel effect
        Hmatrix = eye(N);

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


        %% RECEIVER

        %rcvdsig = ypluswt(CP+1:length(ypluswt)); % removing cyclic prefix
        rcvddemod = fft(ypluswt,N);

        
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
         
        
        if interleaved == 1
            A = revImat*dftmtx(N)*Hmatrix*ifftmat*Imat;
        else
            A = dftmtx(N)*Hmatrix*ifftmat;
        end
        
        
        mln = zeros(1,N);
        grouplen = N*M/g;
        s = pskmod(0:M-1,M);
        sf = s;
        %yf = fft(ypluswt);
        yf = fft(ypluswt);
        
        if interleaved == 1
            yf = revImat*yf';
        end
        

        % choose one posible set of indices and generate a M^2-length vector and
        % determine which is the minimum

        % generate a vector containing length(indeq) elements and determine which
        % is the minumum

        %possibleind = zeros(1,k);
        rcvdactive = [];
        xrcvdml = zeros(1,N);
        %possibleindcomb = permn(1:M,k); % all possible combinations of symbols of
                                        % QAM symbols for k active subcarriers
        num_combs = 2^(SE*n);
        bitsml =[];

        for i=1:g
            minresid = zeros(1,sum(num_SAPs));
            minresidpos = zeros(sum(num_SAPs),max(K));
            for i5=1:length(K) %for every k
                if i5 == 1
                	LUT = LUT1;
                    k=K(1);
                elseif i5 == 2
                	LUT = LUT2;
                    k=K(2);
                elseif i5 == 3
                	LUT = LUT3;
                    k=K(3);
                elseif i5 == 4
                	LUT = LUT4;
                    k=K(4);
                elseif i5 == 5
                	LUT = LUT5;
                    k=K(5);
                else
                	LUT = LUT6;
                    k=K(6);
                end
                minresidposdec = zeros(1,k);
                possibleindcomb = permn(1:M,k);
                
                for i3=1:length(LUT)
                    possiblesymb = zeros(1,M^k); % contains residual norm for each combination
                    residml = yf((i-1)*n+1:i*n);
                    Assample = zeros(n,k,M); %a matrix composed of column vectors
                    possibleind = LUT(i3,:) + (i-1)*n;
                    for i4=1:k %for every active index
                        for i1=1:M %for every possible symbol per subcarrier
                            Assample(:,i4,i1) = A((i-1)*n+1:i*n,possibleind(i4))*sf(i1); %Nx1 vector
                        end
                    end

                    % generate all possible combinations
                    for i2=1:M^k
                        residml = yf((i-1)*n+1:i*n);
                        for i5=1:k
                            residml = residml - Assample(:,i5,possibleindcomb(i2,i5));
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
            end

            [val post] = min(minresid);
            if post <= sum(num_SAPs(1))
                LUT = LUT1;
                k=K(1);
                post_norm = post;
                start = p1_range(1,1);
            elseif post <= sum(num_SAPs(1:2))
                LUT = LUT2;
                k=K(2);
                post_norm = post - sum(num_SAPs(1));
                start = p1_range(2,1);
            elseif post <= sum(num_SAPs(1:3))
                LUT = LUT3;
                k=K(3);
                post_norm = post - sum(num_SAPs(1:2));
                start = p1_range(3,1);
            elseif post <= sum(num_SAPs(1:4))
                LUT = LUT4;
                k=K(4);
                post_norm = post - sum(num_SAPs(1:3));
                start = p1_range(4,1);
            elseif post <= sum(num_SAPs(1:5))
                LUT = LUT5;
                k=K(5);
                post_norm = post - sum(num_SAPs(1:4));
                start = p1_range(5,1);
            else
                LUT = LUT6;
                k=K(6);
                post_norm = post - sum(num_SAPs(1:5));
                start = p1_range(6,1);
            end
            
            
            rcvdactive = [rcvdactive LUT(post_norm,:)+(i-1)*n];
            xrcvdml(LUT(post_norm,:) + (i-1)*n) = minresidpos(post,1:k);
            subblock = zeros(1,n);
            subblock(LUT(post_norm,:)) = minresidpos(post,1:k);
            num_bits = SE*n - k*log2(M);
            [bits_subblock,bitsind,bitsqam] = subblockrec(subblock,post_norm,LUT(post_norm,:),n,k,M,num_bits,start);

            bitsml = [bitsml bits_subblock];
        end
        
        rcvdactive = sort(rcvdactive);
        %timerML = toc;
        
        %%
        origactive;

        %% Error counting
        mlerror = numerror(origactive,rcvdactive,'IER')
        
        % QAM error counting
        mlbe = numerror(databins_beforeintlv,xrcvdml,'SER')

        % bitstream error counting
        mlbite = numerror(intlvddatastream,bitsml,'BER')
       
        aveierml(snrctr) = aveierml(snrctr)+mlerror
        aveserml(snrctr) = aveserml(snrctr)+mlbe
        aveberml(snrctr) = aveberml(snrctr)+mlbite
        
        total_indices = total_indices + length(origactive);
        if rpt == repetitions
            aveierml(snrctr) = aveierml(snrctr)/total_indices;
            total_indices = 0;
        end
        
        repetctr = repetctr + 1
    end
end

aveberml = aveberml/total
aveserml = aveserml/(N*repetitions)
aveierml
timesim = toc
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

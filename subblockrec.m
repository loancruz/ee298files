function [datastream,bitsind,bitsqam] = streamrec(x,index,support,n,k,M,num_bits,start)
%%datastream = streamrec(x,support,lookup,n,k,g,M)
%x = all recovered symbols
%support = recovered support
%lutable = lookup table (row index corresponds to which bits are mapped to
%an index set
%n, k, g characterizing IM
%M characterizing QAM

% RECOVER IM BITS
bitsind = de2bi(index+start-1,num_bits,'left-msb');

% RECOVER QAM BITS
symbols = zeros(1,k);
symbolsbuffer = zeros(k,log2(M));
bitsqam = zeros(1,k*log2(M));
x_supp = x(support);
symbols = pskdemod(x_supp,M);
symbolsbuffer = de2bi(symbols,log2(M),'left-msb');
bitsqam = reshape(transpose(symbolsbuffer),1,[]);

% get bitstream
bitsmat = [bitsind bitsqam];
datastream = reshape(transpose(bitsmat),1,[]);

end
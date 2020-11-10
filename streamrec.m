function datastream = streamrec(x,support,lutable,n,k,g,M)
%%datastream = streamrec(x,support,lookup,n,k,g,M)
%x = all recovered symbols
%support = recovered support
%lutable = lookup table (row index corresponds to which bits are mapped to
%an index set
%n, k, g characterizing IM
%M characterizing QAM

% RECOVER IM BITS
p1 = floor(log2(nchoosek(n,k)));
bitsind = zeros(g,p1);
[possible,active] = size(lutable);

for i=1:g
    suppg = support((i-1)*k+1:i*k);    
    bitsind(i,:) = bitrec(mod(suppg-1,n)+1,lutable);
end

% RECOVER QAM BITS
symbols = zeros(g,k);
symbolsbuffer = zeros(k,log2(M));
bitsqam = zeros(g,k*log2(M));
x_supp = x(support);
for i=1:g
    symbols(i,:) = qamdemod(x_supp((i-1)*k+1:i*k),M);
    symbolsbuffer = de2bi(symbols(i,:),log2(M),'left-msb');
    bitsqam(i,:) = reshape(transpose(symbolsbuffer),1,[]);
end

% get bitstream
bitsmat = [bitsind bitsqam];
datastream = reshape(transpose(bitsmat),1,[]);

end
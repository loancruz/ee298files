function bits = bitrec(support,lookuptable)
%support = [1 4];
%lookup = [1,2;2,3;3,4;1,4];

[possible,active] = size(lookuptable);
support = sort(support);
lookupi = zeros(1,active);

for index = 1:possible
    % assume that both set of indices are in ascending order
    lookupi = lookuptable(index,:);
    if (lookupi - support == 0)
        break;
    end
end

% assume that support is in the lookup table
bits = de2bi(index-1,log2(possible),'left-msb');
end
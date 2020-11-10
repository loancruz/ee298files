function [number,errorv] = numerror(v1,v2,mode)

    number = 0;

    if mode == 'IER'
        for i = 1:length(v1)
            indic = 0;
            j = 1;
            while indic == 0 && j <= length(v2)
                if v1(i) == v2(j)
                    indic = 1;
                else
                    if j == length(v2)
                        number = number + 1;
                    end
                    j = j+1;
                end
            end
        end
    else
        errorv = v1-v2;
        for i = 1:length(v1)
            if abs(errorv(i)) > 1e-10
                number = number+1;
            end
        end
    end

end
% Magnitude Conversion 
function Mw = Ml2Mw(Ml)	
    c0 = 0.742958756;
    c1 = 1.020103103;
    c2 = 0.742974226;
    c3 = 0.11847815;
    c4 = 4.022121888;
    if Ml <= c4
        Mw = c0 * Ml + c1;
    else
        Mw = c0 * c4 + c1 + c2 * (Ml-c4) + c3 * power(Ml-c4,2);
    end
end


    
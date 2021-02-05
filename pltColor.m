function pltColor =pltColor(chNr)
    switch chNr
        case 1
            pltColor= [0.1,0.1,0.7];
        case 2
            pltColor= [0,1,1];
        case 3
            pltColor= [0.3, 0, 0.6];
        case 4
            pltColor= [0.1, 0.4, 0];
        otherwise
            pltColor= [0.0, 0, 0];
    end
end
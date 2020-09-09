function strategy = OperatorSelector(P)
    select = rand;
    if select < P(1)
        strategy=1;
    elseif select < P(1)+P(2)
        strategy=2;
    elseif select < P(1)+P(2)+P(3)
        strategy=3;
    else
        strategy=4;
    end
end


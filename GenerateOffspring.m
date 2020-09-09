function Offspring = GenerateOffspring(strategy,x1,x2,x3,x4,x5)
    switch(strategy)        
        case 1
            offspring = GAhalf([x1;x2]);
        case 2
            offspring = DE(x1,x2,x3);           
        otherwise
            offspring = CGDE(x1,x2,x3,x4,x5);
    end
    
    Offspring = INDIVIDUAL(offspring,x1);
end


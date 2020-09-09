function [Nprobability,Nquality] = CalcProbability(Lquality,performance1,performance2,p_min,alpha)
% 根据PM规则计算质量
    performance1(performance1<performance2)=performance2(performance1<performance2);
    Nquality = Lquality+alpha*(performance1-Lquality);
    Nprobability = ones(1,4)*p_min+ones(1,4)*(1-p_min*4).*(Nquality/sum(Nquality));    
end


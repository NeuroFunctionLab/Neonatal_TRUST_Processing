function Ind = GetInd(eTE,LCflag,repTime)
    % Label: LCflag=1, Control: LCflag =2   
    eTENum = 4;
    Ind = (repTime-1)*(eTENum*2) + (eTE-1)*2 + LCflag;
end
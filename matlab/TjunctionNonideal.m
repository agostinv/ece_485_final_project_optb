function capT = tjunctionCap(Z2, width1)
    capT = width1*((100/(tanh(Z2*0.0072)))+(0.64*Z2) - 261);

end

function L1 = tjuctionInduct1(width2, width1, h, Z1, epsilon)
    lw = tjuncLw(Z1, epsilon);
    w1 = width1/h;
    w2 = width2/h;
    L1 = -1*w2*h*(((w2)*((-0.016*w1)+0.064))+(0.016/w1))*lw;
end

function L2 = tjunctionInduct2(width2, width1, h, Z2, epsilon)
    lw = tjunctLw(Z2, epsilon);
    w1 = width1/h;
    w2 = width2/h;
    L2 = {w2*((0.12*w1)-0.47)+(0.195*w1)-0.357+(0.0283*sin((pi*w1)-(0.75*pi)}*lw;
end


function Lw = tjunctLw(Zo, epsilR)
    C = 3*10^8;
    Lw = Zo*sqrt(epsilR)/C;
end
freqCut = 2*10^9;
frequencies = (0:0.005:4)*10^9;
len45 = (3*10^8)/(8*freqCut);
%freq =2*10^9;
Zo = 50;
stub1 = 143.6;
stub1len = len45;
tline12 =76.71;
tline12len=len45;
stub2 = 31.2;
stub2len=len45;
tline23=91.9;
tline23len=len45;
stub3 = 25.3;
stub3len=len45;
tline34=91.9;
tline34len=len45;
stub4 = 31.2;
stub4len=len45;
tline45 =76.71;
tline45len=len45;
stub5 = 143.6;
stub5len =len45;
Zl = 50;


stub1ni = 138.505000;
stub1lenni = 0.00076200508;
tline12ni =22.627400;
tline12lenni=0.00076272898;
stub2ni = 51.835300;
stub2lenni=0.0085134196;
tline23ni=102.441000;
tline23lenni=0.0031803848;
stub3ni = 58.015500;
stub3lenni=0.0184534556;
tline34ni=tline23ni;
tline34lenni=tline23lenni;
stub4ni = stub2ni;
stub4lenni=stub2lenni;
tline45ni =tline12ni;
tline45lenni=tline12lenni;
stub5ni = stub1ni;
stub5lenni =stub1lenni;

%Inductor and Capacitor circuit (LC)
L1 = 5.329611747263777e-09;
C1 = 2.127914689331038e-12;
L2 = 8.618583238843864e-09;
C2 = 2.127914689331038e-12;
L3 = 5.329611747263777e-09;

%inductor and Capacitor circuit (CL)
C1b = 1.825249*10^(-12); %Farad
L1b = 5.45587118*10^(-9);
C2b = 3.143354017*10^(-12);
L2b = 5.45587118*10^(-9);
C3b = 1.825249*10^(-12);

%T junction Parameters
t0lW =121.215748*2.54e-5;
s1W =9.59538*2.54e-5;
t12W =400.541*2.54e-5;
s2W =122.003*2.54e-5;
t23W =27.9105*2.54e-5;
s3W =100.14*2.54e-5;
t34W =t23W;
s4W =s2W;
t45W =t12W;
s5W =s1W;
er =4.29;
h=0.0015748;
cT1 = tjunctionCap(stub1ni, s1W);
cT2 = tjunctionCap(stub2ni, s2W);
cT3 = tjunctionCap(stub3ni, s3W);
cT4 = tjunctionCap(stub4ni, s4W);
cT5 = tjunctionCap(stub5ni, s5W);

Ltline0 = tjunctionInduct1(s1W, t0lW, h, Zo, er);
Ltline12 = tjunctionInduct1(s1W, t12W, h, tline12ni, er);
Lstub1 = tjunctionInduct2(s1W, t0lW, h, stub1ni, er);

Ltline21 = tjunctionInduct1(s2W, t12W, h, tline12ni, er);
Ltline23 = tjunctionInduct1(s2W, t23W, h, tline23ni, er);
Lstub2 = tjunctionInduct2(s2W, t12W, h, stub2ni, er);

Ltline32 = tjunctionInduct1(s3W, t23W, h, tline23ni, er);
Ltline34 = tjunctionInduct1(s3W, t34W, h, tline34ni, er);
Lstub3 = tjunctionInduct2(s3W, t23W, h, stub3ni, er);

Ltline43 = tjunctionInduct1(s4W, t34W, h, tline34ni, er);
Ltline45 = tjunctionInduct1(s4W, t45W, h, tline45ni, er);
Lstub4 = tjunctionInduct2(s4W, t34W, h, stub4ni, er);

Ltline54 = tjunctionInduct1(s5W, t45W, h, tline45ni, er);
Ltline50 = tjunctionInduct1(s5W, t0lW, h, Zo, er);
Lstub5 = tjunctionInduct2(s5W, t45W, h, stub5ni, er);

S11=[];
S12=[];
S21=[];
S22=[];

S11ni=[];
S12ni=[];
S21ni=[];
S22ni=[];

S11nie=[];
S12nie=[];
S21nie=[];
S22nie=[];

S11LC = [];
S12LC = [];
S21LC = [];
S22LC = [];

S11CL = [];
S12CL = [];
S21CL = [];
S22CL = [];

for freq = frequencies
    %non ideal Tlines
    AbcdNonIdealS1 = inductorSer(Ltline0,freq)*capParallel(cT1,freq)*abcdParamNonIdealStub(stub1ni,freq,stub1len,Lstub1)*inductorSer(Ltline12,freq);
    AbcdNonIdealS2 = abcdParamSeries(tline12ni,freq,tline12len)*inductorSer(Ltline21,freq)*capParallel(cT2,freq)*abcdParamNonIdealStub(stub2ni,freq,stub2len, Lstub2)*inductorSer(Ltline23,freq);
    AbcdNonIdealS3 = abcdParamSeries(tline23ni,freq,tline23len)*inductorSer(Ltline32,freq)*capParallel(cT3,freq)*abcdParamNonIdealStub(stub3ni,freq,stub3len,Lstub3)*inductorSer(Ltline34,freq);
    AbcdNonIdealS4 = abcdParamSeries(tline34ni,freq,tline34len)*inductorSer(Ltline43,freq)*capParallel(cT4,freq)*abcdParamNonIdealStub(stub4ni,freq,stub4len,Lstub4)*inductorSer(Ltline45,freq);
    AbcdNonIdealS5 = abcdParamSeries(tline45ni,freq,tline12len)*inductorSer(Ltline54,freq)*capParallel(cT5,freq)*abcdParamNonIdealStub(stub5ni,freq,stub5len,Lstub5)*inductorSer(Ltline50,freq);
    AbcdNonIdeal = AbcdNonIdealS1*AbcdNonIdealS2*AbcdNonIdealS3*AbcdNonIdealS4*AbcdNonIdealS5;
    
    %ideal Tline
    AbcdIdeal = abcdParamStub(stub1,freq,stub1len)*abcdParamSeries(tline12,freq,tline12len)*abcdParamStub(stub2,freq,stub2len)*abcdParamSeries(tline23,freq,tline23len)*abcdParamStub(stub3,freq,stub3len)*abcdParamSeries(tline34,freq,tline34len)*abcdParamStub(stub4,freq,stub4len)*abcdParamSeries(tline45,freq,tline12len)*abcdParamStub(stub5,freq,stub5len);
    
    AbcdNonIdealEst = abcdParamStub(stub1ni,freq,stub1lenni)*abcdParamSeries(tline12ni,freq,tline12lenni)*abcdParamStub(stub2ni,freq,stub2lenni)*abcdParamSeries(tline23ni,freq,tline23lenni)*abcdParamStub(stub3ni,freq,stub3lenni)*abcdParamSeries(tline34ni,freq,tline34lenni)*abcdParamStub(stub4ni,freq,stub4lenni)*abcdParamSeries(tline45ni,freq,tline12lenni)*abcdParamStub(stub5ni,freq,stub5lenni);

    %LC circuit
    AbcdLC = inductorSer(L1, freq)*capParallel(C1,freq)*inductorSer(L2,freq)*capParallel(C2,freq)*inductorSer(L3,freq);
    
    %CL circuit
    AbcdCL = capParallel(C1b, freq)*inductorSer(L1b,freq)*capParallel(C2b,freq)*inductorSer(L2b, freq)*capParallel(C3b,freq);
    
    %conversion to S parameters
    sLCIdeal=ABCDtoS(AbcdLC(1,1),AbcdLC(1,2),AbcdLC(2,1),AbcdLC(2,2), Zo);
    sCLIdeal=ABCDtoS(AbcdCL(1,1),AbcdCL(1,2),AbcdCL(2,1),AbcdCL(2,2), Zo);
    sIdeal=ABCDtoS(AbcdIdeal(1,1),AbcdIdeal(1,2),AbcdIdeal(2,1),AbcdIdeal(2,2), Zo);
    sNonIdeal=ABCDtoS(AbcdNonIdeal(1,1),AbcdNonIdeal(1,2),AbcdNonIdeal(2,1),AbcdNonIdeal(2,2),Zo);
    sNonIdealEst = ABCDtoS(AbcdNonIdealEst(1,1), AbcdNonIdealEst(1,2), AbcdNonIdealEst(2,1), AbcdNonIdealEst(2,2),Zo);
    
    S11LC = [S11LC,sLCIdeal(1,1)];
    S12LC = [S12LC,sLCIdeal(1,2)];
    S21LC = [S21LC,sLCIdeal(2,1)];
    S22LC = [S22LC,sLCIdeal(2,2)];
    
    S11CL = [S11CL,sCLIdeal(1,1)];
    S12CL = [S12CL,sCLIdeal(1,2)];
    S21CL = [S21CL,sCLIdeal(2,1)];
    S22CL = [S22CL,sCLIdeal(2,2)];
    
    S11 = [S11,sIdeal(1,1)]; 
    S12 = [S12,sIdeal(1,2)];
    S21 = [S21,sIdeal(2,1)];
    S22 = [S22,sIdeal(2,2)];
    
    S11ni = [S11ni, sNonIdeal(1,1)];
    S12ni = [S12ni, sNonIdeal(1,2)];
    S21ni = [S21ni, sNonIdeal(2,1)];
    S22ni = [S22ni, sNonIdeal(2,2)];
    
    S11nie = [S11nie, sNonIdealEst(1,1)];
    S12nie = [S12nie, sNonIdealEst(1,2)];
    S21nie = [S21nie, sNonIdealEst(2,1)];
    S22nie = [S22nie, sNonIdealEst(2,2)];
end

S12db = 20*log10(abs(S12ni));

S12niedb =20*log10(abs(S12nie));
plot(frequencies, S12niedb);
function [ABCD] = abcdParamStub(zstub,freq,len)
    lambda = (3*10^8)/freq;
    beta = 2*pi/lambda;
    Zin =  -1i*zstub*cot(beta*len);
    Yin =1/Zin;
    ABCD = [1,0;Yin,1];
end

function [ABCD] =abcdParamNonIdealStub(zstub,freq,len,L)
    lambda = (3*10^8)/freq;
    beta = 2*pi/lambda;
    Zin =  -1i*zstub*cot(beta*len);
    Znew = (1i*freq*L)+Zin;
    Yin = 1/Znew;
    ABCD = [1,0;Yin,1];
end

function [ABCD] =abcdParamSeries(Zser, freq, len)
    lambda = (3*10^8)/freq;
    beta = 2*pi/lambda;
    ABCD = [cos(beta*len),1i*Zser*sin(beta*len);1i*(1/Zser)*sin(beta*len),cos(beta*len)];
end

function [ABCD] =capParallel(C, freq)
    Yin = 1i*freq*2*pi*C;
    ABCD = [1,0;Yin,1];
end

function [ABCD] = inductorSer(L, freq)
    Z = 1i*freq*2*pi*L;
    ABCD=[1,Z;0,1];
end

function [Sparam] = ABCDtoS(A,B,C,D, Zo)
 S11 = (A + (B/Zo)-(C*Zo)-D)/(A+(B/Zo)+(C*Zo)+D);
 S12 = 2*(A*D - B*C)/(A+(B/Zo)+(C*Zo)+D);
 S21 = 2/(A+(B/Zo)+(C*Zo)+D);
 S22 = (-A +(B/Zo)-(C*Zo)+D)/(A+(B/Zo)+(C*Zo)+D);
 Sparam = [S11,S12;S21,S22];
% function Ad = ADistrib(lambda, l)
%     beta = 2*pi/lambda;
%     Aser = cos(beta*l);
% end
% 
% function Bd = BDistrib(lambda

end


%non Ideal T junctions
function capT = tjunctionCap(Z2, width1)
    capT = width1*((100/(tanh(Z2*0.0072)))+(0.64*Z2) - 261);

end

function L1 = tjunctionInduct1(width2, width1, h, Z1, epsilon)
    lw = tjunctLw(Z1, epsilon);
    w1 = width1/h;
    w2 = width2/h;
    L1 = -1*w2*h*(((w2)*((-0.016*w1)+0.064))+(0.016/w1))*lw;
end

function L2 = tjunctionInduct2(width2, width1, h, Z2, epsilon)
    lw = tjunctLw(Z2, epsilon);
    w1 = width1/h;
    w2 = width2/h;
    L2 = (w2*((0.12*w1)-0.47)+(0.195*w1)-0.357+(0.0283*sin((pi*w1)-(0.75*pi))))*lw;
end


function Lw = tjunctLw(Zo, epsilR)
    C = 3*10^8;
    Lw = Zo*sqrt(epsilR)/C;
end
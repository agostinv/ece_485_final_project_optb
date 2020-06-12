
S21data =table2array(T4S21);
Freq = S21data(:,1);
S21Dat = 20*log10(abs(S21data(:,3)));
plot(Freq,S21Dat);
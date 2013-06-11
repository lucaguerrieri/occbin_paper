kss = ((1/BETA-1+DELTAK)/ALPHA)^(1/(ALPHA-1));
css = -DELTAK*kss +kss^ALPHA;
iss = DELTAK*kss;
uss = (css^(1-GAMMAC)-1)/(1-GAMMAC);
vss = uss/(1-BETA);

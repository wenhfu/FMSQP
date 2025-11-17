function c = func(x)
c = [-x;x-1];
h = funh(x);
c = [c;h;-h];

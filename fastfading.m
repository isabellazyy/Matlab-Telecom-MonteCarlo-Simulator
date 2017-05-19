% fast fading

c = rayleighchan(1/10000,100);
lin = 10.^(RSSs/10);

RSSs = reshape(RSSs,[1 15]);


plot(RSSs)
lin = reshape(lin,[1,15]);
s = 0;
for i = 1 : 3
   
    out = filter(c,lin);
    power = 10*log10(abs(out))
    s = power + s;
    hold on
    plot(power)

end
hold on
plot(s)
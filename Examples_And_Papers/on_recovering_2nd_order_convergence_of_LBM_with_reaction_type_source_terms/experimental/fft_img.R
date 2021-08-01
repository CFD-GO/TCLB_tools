n = 100
L = 2
f = matrix(0,n,n)
x = L*(row(f)-0.5)/n
y = L*(col(f)-0.5)/n

f = ifelse((x-1)^2+(y-1)^2<0.5^2,1,0)
image(f)

fm = fft(f)

k1 = ((row(f)-1+n/2) %% n-n/2)
k2 = ((col(f)-1+n/2) %% n-n/2)
k1 = k1*2*pi/L
k2 = k2*2*pi/L

t = 0.01
dif = 1
fmt = fm * exp(-t * dif * (k1^2 + k2^2))
ft = fft(fmt,inverse = TRUE)

image(Re(ft))

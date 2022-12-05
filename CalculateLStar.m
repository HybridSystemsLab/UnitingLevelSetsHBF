function Lstar = CalculateLStar()
global z1Star I

z1Star = zeros(100,1);

I = eye(100);

Lstar = 0.25.*((z1Star - z1Star).'*I*(z1Star - z1Star));
end
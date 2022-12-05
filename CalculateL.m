function L = CalculateL(z1)
global z1Star I

L = 0.25.*((z1 - z1Star).'*I*(z1 - z1Star));

end
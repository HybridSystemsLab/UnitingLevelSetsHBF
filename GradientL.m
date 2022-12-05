function GradL = GradientL(z1)
global z1Star I

GradL = (2*0.25).*(I + I.')*(z1 - z1Star);

end
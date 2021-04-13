within DoubleSkinFacades.Media.Materials;
record Glass "Glass(lambda=1,rho=2200,cp=795)"
  extends BaseGlass(lambda = 0.8, rho = 2200, cp = 795, epsilon = 0.84, epsilon_glazing = 0.037);
end Glass;

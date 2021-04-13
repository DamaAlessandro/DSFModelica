within DoubleSkinFacades.Media.Gases;
record Air
  extends BaseGas(
    cp_a=1002.737,
    cp_b=1.2324*10^(-2),
    mu_a=3.723*10^(-6),
    mu_b=4.94*10^(-8),
    lambda_a=2.873*10^(-3),
    lambda_b=7.76*10^(-5),
    Mm=28.97);
end Air;

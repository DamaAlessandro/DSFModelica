within DoubleSkinFacades.Media.Gases;
record Argon
  extends BaseGas(
    cp_a= 521.9285,
    cp_b=0,
    mu_a=3.379*10^(-6),
    mu_b=6.451*10^(-8),
    lambda_a=2.285*10^(-3),
    lambda_b=5.149*10^(-5),
    Mm=39.948);
end Argon;

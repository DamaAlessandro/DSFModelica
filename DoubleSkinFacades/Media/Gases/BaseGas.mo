within DoubleSkinFacades.Media.Gases;
record BaseGas
  extends Modelica.Icons.Record;
  parameter Modelica.SIunits.PerUnit cp_a "Specific heat capacity coefficient";
  parameter Modelica.SIunits.PerUnit cp_b "Specific heat capacity coefficient";
  parameter Modelica.SIunits.PerUnit mu_a "Dynamic viscosity coefficient";
  parameter Modelica.SIunits.PerUnit mu_b "Dynamic viscosity coefficient";
  parameter Modelica.SIunits.PerUnit lambda_a "Thermal conductivity coefficient";
  parameter Modelica.SIunits.PerUnit lambda_b "Thermal conductivity coefficient";
  parameter Modelica.SIunits.MolarMass Mm "Molecular mass";
  final parameter Modelica.SIunits.SpecificEntropy R = 8314 / Mm "Specific ideal gas constant";

end BaseGas;

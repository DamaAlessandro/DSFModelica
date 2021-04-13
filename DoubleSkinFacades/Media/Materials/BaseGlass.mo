within DoubleSkinFacades.Media.Materials;
record BaseGlass
  extends Modelica.Icons.Record;
  parameter Modelica.SIunits.ThermalConductivity lambda "Thermal conductivity";
  parameter Modelica.SIunits.SpecificHeatCapacity cp "Specific heat capacity";
  parameter Modelica.SIunits.Density rho "Mass density";
  parameter Modelica.SIunits.PerUnit epsilon "Emissivity";
  parameter Modelica.SIunits.PerUnit epsilon_glazing "Emissivity of the glazing";

end BaseGlass;

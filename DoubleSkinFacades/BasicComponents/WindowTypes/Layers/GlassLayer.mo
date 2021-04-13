within DoubleSkinFacades.BasicComponents.WindowTypes.Layers;
model GlassLayer
  outer Modelica.SIunits.Irradiance Gdir "Direct irradiation on the façade";
  outer  Modelica.SIunits.Irradiance Gdiff "Diffuse irradiation on the façade";
  parameter Modelica.SIunits.Length s_glass "Thickness of the glass layer";
  parameter Modelica.SIunits.Height H "Height of the window";
  parameter Modelica.SIunits.Height L "Length of the window";
  parameter Media.Materials.BaseGlass material "Type of glass" annotation(choicesAllMatching = true);
  parameter Modelica.SIunits.Density rho = material.rho "Density of the glass";
  parameter Modelica.SIunits.SpecificHeatCapacity cp = material.cp "Specific heat capacity of the glass";
  parameter Modelica.SIunits.ThermalConductivity lambda = material.lambda "Heat transfer coefficient of the glass";
  parameter Modelica.SIunits.Temperature Tstart = 273.15 + 18 "Start temperature of the glass probar Tsurf_a Tsurf_b?";
  parameter Modelica.SIunits.PerUnit alpha_diff "Diffuse absortivity of the glass";
  final parameter Modelica.SIunits.Area A = H * L "Area of the window surface";
  final parameter Modelica.SIunits.Mass M = rho * L * H * s_glass "Mass of the glass layer";

  // Variables
  Modelica.SIunits.PerUnit alpha_dir;
  Modelica.SIunits.PerUnit shade_factor;
  Modelica.SIunits.Temperature T_glass(start = Tstart) "Mean internal temperature in the glass";
  Modelica.SIunits.Temperature Tsurf_a "Surface temperature of glass a";
  Modelica.SIunits.Temperature Tsurf_b "Surface temperature of glass b";
  Modelica.SIunits.HeatFlowRate Qa "Power into glassSurface_a [W]";
  Modelica.SIunits.HeatFlowRate Qb "Power to glassSurface_b [W]";
  Modelica.SIunits.HeatFlowRate Qint "Power stored in the glass [W]";
  Modelica.SIunits.HeatFlowRate Qabs "Power absorbed from solar radiation [W]";
 // Modelica.SIunits.HeatFlux qa "Heat flux into glassSurface_a [W/m2]";
  //Modelica.SIunits.HeatFlux qb "Heat flux to glassSurface_b [W/m2]";
  //Modelica.SIunits.HeatFlux qabs "Heat flux absorbed from solar radiation [W/m2]";
  final parameter Modelica.SIunits.ThermalResistance Rcond = 1*(s_glass/2)/ (lambda * A) "Conductive thermal resistance";

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a glassSurface_a
    annotation (Placement(transformation(extent={{-52,-10},{-32,10}}),
        iconTransformation(extent={{-52,-10},{-32,10}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b glassSurface_b
    annotation (Placement(transformation(extent={{46,-8},{66,12}}),
        iconTransformation(extent={{32,-10},{52,10}})));

equation

  // Constitutive equations
  Qa + Qabs - Qint - Qb = 0;
  Qa = 1/Rcond * (Tsurf_a - T_glass);
  Qb = 1/Rcond * (T_glass - Tsurf_b);
  Qint = cp * M * der(T_glass);
  Qabs = A*(Gdir * alpha_dir * (1 - shade_factor) + Gdiff * alpha_diff*(1-shade_factor));

  // Boundary variables
  Qa = glassSurface_a.Q_flow;
  Qb = -glassSurface_b.Q_flow;
  Tsurf_a = glassSurface_a.T;
  Tsurf_b = glassSurface_b.T;

  annotation (Icon(graphics={Rectangle(
          extent={{-40,80},{40,-80}},
          lineColor={28,108,200},
          fillColor={0,255,255},
          fillPattern=FillPattern.Solid)}));
end GlassLayer;

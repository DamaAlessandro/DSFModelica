within DoubleSkinFacades.BoundaryHeatTransfer;
model ExternalSurface

  final constant Real sigma(final unit = "W/(m2.K4)") = 5.670400e-8 "Stefan-Boltzmann constant";
  parameter Modelica.SIunits.PerUnit epsilon_ext = 0.84  "Exterior glass emissitivity";
  final constant Real pi = 2*Modelica.Math.asin(1.0);
  parameter Modelica.SIunits.Length H "Height of window";
  parameter Modelica.SIunits.Length L "Width of the window";
  final parameter Modelica.SIunits.Area A = H*L "Area of the window";
  outer Modelica.SIunits.Velocity v_wind;
  outer Modelica.SIunits.Angle theta_wind;
  parameter Modelica.SIunits.Temperature T_ground = 273.15 + 16 "Average temperature of the ground for sky temperature calculation";
  parameter Real orientation(final unit="deg") "Façade orientation angle in degrees" annotation (Dialog(group="Façade specifications"));
  Modelica.SIunits.Temperature Trad_ext "Radiative temperature of the exterior";
  Modelica.SIunits.Temperature T_ext "Radiative temperature of the exterior";
  Modelica.SIunits.Temperature Tsurf "Surface temperature of the exterior glass";
  Modelica.SIunits.Temperature T_sky "Temperature of the sky";
  Modelica.SIunits.Temperature Tm "Mean radiant temperature";
  Modelica.SIunits.CoefficientOfHeatTransfer hcv_ext "Convective heat transfer coefficient exterior";
  Modelica.SIunits.CoefficientOfHeatTransfer hr_ext "Radiative heat transfer coefficient exterior";
  Modelica.SIunits.CoefficientOfHeatTransfer hcr_ext "Convective-radiative heat transfer coefficient exterior";
  Modelica.SIunits.PerUnit vs1 "Vs1 for Swinbank calcultion";
  Modelica.SIunits.PerUnit vs2 "Vs2 for Swinbank calculation";
  //Modelica.SIunits.PerUnit vs "Vs for Swinbank calculation";
  Modelica.SIunits.Power Qsurf "Power from external surface";
  Modelica.SIunits.Power Qair "Power released to the atmosphere";
  Real cosa;

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a airExterior annotation (
      Placement(transformation(extent={{-110,-10},{-90,10}}),
        iconTransformation(extent={{-110,-10},{-90,10}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b glassSurface annotation (
      Placement(transformation(extent={{90,-10},{110,10}}), iconTransformation(
          extent={{90,-10},{110,10}})));

equation

  //Radiation
  T_sky = 0.05532 * T_ext^(1.5) "Swinbank model (Stanzel 1994)";
  Trad_ext = (T_sky + T_ground)/2;
  Tm = (Tsurf + Trad_ext)/2;
  hr_ext = 4 * sigma * epsilon_ext * Tm^3;
  cosa = (cos(theta_wind -orientation*Modelica.Constants.pi/180- Modelica.Constants.pi));

  //Convection
  vs1 = 1.8 - 1.6 * (cos(theta_wind -orientation*Modelica.Constants.pi/180- Modelica.Constants.pi));
  vs2 = 0.2 + 1.5 * (cos(theta_wind -orientation*Modelica.Constants.pi/180- Modelica.Constants.pi));
  hcv_ext = v_wind * vs1 + vs2;
  //hcv_ext = (1.7 * vs + 5.1);
  hcr_ext = hcv_ext + hr_ext;

  Qair = -Qsurf;
  Qsurf = hcr_ext * A * (Tsurf - T_ext);

  Tsurf = glassSurface.T;
  T_ext = airExterior.T;
  Qair = airExterior.Q_flow;
  Qsurf = glassSurface.Q_flow;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{20,80},{52,-80}},
          lineColor={28,108,200},
          fillColor={216,255,251},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-22,44},{4,38}},
          lineColor={238,46,47},
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-16,52},{-48,40},{-16,30},{-22,40},{-16,52}},
          lineColor={238,46,47},
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-20,2},{6,-4}},
          lineColor={238,46,47},
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-14,10},{-46,-2},{-14,-12},{-20,-2},{-14,10}},
          lineColor={238,46,47},
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-10,-36},{-42,-48},{-10,-58},{-16,-48},{-10,-36}},
          lineColor={238,46,47},
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-16,-44},{10,-50}},
          lineColor={238,46,47},
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid)}),                      Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end ExternalSurface;

within DoubleSkinFacades.BoundaryHeatTransfer;
model InternalSurface

  final constant Real sigma(final unit = "W/(m2.K4)") = 5.670400e-8 "Stefan-Boltzmann constant";
  parameter Modelica.SIunits.PerUnit epsilon_int = 0.8 "Interior glass emissitivity";
  parameter Modelica.SIunits.Length H "Height of window";
  parameter Modelica.SIunits.Length L "Width of the window";
  final parameter Modelica.SIunits.Area A = H*L "Area of the window";
  parameter Modelica.SIunits.Velocity v_int= 0.2 "Velocity of air in the interior";
  parameter Modelica.SIunits.PerUnit alpha = 1.9*10^(-5);
  parameter Modelica.SIunits.SpecificHeatCapacity cp= 1010;
  parameter Modelica.SIunits.Acceleration g = 9.806 "Gravity acceleration";

  Modelica.SIunits.Temperature Tint "Inside temperature";
  Modelica.SIunits.Temperature Tsurf "Surface temperature of the interior glass";
  Modelica.SIunits.Temperature Tmf "Mean radiant temperature";

  Modelica.SIunits.CoefficientOfHeatTransfer hcr_int "Radiative heat transfer coefficient interior";

  Modelica.SIunits.CoefficientOfHeatTransfer hr_int "Radiative heat transfer coefficient interior";
  Modelica.SIunits.Power Qsurf "Power from external surface";
  Modelica.SIunits.Power Qair "Power released to the atmosphere";

  Modelica.SIunits.CoefficientOfHeatTransfer hcv_int "Convective heat transfer coefficient interior";
  outer Modelica.SIunits.Pressure Patm;

  Modelica.SIunits.Temperature deltaT;

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b airInterior annotation (
      Placement(transformation(extent={{90,-10},{110,10}}),
        iconTransformation(extent={{90,-10},{110,10}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a glassSurface annotation (
      Placement(transformation(extent={{-110,-10},{-90,10}}),
                                                            iconTransformation(
          extent={{-110,-10},{-90,10}})));

equation
  //Radiation
  hr_int = 4 * sigma * epsilon_int * Tmf^3;

  if abs(Tsurf-Tint)<0.01 then deltaT=0.01;
  else deltaT=abs(Tsurf-Tint);
  end if;

  Tmf = Tint + 0.25*(Tsurf - Tint);

  hcv_int = ((1.5*(deltaT/H)^(0.25)) + (1.23*(deltaT)^(0.33))^6)^(1/6);

  hcr_int = hcv_int + hr_int;
  Qair = -Qsurf;
  Qair = hcr_int * A * (Tint - Tsurf);

  Tsurf = glassSurface.T;
  Tint = airInterior.T;
  Qair = airInterior.Q_flow;
  Qsurf = glassSurface.Q_flow;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-20,80},{-52,-80}},
          lineColor={28,108,200},
          fillColor={216,255,251},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{36,50},{10,44}},
          lineColor={238,46,47},
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{30,58},{62,46},{30,36},{36,46},{30,58}},
          lineColor={238,46,47},
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{36,2},{10,-4}},
          lineColor={238,46,47},
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{30,10},{62,-2},{30,-12},{36,-2},{30,10}},
          lineColor={238,46,47},
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{32,-38},{64,-50},{32,-60},{38,-50},{32,-38}},
          lineColor={238,46,47},
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{38,-46},{12,-52}},
          lineColor={238,46,47},
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid)}),                      Diagram(
        coordinateSystem(preserveAspectRatio=false)),
              Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end InternalSurface;

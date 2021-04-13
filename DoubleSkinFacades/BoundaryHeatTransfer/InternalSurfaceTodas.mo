within DoubleSkinFacades.BoundaryHeatTransfer;
model InternalSurfaceTodas

  final constant Real sigma(final unit = "W/(m2.K4)") = 5.670400e-8 "Stefan-Boltzmann constant";
  parameter Modelica.SIunits.PerUnit epsilon_int = 0.8 "Interior glass emissitivity";
  parameter Modelica.SIunits.Length H "Height of window";
  parameter Modelica.SIunits.Length L "Width of the window";
  final parameter Modelica.SIunits.Area A = H*L "Area of the window";
  parameter Modelica.SIunits.Velocity v_int= 0.2 "Velocity of air in the interior";
  parameter Modelica.SIunits.PerUnit nu = 1.8*10^(-5);
  parameter Modelica.SIunits.PerUnit alpha = 1.9*10^(-5);
  parameter Modelica.SIunits.SpecificHeatCapacity cp= 1010;
  parameter Modelica.SIunits.Acceleration g = 9.806 "Gravity acceleration";

  Modelica.SIunits.Temperature Tint "Inside temperature";
  Modelica.SIunits.Temperature Tsurf "Surface temperature of the interior glass";
  Modelica.SIunits.Temperature Tm "Mean radiant temperature";
  Modelica.SIunits.CoefficientOfHeatTransfer hcv_int "Convective heat transfer coefficient interior";
  Modelica.SIunits.CoefficientOfHeatTransfer hcv_int2 "Convective heat transfer coefficient interior";
  Modelica.SIunits.CoefficientOfHeatTransfer hcv_int3 "Convective heat transfer coefficient interior";
  Modelica.SIunits.CoefficientOfHeatTransfer hcr_int2 "Convective-radiative heat transfer coefficient interior";

  Modelica.SIunits.CoefficientOfHeatTransfer hr_int "Radiative heat transfer coefficient interior";
  //Modelica.SIunits.CoefficientOfHeatTransfer hcr_int "Convective-radiative heat transfer coefficient interior";
  Modelica.SIunits.Power Qsurf "Power from external surface";
  Modelica.SIunits.Power Qair "Power released to the atmosphere";
  Modelica.SIunits.PerUnit lambda;
  Modelica.SIunits.PerUnit Gr(start=1) "Rayleigh number";
  Modelica.SIunits.PerUnit Nu "Nusselt number";

    Modelica.SIunits.Temperature Tmf "Inside temperature";
  Modelica.SIunits.PerUnit Nu4 "Nusselt number";
   Modelica.SIunits.PerUnit nu4;
   Modelica.SIunits.PerUnit Ra_h;
   Modelica.SIunits.PerUnit Ra_cv;
  Modelica.SIunits.Density rho;
     Modelica.SIunits.PerUnit lambda4;
  Modelica.SIunits.CoefficientOfHeatTransfer hcv_int4 "Convective heat transfer coefficient interior";
  outer Modelica.SIunits.Pressure Patm;

      Modelica.SIunits.PerUnit Pr;
  Modelica.SIunits.PerUnit Nu3;
Modelica.SIunits.PerUnit Ra_1(start=1);

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
  Tm = (Tsurf + Tint)/2;
  hr_int = 4 * sigma * epsilon_int * Tm^3;

  //Convection
  lambda = (2.68482*10^(-3) * Tint^(1.5)) / (Tint + 254.4*10^(-12/Tint));
  Ra_1 = abs((g * 1/Tint * deltaT * H^3)/(nu * alpha));
  Nu3 = (0.825 + 0.325 * (abs(Ra_1)^(1/6)))^2;
  Nu = 0.478 * (abs(Gr)^0.25);
  Gr = (g * H^3 * (Tsurf-Tint))/(nu^2 * Tint);
  hcv_int =  Nu * lambda/H;
  hcv_int3 = Nu3 * lambda /H;
  Pr = 1010*nu/0.026;
  hcv_int2 = 0.026 / H * (0.825 + 0.387*((abs(Ra_1))^(1/6))/((1+((0.492/Pr)^(9/16)))^(4/9)));

  if abs(Tsurf-Tint)<0.01 then deltaT=0.01;
  else deltaT=Tsurf-Tint;
  end if;

  //4th correlation
  rho = Patm/(286.9*Tint);
  Tmf = Tint + 0.25*(Tsurf - Tint);
  lambda4 = 2.873*10^(-3) + 7.76*10^(-5)*Tmf;
  nu4 = 3.723*10^(-6) + 4.94*10^(-8)*Tmf;
  Ra_h = rho^2 * H^3 * g * cp * abs(deltaT)/(Tmf* nu4 * lambda4);
  Ra_cv = 2.5 * 10^5 * (exp(0.72*90))^0.2;

  if Ra_h <= Ra_cv then Nu4 = 0.56 * (Ra_h)^0.25;
  else Nu4 = 0.13 * (Ra_h^(1/3) - Ra_cv^(1/3)) + 0.56 * (Ra_cv)^0.25;
  end if;

  hcv_int4 = Nu4 * lambda4 / H;

  //hcr_int = hcv_int + hr_int;
  hcr_int2 = hcv_int2 + hr_int;
  Qair = -Qsurf;
  Qair = hcr_int2 * A * (Tint - Tsurf);

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
end InternalSurfaceTodas;

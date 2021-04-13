within DoubleSkinFacades.BasicComponents.AirChannelModels;
model AirCavityThermalBalance
  "This component describes the thermal balance in a generic air cavity of a double skin façade"
  outer Modelica.SIunits.Temperature T_ext;
  outer Modelica.SIunits.Pressure Patm;
  outer Modelica.SIunits.Irradiance Gdiff;
  outer Modelica.SIunits.Irradiance Gdir;
  outer Modelica.SIunits.Angle incidence_angle_in;
  outer Modelica.SIunits.PerUnit tau_dir_ext "Direct transmittance of outer glass";
  parameter Modelica.SIunits.Length s "Thickness of the air cavity [m]";
  parameter Modelica.SIunits.SpecificHeatCapacity cp = 1010 "Specific heat capacity of air [J/(kg.K)]";
  parameter Modelica.SIunits.Height H "Height of the cavity [m]";
  parameter Modelica.SIunits.Height L "Length of the cavity [m]";
  parameter Modelica.SIunits.SpecificEntropy R_air = 8314/29 "Ideal gas constant of air [J/kg.K)]";
  parameter Modelica.SIunits.PerUnit epsilon_glass "Emissivity of internal side of exterior skin";
  parameter Modelica.SIunits.DynamicViscosity mu = 1.8*10^(-5) "Dynamic viscosity of air";
  parameter Real sigma(unit="W/(m2.K4)") = 5.6704*10^(-8) "Stefan-Boltzmann constant [W/(m2.K4)]";
  parameter Modelica.SIunits.Acceleration g = 9.814;
  parameter Modelica.SIunits.PerUnit tau_diff_ext "Diffuse transmittance of outer glass";
  parameter Modelica.SIunits.PerUnit Fbv = 0.4 "View factor basement-façade";
  parameter Modelica.SIunits.PerUnit abs_bas = 1 "Absorbance coefficient of the basement";
  parameter Modelica.SIunits.PerUnit q_conv= 0.5 "Fraction of radiation that is absorbed by the basement given to the air by convection";
  parameter Modelica.SIunits.ThermalConductivity lambda = 0.024;
  final parameter Modelica.SIunits.Area A = H*L "Area of the window façade[m2] ";
  //outer Integer upwards(start=1);

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a glassSurface_a annotation (Placement(transformation(extent={{-18,4},{2,24}}),
        iconTransformation(extent={{-110,-8},{-90,12}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b glassSurface_b annotation (Placement(transformation(extent={{4,-12},{24,8}}),
        iconTransformation(extent={{90,-10},{110,10}})));
  Modelica.Blocks.Interfaces.RealInput v_m "Mean velocity of the air inside the cavity (input)" annotation (Placement(
        transformation(extent={{-66,-112},{-26,-72}}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={60,-100})));
  Modelica.Blocks.Interfaces.RealOutput Tf "Mean temperature of the air inside de cavity" annotation (Placement(
        transformation(extent={{32,28},{52,48}}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={-36,-100})));

  //Variables
  Modelica.SIunits.Temperature T_inlet "Inlet temperature of the cooling flow [K]";
  Modelica.SIunits.Temperature Tsurf_a "Surface temperature of glass a [K]";
  Modelica.SIunits.Temperature Tsurf_b "Surface temperature of glass b [K]";
  Modelica.SIunits.Temperature T_outlet "Temperature of air at the cavity outlet [K]";

  Modelica.SIunits.HeatFlowRate Qa "Power into glassSurface_a [W]";
  Modelica.SIunits.HeatFlowRate Qb "Power to glassSurface_b [W]";
  Modelica.SIunits.HeatFlowRate Qv "Cooling power of the air flow [W]";
  Modelica.SIunits.HeatFlux qa "Heat flux into glassSurface_a [W/m2]";
  Modelica.SIunits.HeatFlux qb "Heat flux to glassSurface_b [W/m2]";
  Modelica.SIunits.HeatFlux qv "Heat flux extracted by ventilation [W/m2]";
  Modelica.SIunits.HeatFlux qconv_a "Heat flux removed by convecton from glass a [W/m2]";
  Modelica.SIunits.HeatFlux qconv_b "Heat flux removed by convecton from glass b [W/m2]";
  Modelica.SIunits.HeatFlux qrad "Heat flux through radiation between the glasses [W/m2]";
  Modelica.SIunits.ThermalInsulance Rv(start = 5) "Cavity cooling power fictitious resistance [m2.K/W]";
  Modelica.SIunits.CoefficientOfHeatTransfer hcv_a "Convective heat transfer coefficient of glass a [W/(m2.K)]";
  Modelica.SIunits.CoefficientOfHeatTransfer hcv_b "Convective heat transfer coefficient of glass b [W/(m2.K)]";
  Modelica.SIunits.CoefficientOfHeatTransfer hr "Radiative heat transfer coefficient between glasses [W/(m2.K)]";
  Modelica.SIunits.Density rho_m(start = 1) "Density of air inside the cavity [kg/m3]";
  Modelica.SIunits.Irradiance G_bas "Radiation absorbed by basement [W/m2]";
  Modelica.SIunits.PerUnit Nu_a;
    Modelica.SIunits.PerUnit Nu_b;
    parameter Modelica.SIunits.PerUnit   Pr = cp*mu/lambda;

  Modelica.SIunits.ThermalDiffusivity kappa "Themal diffusivity";
  Modelica.SIunits.KinematicViscosity nu;

  Modelica.SIunits.PerUnit Nu_for;
  Modelica.SIunits.PerUnit Nu_nat_a;
  Modelica.SIunits.PerUnit Nu_nat_b;
  Modelica.SIunits.PerUnit Nu_nat_a_0;
  Modelica.SIunits.PerUnit Nu_nat_b_0;

  Modelica.SIunits.PerUnit Ra_a;
  Modelica.SIunits.PerUnit Ra_a_0;
    Modelica.SIunits.PerUnit Ra_b;

  Modelica.SIunits.PerUnit Re;
  //Modelica.SIunits.PerUnit Gr;
  Real beta(unit = "K-1");

  Real Ma(unit="kg/(m.s)") "Air flow rate per meter of section [kg/(m.s)]";
  Real K(unit="m-1") "Paramater K for simplification of fictitious resistance [m-1]";
  //Modelica.SIunits.Height H0 "Characteristic height [m]";
  //Modelica.SIunits.Temperature Tmr "Inlet temperature of the cooling flow [K]";
  Modelica.SIunits.TemperatureDifference deltaTf;
  Modelica.SIunits.TemperatureDifference deltaTf2;
  Modelica.SIunits.TemperatureDifference deltaTinf;

equation

  Ma = rho_m * s * v_m;
  kappa = lambda / (rho_m * cp);
  nu = mu / rho_m;
  beta = 1/Tf;

  // Constitutive equations
  Qa = qa * A;
  Qb = qb * A;
  Qv = qv * A;
  Qa = Qb + Qv;

  // Kirchoff's law
  qa = qrad + qconv_a;
  qconv_a = qconv_b + qv;

  qconv_a = hcv_a * (Tsurf_a - Tf);
  qconv_b = hcv_b * (Tf - Tsurf_b);
  qv = 1/Rv * (Tf - T_inlet);

  qrad = hr * (Tsurf_a - Tsurf_b);

  qv = Ma * cp * (T_outlet - T_inlet) / H;

  //MOD
  Ra_a_0 = abs( H^3 * g * beta * deltaTinf / (nu*kappa));
  Ra_b = abs( H^3 * g * beta * deltaTinf / (nu*kappa));

  //Gr = rho_m^2 *g * beta * abs(Tsurf_a-Tf)*H^3/(mu^2);

  Re = rho_m * v_m * H / mu;

 // if abs(Tsurf_a - Tf) < 0.01 then deltaTf = 0.01;
 // else deltaTf = abs(Tsurf_a - Tf);
 //end if;
 //deltaTf2 = Tsurf_a - Tf;

 //if abs(Tsurf_a - T_inlet) < 0.01 then deltaTinf = 0.01;
  //else deltaTinf = abs(Tsurf_a - T_inlet);
  //end if;

  if Ra_a_0 <0.01 then Ra_a = 0.01;
  else Ra_a = Ra_a_0;
  end if;

  Nu_nat_a  = (0.825 + 0.325* (abs(Ra_a))^(1/6))^2;
  Nu_nat_b  = (0.825 + 0.325* (abs(Ra_b))^(1/6))^2;
  Nu_nat_a_0 = Nu_nat_a * abs((Tsurf_a - T_inlet))/(deltaTf);
  Nu_nat_b_0 = Nu_nat_b * abs((Tsurf_b - T_inlet))/(deltaTf);

  //Nu_for = 0.0296 * Pr^(1/3) * abs(Re)^(4/5);
  //Nu_for=0;
  //Nu_for = 0.664 * Pr^(1/3) * Re^(1/2);
  Nu_for = 7.541;
  Nu_a = (abs(Nu_nat_a_0)^3 + abs(Nu_for)^3)^(1/3);
  Nu_b = (abs(Nu_nat_b_0)^3 + abs(Nu_for)^3)^(1/3);

//CORREZIONE rescalare con Tf comprobar el deltaT del rayleigh!!
  hcv_a = Nu_a * lambda / H;
  hcv_b = Nu_b * lambda / H;

  // Radiation
  hr = 4 * sigma * (1 / (1/epsilon_glass + 1/epsilon_glass - 1)) * ((Tsurf_a + Tsurf_b)/2)^3;

  // Fictitious resistance
  rho_m = Patm/(R_air*Tf);
  K = (hcv_a + hcv_b)/(s*v_m*rho_m*cp);
  Rv = (K*H - 1 + exp(-K*H))/((hcv_a+hcv_b)*(1-exp(-K*H)));

  T_inlet = T_ext + (G_bas * q_conv * s/(Ma*cp));

  G_bas = abs_bas * (tau_dir_ext * Gdir * tan(incidence_angle_in) + tau_diff_ext * Gdiff * Fbv);

  // Boundary variables
  Qa = glassSurface_a.Q_flow;
  Qb = -glassSurface_b.Q_flow;
  Tsurf_a = glassSurface_a.T;
  Tsurf_b = glassSurface_b.T;

   annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200}),
          Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}),
      Diagram(coordinateSystem(preserveAspectRatio=false)));
end AirCavityThermalBalance;

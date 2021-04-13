within DoubleSkinFacades.BasicComponents.AirChannelModels;
model AirCavityThermalBalanceChuChu
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
  Modelica.Blocks.Interfaces.RealOutput Tf(final unit = "K", displayUnit= "deg") "Mean temperature of the air inside de cavity" annotation (Placement(
        transformation(extent={{32,28},{52,48}}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={-36,-100})));

  //Variables
  Modelica.SIunits.Temperature T_inlet "Inlet temperature of the cooling flow [K]";
  Modelica.SIunits.Temperature Tsurf_a "Surface temperature of glass a [K]";
  Modelica.SIunits.Temperature Tsurf_b "Surface temperature of glass b [K]";
  Modelica.SIunits.Temperature T_outlet "Temperature of air at the cavity outlet [K]";
  Modelica.SIunits.Temperature T_inf "Temperature of air at the cavity outlet [K]";

  Modelica.SIunits.HeatFlowRate Qa "Power into glassSurface_a [W]";
  Modelica.SIunits.HeatFlowRate Qb "Power to glassSurface_b [W]";
  Modelica.SIunits.HeatFlowRate Qv "Cooling power of the air flow [W]";
  Modelica.SIunits.HeatFlux qa "Heat flux into glassSurface_a [W/m2]";
  Modelica.SIunits.HeatFlux qb "Heat flux to glassSurface_b [W/m2]";
  Modelica.SIunits.HeatFlux qv "Heat flux extracted by ventilation [W/m2]";
  Modelica.SIunits.HeatFlux qconv_a "Heat flux removed by convecton from glass a [W/m2]";
  Modelica.SIunits.HeatFlux qconv_b "Heat flux removed by convecton from glass b [W/m2]";
  Modelica.SIunits.HeatFlux qrad "Heat flux through radiation between the glasses [W/m2]";
  //Modelica.SIunits.ThermalInsulance Rv(start = 5) "Cavity cooling power fictitious resistance [m2.K/W]";
  Modelica.SIunits.CoefficientOfHeatTransfer hcv_1 "Convective heat transfer coefficient of glass a related to inlet [W/(m2.K)]";
  Modelica.SIunits.CoefficientOfHeatTransfer hcv_2 "Convective heat transfer coefficient of glass b related to inlet [W/(m2.K)]";
  Modelica.SIunits.CoefficientOfHeatTransfer hcv_a "Convective heat transfer coefficient of glass a [W/(m2.K)]";
  Modelica.SIunits.CoefficientOfHeatTransfer hcv_b "Convective heat transfer coefficient of glass b [W/(m2.K)]";
  Modelica.SIunits.CoefficientOfHeatTransfer hr "Radiative heat transfer coefficient between glasses [W/(m2.K)]";
  Modelica.SIunits.Density rho_m(start = 1) "Density of air inside the cavity [kg/m3]";
  Modelica.SIunits.Irradiance G_bas "Radiation absorbed by basement [W/m2]";
  parameter Modelica.SIunits.PerUnit   Pr = cp*mu/lambda;

  Modelica.SIunits.ThermalDiffusivity kappa "Themal diffusivity";
  Modelica.SIunits.KinematicViscosity nu;
  Modelica.SIunits.PerUnit Re;
  Modelica.SIunits.PerUnit Ra_a;
  Modelica.SIunits.PerUnit Ra_b;
  Modelica.SIunits.PerUnit Nu_a;
  Modelica.SIunits.PerUnit Nu_b;
  Modelica.SIunits.PerUnit Nu_for;
  Modelica.SIunits.PerUnit Nu_nat_a;
  Modelica.SIunits.PerUnit Nu_nat_b;

  Real beta(unit = "K-1");
  Real Ma(unit="kg/(m.s)") "Air flow rate per meter of section [kg/(m.s)]";
  Real c(unit="m-1") "Paramater c for simplification of fictitious resistance [m-1]";
  Real K(unit="m-1") "Paramater K for simplification of fictitious resistance [m-1]";
  Modelica.SIunits.ThermalInsulance Rv(start = 0.2) "Cavity cooling power fictitious resistance [m2.K/W]";

initial equation

equation

  Ma = rho_m * s * v_m;
  kappa = lambda / (rho_m * cp);
  nu = mu / rho_m;
  beta = 1/Tf;
  rho_m = Patm/(R_air*Tf);
  T_inf= (hcv_1*Tsurf_a + hcv_2*Tsurf_b)/(hcv_1+hcv_2);

  // Constitutive equations
  Qa = qa * A;
  Qb = qb * A;
  Qv = qv * A;

  // Kirchoff's law
  qa = -qrad + qconv_a;
  qb = qrad + qconv_b; // added
  qconv_a + qconv_b + qv = 0;
  qconv_a = hcv_1*(T_inlet - Tsurf_a);
  qconv_b = hcv_2*(T_inlet - Tsurf_b);
  qrad = hr * (Tsurf_a - Tsurf_b);
  qv = Ma * cp * (T_outlet - T_inlet) / H;

  // Fictitious resistance
  qv = 1/Rv * (Tf - T_inlet);
  c = (T_inlet-T_inf)/(Tf-T_inf);
  K = c*(hcv_1 + hcv_2)/(s*v_m*rho_m*cp);
  Rv = (K*H - 1 + exp(-K*H))/(c*(hcv_1+hcv_2)*(1-exp(-K*H)));

  //Correlazione di Churchill e Chu per lastra piana
  Ra_a = abs( H^3 * g * beta * (Tsurf_a-T_inlet) / (nu*kappa));
  Ra_b = abs( H^3 * g * beta * (Tsurf_b-T_inlet) / (nu*kappa));
  Re = rho_m * v_m * H / mu;

   // Calcolo in naturale per vetro esterno
   Nu_nat_a  = (0.825 + 0.325* (Ra_a)^(1/6))^2;

   //Nu_nat_a  = DoubleSkinFacades.AuxiliaryBlocks.regStep(
     //        y1=(0.825 + 0.325*(Ra_a)^(1/6))^2,
       //      y2=0.68+0.515*(Ra_a)^(1/4),
         //    x=Ra_a-10^9,
           //  x_small=10^8);

   // Calcolo in naturale per vetro interno
   Nu_nat_b  = (0.825 + 0.325* (Ra_b)^(1/6))^2;

   //Nu_nat_b  = DoubleSkinFacades.AuxiliaryBlocks.regStep(
     //        y1=(0.825 + 0.325*(Ra_b)^(1/6))^2,
       //      y2=0.68+0.515*(Ra_b)^(1/4),
         //    x=Ra_b-10^9,
           //  x_small=10^8);

  // Calcolo in forzata
  Nu_for = 0.664 * Pr^(1/3) * Re^(1/2);

  // Convezione mista
  Nu_a = (abs(Nu_nat_a)^3 + abs(Nu_for)^3)^(1/3);
  Nu_b = (abs(Nu_nat_b)^3 + abs(Nu_for)^3)^(1/3);

  hcv_1 = Nu_a * lambda / H;
  hcv_2 = Nu_b * lambda / H;
  hcv_a = c*hcv_1;
  hcv_b = c*hcv_2;
  // Radiation
  hr = 4 * sigma * (1 / (1/epsilon_glass + 1/epsilon_glass - 1)) * ((Tsurf_a + Tsurf_b)/2)^3;

  // preheating
  T_inlet = T_ext + (G_bas * q_conv * s/(Ma*cp));
  G_bas = abs_bas * (tau_dir_ext * Gdir * tan(incidence_angle_in) + tau_diff_ext * Gdiff * Fbv);

  // Boundary variables
  Qa = -glassSurface_a.Q_flow;
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
end AirCavityThermalBalanceChuChu;

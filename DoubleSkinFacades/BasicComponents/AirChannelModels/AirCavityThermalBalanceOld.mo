within DoubleSkinFacades.BasicComponents.AirChannelModels;
model AirCavityThermalBalanceOld
  "This component describes the thermal balance in a generic air cavity of a double skin façade"
  outer Modelica.SIunits.Temperature T_ext;
  outer Modelica.SIunits.Pressure Patm;
  outer Modelica.SIunits.Irradiance Gdiff;
  outer Modelica.SIunits.Irradiance Gdir;
  outer Modelica.SIunits.Angle incidence_angle;
  outer Modelica.SIunits.PerUnit tau_dir_ext "Direct transmittance of outer glass";
  parameter Modelica.SIunits.Length s "Thickness of the air cavity [m]";
  parameter Modelica.SIunits.Height H "Height of the cavity [m]";
 // Modelica.SIunits.Height H0 "Characteristic height [m]";
  parameter Modelica.SIunits.Height L "Length of the cavity [m]";
  parameter Modelica.SIunits.Acceleration g = 9.814 "Acceleration gravity";
  parameter Media.Gases.BaseGas material_gas=
      DoubleSkinFacades.Media.Gases.Air()
    "Type of internal gas between the glasses"
    annotation (choicesAllMatching=true);
  parameter Real sigma(unit="W/(m2.K4)") = 5.6704*10^(-8) "Stefan-Boltzmann constant";
  final parameter Modelica.SIunits.Area A = H * L "Area of the window";
  parameter Modelica.SIunits.PerUnit cp_a = material_gas.cp_a "Specific heat capacity coefficient";
  parameter Modelica.SIunits.PerUnit cp_b = material_gas.cp_b "Specific heat capacity coefficient";
  parameter Modelica.SIunits.PerUnit mu_a = material_gas.mu_a "Dynamic viscosity coefficient";
  parameter Modelica.SIunits.PerUnit mu_b = material_gas.mu_b "Dynamic viscosity coefficient";
  parameter Modelica.SIunits.PerUnit lambda_a = material_gas.lambda_a "Thermal conductivity coefficient";
  parameter Modelica.SIunits.PerUnit lambda_b = material_gas.lambda_b "Thermal conductivity coefficient";
  parameter Modelica.SIunits.MolarMass Mm = material_gas.Mm "Molecular mass";
  parameter Modelica.SIunits.SpecificEntropy R_air = 8314/Mm "Ideal gas constant of air [J/kg.K)]";
  parameter Modelica.SIunits.PerUnit epsilon_glass "Emissivity of internal side of exterior skin";
  parameter Modelica.SIunits.PerUnit tau_diff_ext "Diffuse transmittance of outer glass";
  parameter Modelica.SIunits.PerUnit Fbv = 0.4 "View factor basement-façade";
  parameter Modelica.SIunits.PerUnit abs_bas = 1 "Absorbance coefficient of the basement";
  parameter Modelica.SIunits.PerUnit q_conv= 0.5 "Fraction of radiation that is absorbed by the basement given to the air by convection";
  outer Integer upwards;

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
  Modelica.SIunits.Temperature Tinlet "Inlet temperature of the cooling flow [K]";
  Modelica.SIunits.Temperature Tsurf_a "Surface temperature of glass a [K]";
  Modelica.SIunits.Temperature Tsurf_b "Surface temperature of glass b [K]";
  Modelica.SIunits.Temperature Toutlet "Temperature of air at the cavity outlet [K]";
  Modelica.SIunits.Temperature Tm "Mean temperatue of glasses [K]";
  Modelica.SIunits.HeatFlowRate Qa "Power into glassSurface_a [W]";
  Modelica.SIunits.HeatFlowRate Qb "Power to glassSurface_b [W]";
  Modelica.SIunits.HeatFlowRate Qv "Cooling power of the air flow [W]";
  Modelica.SIunits.HeatFlux qa "Heat flux into glassSurface_a [W/m2]";
  Modelica.SIunits.HeatFlux qb "Heat flux to glassSurface_b [W/m2]";
  Modelica.SIunits.HeatFlux qv "Heat flux extracted by ventilation [W/m2]";
  Modelica.SIunits.HeatFlux qconv_a "Heat flux removed by convecton from glass a [W/m2]";
  Modelica.SIunits.HeatFlux qconv_b "Heat flux removed by convecton from glass b [W/m2]";
  Modelica.SIunits.HeatFlux qrad "Heat flux through radiation between the glasses [W/m2]";
  Modelica.SIunits.CoefficientOfHeatTransfer hcv "Convective heat transfer coefficient of closed cavity [W/(m2.K)]";
  Modelica.SIunits.CoefficientOfHeatTransfer hcv_a "Convective heat transfer coefficient of glass a [W/(m2.K)]";
  Modelica.SIunits.CoefficientOfHeatTransfer hcv_b "Convective heat transfer coefficient of glass b [W/(m2.K)]";
  Modelica.SIunits.CoefficientOfHeatTransfer hr "Radiative heat transfer coefficient between glasses [W/(m2.K)]";
  Modelica.SIunits.PerUnit Ra(start=10000) "Grashof number";
  Modelica.SIunits.PerUnit Nu "Nusselt number";
  Modelica.SIunits.PerUnit Nu1 "Nusselt number";
  Modelica.SIunits.PerUnit Nu2 "Nusselt number";
  Modelica.SIunits.Density rho_m(start = 1) "Density of air inside the cavity [kg/m3]";
  Modelica.SIunits.Irradiance G_bas "Radiation absorbed by basement [W/m2]";
  Real Ma(unit="kg/(m.s)") "Air flow rate per meter of section [kg/(m.s)]";
  Modelica.SIunits.SpecificHeatCapacity cp;
  Modelica.SIunits.PerUnit lambda;
  Modelica.SIunits.DynamicViscosity mu "Dynamic viscosity of internal gas";
  Modelica.SIunits.TemperatureDifference deltaT(start=0.01);
  Real beta(final unit "1/K");
  Real K;
  Real Rv;

  parameter Modelica.SIunits.PerUnit A_cavity = H/s;

equation

  Ma = rho_m * s * v_m;

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
  qrad = hr * (Tsurf_a - Tsurf_b);

  qv = 1/Rv * (Tf - Tinlet);
  K = (hcv_a + hcv_b)/(s*v_m*rho_m*cp);
  Rv = (K*H - 1 + exp(-K*H))/((hcv_a+hcv_b)*(1-exp(-K*H)));

  qv = Ma * cp * (Toutlet - Tinlet) / H;

 // Convective coefficients
  Tm = (Tsurf_a + Tsurf_b)/2;

  //Gas properties
  rho_m = Patm/(R_air*Tf);
  cp = cp_a + cp_b*Tm;
  mu = mu_a + mu_b*Tm;
  lambda = lambda_a + lambda_b*Tm;

  //Convective coefficient

    if (abs(Tsurf_a - Tsurf_b) <= 0.01) then
        deltaT = 0.01;
    else
        deltaT = abs(Tsurf_a - Tsurf_b);
    end if;

  beta = 1 / Tm;
  Ra = abs( rho_m^2 * s^3 * g * beta * cp * deltaT / (mu*lambda));
  if Ra>= 5*10^4 then Nu1 = 0.067383 * abs(Ra)^(1/3);
  elseif Ra< 5*10^4 and Ra> 10^4 then Nu1 = 0.028154 * abs(Ra)^0.413;
  else Nu1 = 1 + 1.7596678*10^(-10) * abs(Ra)^2.2984755;
  end if;
  Nu2 = 0.242 * abs(Ra/A_cavity)^0.272;

  Nu = max(Nu1,Nu2);
  hcv =  Nu * lambda / H;

  hcv_a = hcv + 4*v_m;
  hcv_b = hcv + 4*v_m;

  // Radiation
  hr = 4 * sigma * (1 / (1/epsilon_glass + 1/epsilon_glass - 1)) * (Tm)^3;

  //Temperature profile
  //H0 = rho_m * cp * s * v_m / (hcv_a + hcv_b);
  //Toutlet = Tm - (Tm - Tinlet)* exp(-H0/H);
  //Tf = Tm - H0/H * (Toutlet - Tinlet);

  //Inlet temperature calculation
  G_bas = abs_bas * (tau_dir_ext * Gdir * tan(incidence_angle) + tau_diff_ext * Gdiff * Fbv);
  Tinlet = T_ext + upwards*(G_bas * q_conv * s/(Ma*cp));

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
end AirCavityThermalBalanceOld;

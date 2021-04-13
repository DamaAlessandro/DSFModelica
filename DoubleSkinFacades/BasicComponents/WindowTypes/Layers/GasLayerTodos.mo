within DoubleSkinFacades.BasicComponents.WindowTypes.Layers;
model GasLayerTodos "Double glass window with an intermediate gas gap"
  parameter Modelica.SIunits.Height H "Height of the window";
  parameter Modelica.SIunits.Height L "Length of the window";
  parameter Modelica.SIunits.Length s_gap "Thickness of the air gap";
  parameter Modelica.SIunits.Acceleration g = 9.814 "Acceleration gravity";
  parameter Media.Materials.BaseGlass material_glass "Type of glass" annotation(choicesAllMatching = true);
  parameter Modelica.SIunits.PerUnit epsilon_glass "Emissivity of the glass surface";
  parameter Modelica.SIunits.PerUnit epsilon_glazing "Emissivity of the glazing";
  parameter Media.Gases.BaseGas material_gas
    "Type of internal gas between the glasses"
    annotation (choicesAllMatching=true);
  parameter Modelica.SIunits.SpecificHeatCapacity cp_gas = material_gas.cp "Specific heat capacity of internal gas";
  parameter Modelica.SIunits.ThermalConductivity lambda_gas = material_gas.lambda "Thermal conductivity of internal gas";
  parameter Modelica.SIunits.SpecificEntropy R_gas = material_gas.R "Gas constant of the gas inside the gap";
  parameter Real sigma(unit="W/(m2.K4)") = 5.6704*10^(-8) "Stefan-Boltzmann constant";
  parameter Modelica.SIunits.Temperature Tstart = 273.15 + 25 "Initial temperature, all layers";
  parameter Modelica.SIunits.Pressure P_gap = 101325 "Pressure inside the gas gap";
  final parameter Modelica.SIunits.Area A = H * L "Area of the window";
  Modelica.SIunits.PerUnit Pr "Prandtl number";
  Modelica.SIunits.DynamicViscosity mu "Dynamic viscosity of internal gas";
  Modelica.SIunits.TemperatureDifference deltaT(start=0.01);
  Modelica.SIunits.Mass M_gas "Mass of the gas inside the gap";
  Modelica.SIunits.Density rho_gas(start = 1.225) "Density of the gas  inside the gap";
  Modelica.SIunits.HeatFlowRate Qa "Power into glassSurface_a [W]";
  Modelica.SIunits.HeatFlowRate Qb "Power to glassSurface_b [W]";
  Modelica.SIunits.HeatFlux qa "Heat flux into glassSurface_a [W/m2]";
  Modelica.SIunits.HeatFlux qb "Heat flux to glassSurface_b [W/m2]";
  Modelica.SIunits.Temperature Tmr(start = Tstart) "Mean radiant temperature inside the gas gap";
  Modelica.SIunits.Temperature Tsurf_a(start = 295) "Inner surface temperature of glass a";
  Modelica.SIunits.Temperature Tsurf_b(start = 290) "Inner surface temperature of glass b";
  Modelica.SIunits.PerUnit Gr "Grashof number";
  Modelica.SIunits.PerUnit Nu "Nusselt number";
  Modelica.SIunits.CoefficientOfHeatTransfer hcv_gap "Convective heat transfer coefficient in the gas gap";
  Modelica.SIunits.CoefficientOfHeatTransfer hr_gap "Radiative heat transfer coefficient in the gas gap";
  Modelica.SIunits.CoefficientOfHeatTransfer hcr_gap "Convective-radiative heat transfer coefficient in the gas gap";
  Modelica.SIunits.ThermalResistance Rcr_gap;

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a glassSurface_a annotation (
      Placement(transformation(extent={{-110,-10},{-90,10}}),
        iconTransformation(extent={{-110,-10},{-90,10}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b glassSurface_b annotation (
      Placement(transformation(extent={{90,-10},{110,10}}), iconTransformation(
          extent={{90,-10},{110,10}})));

equation

  // Pressure inside de gap assumed the same as outside pressure

    if (abs(Tsurf_a - Tsurf_b) <= 0.01) then
        deltaT = 0.01;
    else
        deltaT = abs(Tsurf_a - Tsurf_b);
    end if;

  Tmr = (Tsurf_a + Tsurf_b)/2;

  //Gas properties
  rho_gas = P_gap / (R_gas * Tmr);
  M_gas = rho_gas * L * H * s_gap;

  //Convective coefficient
  mu = (1.458*(10^(-6))*(Tmr^1.5))/(Tmr + 110.4);
  Gr = rho_gas^2 * s_gap^3 * g * deltaT / Tmr / mu^2;
  Pr = cp_gas * mu / lambda_gas;
  Nu = 0.035 * (Gr * Pr)^0.38;
  hcv_gap = 1.4*max(Nu,1) * lambda_gas / s_gap;

  //Radiative coefficient
  hr_gap = 4 * sigma * (1 / (1/epsilon_glass + 1/epsilon_glazing - 1)) * (Tmr)^3;

  //Radiative convective coefficient
  hcr_gap = hcv_gap + hr_gap;
  Rcr_gap = 1/(hcr_gap*A);

  //Thermal balance

  Qa = Qb;
  Qa = 1/Rcr_gap * (Tsurf_a - Tsurf_b);
  qa = Qa / A;
  qb = Qb / A;

  //Constitutive equations
  Qa = glassSurface_a.Q_flow;
  Qb = - glassSurface_b.Q_flow;
  Tsurf_a = glassSurface_a.T;
  Tsurf_b = glassSurface_b.T;

  annotation (
  Icon(graphics={  Rectangle(extent = {{-30, 100}, {30, -100}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
            fillPattern =  FillPattern.Solid), Rectangle(extent = {{-40, 100}, {-20, -100}}, lineColor = {0, 0, 0}, fillColor = {213, 245, 255},
            fillPattern =  FillPattern.Solid), Rectangle(extent = {{20, 100}, {40, -100}}, lineColor = {0, 0, 0}, fillColor = {213, 245, 255},
            fillPattern =  FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}})));
end GasLayerTodos;

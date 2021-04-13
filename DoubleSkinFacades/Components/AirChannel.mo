within DoubleSkinFacades.Components;
model AirChannel
  parameter Modelica.SIunits.Length s "Thickness of the air cavity";
  parameter Modelica.SIunits.Height H "Height of the channel";
  parameter Modelica.SIunits.Height L "Length of the channel";
  parameter Modelica.SIunits.Area A_in "Area of air cavity inlet";
  parameter Modelica.SIunits.Area A_out "Area of air cavity outlet";
  parameter Modelica.SIunits.DynamicViscosity mu = 2.27*10^(-5) "Dynamic viscosity of air";
  parameter Modelica.SIunits.PerUnit epsilon_glass "Emissivity of glass";
  parameter Modelica.SIunits.PerUnit tau_diff_ext;
  parameter Real sigma(unit="W/(m2.K4)") = 5.6704*10^(-8) "Stefan-Boltzmann constant";
  parameter Real orientation(final unit="deg") "Façade orientation angle in degrees" annotation (Dialog(group="Façade specifications"));
  inner Integer upwards;

  BasicComponents.AirChannelModels.AirCavityThermalBalancePre
    AirCavityThermalBalance(
    H=H,
    s=s,
    L=L,
    epsilon_glass=epsilon_glass,
    tau_diff_ext=tau_diff_ext)
    annotation (Placement(transformation(extent={{-10,12},{10,32}})));

  replaceable
    BasicComponents.AirChannelModels.VentilationType.NaturalVentilation
    ventilationType(
    A_in=A_in,
    A_out=A_out,
    s=s,
    H=H,
    cpValuesFile=
        "C:/Users/damaa/Documents/Ricerca/DSF/DSF_modelica/DSF_Library/CpValues.txt",
    orientation=orientation) annotation (choices(
      choice(redeclare
          DoubleSkinFacades.BasicComponents.AirChannelModels.VentilationType.NaturalVentilation
          ventilationType),
      choice(redeclare
          DoubleSkinFacades.BasicComponents.AirChannelModels.VentilationType.ForzedVentilation
          ventilationType),
      Dialog(group="Ventilation characteristics")), Placement(
        transformation(extent={{-10,-34},{10,-14}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));

equation
  upwards=ventilationType.upwards;

  connect(ventilationType.Tf, AirCavityThermalBalance.Tf) annotation (Line(
        points={{-4.4,-13.8},{-4.4,2.1},{-3.6,2.1},{-3.6,12}},color={0,0,127}));
  connect(AirCavityThermalBalance.glassSurface_a, port_a) annotation (Line(
        points={{-10,22.2},{-56,22.2},{-56,0},{-100,0}}, color={191,0,0}));
  connect(AirCavityThermalBalance.glassSurface_b, port_b) annotation (Line(
        points={{10,22},{54,22},{54,0},{100,0}}, color={191,0,0}));
  connect(ventilationType.v_m, AirCavityThermalBalance.v_m) annotation (
      Line(points={{4.5,-13.9},{4.5,2.05},{6,2.05},{6,12}},color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-100,-100},{100,100}},
          lineColor={28,108,200},
          fillColor={170,255,255},
          fillPattern=FillPattern.Solid)}),                      Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=1294800,
      Interval=600,
      Tolerance=1e-05,
      __Dymola_Algorithm="Esdirk45a"));
end AirChannel;

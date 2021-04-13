within DoubleSkinFacades.Components;
block ExternalWindow
  parameter Modelica.SIunits.Height H "Height of the façade";
  parameter Modelica.SIunits.Height L "Length of the façade";
  parameter Modelica.SIunits.Length s_glass "Thickness of glass layer";
  parameter Modelica.SIunits.Length s_gap "Thickness of gas gap (if doubleglass)";
  parameter Media.Materials.BaseGlass material_glass "Type of glass"   annotation(choicesAllMatching = true);
  parameter Media.Gases.BaseGas material_gas "Type of gas in the gap"
    annotation (choicesAllMatching=true);
  parameter Modelica.SIunits.Temperature T_ground = 273.15 + 16 "Average temperature of the ground for sky temperature calculation";
  parameter Modelica.SIunits.PerUnit epsilon_glass;
  parameter Modelica.SIunits.PerUnit epsilon_glazing;
  parameter Modelica.SIunits.PerUnit alpha1_diff "Diffuse absorptivity of external glass";
  parameter Modelica.SIunits.PerUnit alpha2_diff "Diffuse absorptivity of internal glass";
  parameter Real orientation(final unit="deg") "Façade orientation angle in degrees" annotation (Dialog(group="Façade specifications"));
  outer Modelica.SIunits.PerUnit alpha1_dir_ext "Direct absorptivity of external glass";
  outer Modelica.SIunits.PerUnit alpha2_dir_ext "Direct absorptivity of internal glass";
  inner Modelica.SIunits.PerUnit alpha1_dir "Direct absorptivity of external glass";
  inner Modelica.SIunits.PerUnit alpha2_dir "Direct absorptivity of internal glass";

  replaceable BasicComponents.WindowTypes.DoubleGlass glassType(
    s_gap=s_gap,
    s_glass=s_glass,
    H=H,
    L=L,
    Tstart=290.15,
    alpha1_diff=alpha1_diff,
    alpha2_diff=alpha2_diff,
    material_gas=material_gas,
    material_glass=material_glass,
    epsilon_glass=epsilon_glass,
    epsilon_glazing=epsilon_glazing) annotation (choices(choice(redeclare
          DoubleSkinFacades.BasicComponents.WindowTypes.SingleGlass
          glassType), choice(redeclare
          DoubleSkinFacades.BasicComponents.WindowTypes.DoubleGlass
          glassType)), Placement(transformation(extent={{12,-18},{48,18}})));

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b annotation (
      Placement(transformation(extent={{90,-10},{110,10}}),iconTransformation(
          extent={{90,-10},{110,10}})));

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

  BoundaryHeatTransfer.ExternalSurface outsideSurface(
    H=H,
    L=L,
    T_ground=T_ground,
    orientation=orientation)
    annotation (Placement(transformation(extent={{-64,-22},{-22,22}})));

equation

  alpha1_dir = alpha1_dir_ext;
  alpha2_dir = alpha2_dir_ext;
  glassType.shade_factor = 0;

  connect(glassType.port_b, port_b)
    annotation (Line(points={{40.8,0},{100,0}}, color={191,0,0}));
  connect(outsideSurface.glassSurface, glassType.port_a)
    annotation (Line(points={{-22,0},{19.2,0}}, color={191,0,0}));
  connect(outsideSurface.airExterior, port_a)
    annotation (Line(points={{-64,0},{-80,0},{-80,0},{-100,0}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{
            -120,-100},{120,100}}),                             graphics={
          Rectangle(
          extent={{-30,100},{30,-100}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-88,86},{-64,62}},
          lineColor={255,255,85},
          fillColor={255,255,85},
          fillPattern=FillPattern.Solid),
        Line(points={{-74,56},{-74,46}}, color={255,255,85}),
        Line(points={{-58,74},{-48,74}}, color={255,255,85}),
        Line(points={{-76,90},{-76,98}}, color={255,255,85}),
        Line(points={{-92,72},{-98,72}}, color={255,255,85})}),  Diagram(
        coordinateSystem(preserveAspectRatio=false, extent={{-120,-100},{
            120,100}}), graphics={Text(
          extent={{-52,30},{-18,24}},
          lineColor={0,0,0},
          textString="ExteriorSurface"), Text(
          extent={{14,32},{48,26}},
          lineColor={0,0,0},
          textString="DoubleGlass")}));
end ExternalWindow;

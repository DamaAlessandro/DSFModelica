within DoubleSkinFacades.BasicComponents.WindowTypes;
model DoubleGlass "Double glass window with an intermediate gas gap"
  extends DoubleSkinFacades.BasicComponents.WindowTypes.PartialGlass;
  parameter Modelica.SIunits.Height H;
  parameter Modelica.SIunits.Height L;
  parameter Modelica.SIunits.Length s_gap "Thickness of the air gap";
  parameter Modelica.SIunits.Length s_glass "Thickness of the glass layer";
  parameter Media.Materials.BaseGlass material_glass "Type of glass";
  parameter Media.Gases.BaseGas material_gas
    "Type of internal gas between the glasses";
  parameter Modelica.SIunits.Temperature Tstart = 273.15 + 15 "Start temperature of the glass";
  parameter Modelica.SIunits.PerUnit epsilon_glass;
  parameter Modelica.SIunits.PerUnit epsilon_glazing;
  parameter Modelica.SIunits.PerUnit alpha1_diff "Absorptivity of glass";
  parameter Modelica.SIunits.PerUnit alpha2_diff "Absorptivity of glass";
  outer Modelica.SIunits.PerUnit alpha1_dir;
  outer Modelica.SIunits.PerUnit alpha2_dir;
  Modelica.SIunits.PerUnit shade_factor;

  Layers.GasLayer gasGap(
    s_gap=s_gap,
    H=H,
    L=L,
    material_glass=material_glass,
    material_gas=material_gas,
    epsilon_glass=epsilon_glass,
    epsilon_glazing=epsilon_glazing)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  Layers.GlassLayer glassLayer1(
    s_glass=s_glass,
    H=H,
    L=L,
    material=material_glass,
    alpha_diff=alpha1_diff)
    annotation (Placement(transformation(extent={{-42,-10},{-22,10}})));
  Layers.GlassLayer glassLayer2(
    s_glass=s_glass,
    H=H,
    L=L,
    material=material_glass,
    alpha_diff=alpha2_diff)
    annotation (Placement(transformation(extent={{22,-10},{42,10}})));
equation

  glassLayer1.alpha_dir = alpha1_dir;
  glassLayer2.alpha_dir = alpha2_dir;
  glassLayer1.shade_factor = shade_factor;
  glassLayer2.shade_factor = shade_factor;

  connect(glassLayer1.glassSurface_a, port_a) annotation (Line(points={{-36.2,0},
          {-60,0}},                 color={191,0,0}));
  connect(glassLayer1.glassSurface_b, gasGap.glassSurface_a)
    annotation (Line(points={{-27.8,0},{-10,0}}, color={191,0,0}));
  connect(gasGap.glassSurface_b, glassLayer2.glassSurface_a)
    annotation (Line(points={{10,0},{27.8,0}}, color={191,0,0}));
  connect(port_b, glassLayer2.glassSurface_b)
    annotation (Line(points={{60,0},{36.2,0}}, color={191,0,0}));
  annotation (Icon(graphics={
        Rectangle(
          extent={{-40,80},{-10,-80}},
          lineColor={28,108,200},
          fillColor={216,255,251},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{10,80},{40,-80}},
          lineColor={28,108,200},
          fillColor={216,255,251},
          fillPattern=FillPattern.Solid),
        Rectangle(extent={{-60,100},{60,-100}}, lineColor={28,108,200})}));
end DoubleGlass;

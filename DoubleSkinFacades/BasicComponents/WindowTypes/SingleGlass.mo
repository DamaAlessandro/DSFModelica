within DoubleSkinFacades.BasicComponents.WindowTypes;
block SingleGlass
  extends DoubleSkinFacades.BasicComponents.WindowTypes.PartialGlass;
  parameter Modelica.SIunits.Height H "EY Height of the window";
  parameter Modelica.SIunits.Height L "Length of the window";
  parameter Modelica.SIunits.Length s_glass "Thickness of the glass layer";
  parameter Modelica.SIunits.Length s_gap "Thickness of the gap layer";
  parameter Media.Materials.BaseGlass material_glass "Type of glass" annotation(choicesAllMatching = true);
  parameter Media.Gases.BaseGas material_gas "Type of glass";
  parameter Modelica.SIunits.Temperature Tstart = 273.15 + 15 "Start temperature of the glass";
  parameter Modelica.SIunits.PerUnit epsilon_glass;
  parameter Modelica.SIunits.PerUnit epsilon_glazing;
  parameter Modelica.SIunits.PerUnit alpha1_diff "Absorptivity of glass";
  parameter Modelica.SIunits.PerUnit alpha2_diff "Absorptivity of glass";
  outer Modelica.SIunits.PerUnit alpha1_dir;
  Modelica.SIunits.PerUnit shade_factor;

  Layers.GlassLayer glassLayer1(
    s_glass=s_glass,
    H=H,
    L=L,
    material=material_glass,
    Tstart=Tstart,
    alpha_diff=alpha1_diff)
    annotation (Placement(transformation(extent={{-18,-18},{18,18}})));

equation

  glassLayer1.alpha_dir = alpha1_dir;
  glassLayer1.shade_factor = shade_factor;

  connect(port_a, glassLayer1.glassSurface_a) annotation (Line(points={{-60,0},
          {-7.56,0}},                    color={191,0,0}));

  connect(glassLayer1.glassSurface_b, port_b) annotation (Line(points={{7.56,0},
          {60,0}},                   color={191,0,0}));
  connect(port_a, port_a)
    annotation (Line(points={{-60,0},{-60,0}}, color={191,0,0}));
  annotation (Icon(graphics={Rectangle(extent={{-60,100},{60,-100}},
            lineColor={28,108,200}),
        Rectangle(
          extent={{-16,80},{16,-80}},
          lineColor={28,108,200},
          fillColor={216,255,251},
          fillPattern=FillPattern.Solid)}));
end SingleGlass;

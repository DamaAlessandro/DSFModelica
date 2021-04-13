within DoubleSkinFacades.Components;
model InternalWindow
  parameter Modelica.SIunits.Height H "Height of the façade";
  parameter Modelica.SIunits.Height L "Length of the façade";
  parameter Modelica.SIunits.Length s "Width of air cavity";
  parameter Modelica.SIunits.Length s_glass "Thickness of glass layer";
  parameter Modelica.SIunits.Length s_gap "Thickness of gas gap (if doubleglass)";
  parameter Media.Materials.BaseGlass material_glass "Type of glass"   annotation(choicesAllMatching = true);
  parameter Media.Gases.BaseGas material_gas
    "Type of internal gas between the glasses"
    annotation (choicesAllMatching=true);
  parameter Modelica.SIunits.PerUnit epsilon_glass;
  parameter Modelica.SIunits.PerUnit epsilon_glazing;
  parameter Modelica.SIunits.PerUnit alpha1_diff "Diffuse absorptivity of external glass";
  parameter Modelica.SIunits.PerUnit alpha2_diff "Diffuse absorptivity of internal glass";
  outer Modelica.SIunits.PerUnit alpha1_dir_int "Direct absorptivity of external glass";
  outer Modelica.SIunits.PerUnit alpha2_dir_int "Direct absorptivity of internal glass";
  inner Modelica.SIunits.PerUnit alpha1_dir "Direct absorptivity of external glass";
  inner Modelica.SIunits.PerUnit alpha2_dir "Direct absorptivity of internal glass";
  outer Modelica.SIunits.Angle incidence_angle;
  outer Modelica.SIunits.Angle elevation;
  outer Modelica.SIunits.Angle HSA;
  Real s1;
  Real s2;

  replaceable BasicComponents.WindowTypes.DoubleGlass glassType(
    s_gap=s_gap,
    s_glass=s_glass,
    H=H,
    L=L,
    alpha1_diff=alpha1_diff,
    alpha2_diff=alpha2_diff,
    material_glass=material_glass,
    material_gas=material_gas,
    epsilon_glass=epsilon_glass,
    epsilon_glazing=epsilon_glazing) annotation (choices(choice(redeclare
          DoubleSkinFacades.BasicComponents.WindowTypes.SingleGlass
          glassType), choice(redeclare
          DoubleSkinFacades.BasicComponents.WindowTypes.DoubleGlass
          glassType)), Placement(transformation(extent={{-66,-18},{-30,18}})));

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b annotation (
      Placement(transformation(extent={{90,-10},{110,10}}),iconTransformation(
          extent={{90,-10},{110,10}})));

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

  BoundaryHeatTransfer.InternalSurface insideSurface(H=H, L=L)
    annotation (Placement(transformation(extent={{6,-14},{42,14}})));
equation

  alpha1_dir = alpha1_dir_int;
  alpha2_dir = alpha2_dir_int;

  if cos(incidence_angle)>0.05 then s1=s/H*abs(tan(elevation));
  else s1=0;
  end if;
  if cos(HSA)>0.05 then s2=s*abs(tan(HSA))/L;
  else s2=0;
  end if;

  //glassType.shade_factor = 1-(1-s1)*(1-s2);
   glassType.shade_factor = 0;

  connect(port_a, glassType.port_a)
    annotation (Line(points={{-100,0},{-58.8,0}}, color={191,0,0}));
  connect(insideSurface.glassSurface, glassType.port_b) annotation (Line(points={{6,
          1.77636e-15},{-16,1.77636e-15},{-16,0},{-37.2,0}},
                                            color={191,0,0}));
  connect(insideSurface.airInterior, port_b)
    annotation (Line(points={{42,1.77636e-15},{64,1.77636e-15},{64,0},{100,
          0}},                                              color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-30,100},{30,-100}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid)}),                      Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end InternalWindow;

within DoubleSkinFacades.BasicComponents.AirChannelModels.VentilationType;
block ForzedVentilation
  "This component models the fluid-dynamic behaviour of the naturally ventalited air channel"
  extends PartialVentilation(            final useInputs = false);

  parameter Modelica.SIunits.Length s "Thickness of the air channel";
  parameter Modelica.SIunits.Area A_in "Area of air channel inlet";
  parameter Modelica.SIunits.Area A_out "Area of air channel outlet";
  parameter Modelica.SIunits.Height H "Height of the double-skin facade";
  // parameter Modelica.SIunits.DynamicViscosity mu = 2.27*10^(-5) "Dynamic viscosity of air";
  parameter Modelica.SIunits.Length Dh "Hydraulic diameter of air channel";

  Modelica.Blocks.Sources.Constant v     "Value of air velocity induced by forzed ventilation"
    annotation (Placement(transformation(extent={{-12,-10},{8,10}})));
equation
  connect(v.y, v_m)
    annotation (Line(points={{9,0},{37,0},{37,111}}, color={0,0,127}));
  annotation (Icon(graphics={Rectangle(extent={{-100,100},{100,-100}},
            lineColor={28,108,200}),
        Polygon(
          points={{-2,0},{-22,62},{0,80},{22,66},{
              -2,2},{-2,0}},
          lineColor={28,108,200},
          smooth=Smooth.Bezier,
          fillColor={0,255,255},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-2,-40},{-22,22},{0,40},{22,26},
              {-2,-38},{-2,-40}},
          lineColor={28,108,200},
          smooth=Smooth.Bezier,
          fillColor={0,255,255},
          fillPattern=FillPattern.Solid,
          origin={-40,4},
          rotation=90),
        Polygon(
          points={{-2,-40},{-22,22},{0,40},{22,26},
              {-2,-38},{-2,-40}},
          lineColor={28,108,200},
          smooth=Smooth.Bezier,
          fillColor={0,255,255},
          fillPattern=FillPattern.Solid,
          origin={-4,-36},
          rotation=180),
        Polygon(
          points={{-2,-40},{-22,22},{0,40},{22,26},
              {-2,-38},{-2,-40}},
          lineColor={28,108,200},
          smooth=Smooth.Bezier,
          fillColor={0,255,255},
          fillPattern=FillPattern.Solid,
          origin={36,0},
          rotation=270)}));
end ForzedVentilation;

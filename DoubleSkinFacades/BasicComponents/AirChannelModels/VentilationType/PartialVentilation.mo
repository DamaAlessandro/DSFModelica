within DoubleSkinFacades.BasicComponents.AirChannelModels.VentilationType;
partial block PartialVentilation
  "Partial for defining the number of occupants"
  extends Modelica.Blocks.Icons.Block;
  parameter Boolean useInputs = false
    "=true to use external input";

  Modelica.Blocks.Interfaces.RealInput Tf(start=300) if useInputs
  "Mean temperature of the air inside de channel" annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
        rotation=270,
        origin={-32,108}),            iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={-44,102})));
  Modelica.Blocks.Interfaces.RealOutput v_m(start = 0.2)
  "Mean velocity of the air inside the channel" annotation (Placement(
        transformation(extent={{-21,-21},{21,21}},
        rotation=90,
        origin={37,111}),                            iconTransformation(
        extent={{-21,-21},{21,21}},
        rotation=90,
        origin={45,101})));

  annotation (Icon(graphics,
                   coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
<ul>
<li>
March 28, 2019 by Filip Jorissen:<br/>
Added parameter <code>A</code> 
for <a href=\"https://github.com/open-ideas/IDEAS/issues/998\">#998</a>.
</li>
<li>
July 26, 2018 by Filip Jorissen:<br/>
First implementation
See <a href=\"https://github.com/open-ideas/IDEAS/issues/760\">#760</a>.
</li>
</ul>
</html>"));
end PartialVentilation;

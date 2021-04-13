within DoubleSkinFacades.BasicComponents.WindowTypes;
partial block PartialGlass "Partial for defining the type of glass"

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a
    annotation (Placement(transformation(extent={{-70,-10},{-50,10}}),
        iconTransformation(extent={{-70,-10},{-50,10}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b
    annotation (Placement(transformation(extent={{50,-10},{70,10}}),
        iconTransformation(extent={{50,-10},{70,10}})));

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end PartialGlass;

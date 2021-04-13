within DoubleSkinFacades.Examples;
model GlassPropertiesProba

  Modelica.Blocks.Tables.CombiTable1Ds GlassPropertiesDirect(
    tableOnFile=true,
    tableName="tab1",
    fileName="C:/Users/Jaime/OneDrive/TFM/GlassProperties.txt",
    columns=2:6,
    smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative)
    annotation (Placement(transformation(extent={{-10,56},{10,76}})));
  Modelica.Blocks.Sources.Clock clock
    annotation (Placement(transformation(extent={{-62,52},{-42,72}})));
equation
  connect(GlassPropertiesDirect.u, clock.y)
    annotation (Line(points={{-12,66},{-26,66},{-26,62},{-41,62}}, color={0,0,127}));
  annotation ();
end GlassPropertiesProba;

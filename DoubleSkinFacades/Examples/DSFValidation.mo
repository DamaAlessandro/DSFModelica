within DoubleSkinFacades.Examples;
model DSFValidation
  inner Real month = AtmosphericData.y[1];
  inner Real day = AtmosphericData.y[2];
  inner Real hour = AtmosphericData.y[3];
  inner Modelica.SIunits.Angle theta_wind = AtmosphericData.y[4] * (Modelica.Constants.pi/180.0) "Wind angle of attack";
  inner Modelica.SIunits.Velocity v_wind = AtmosphericData.y[5]*0.1 "Wind velocity";
  inner Modelica.SIunits.Irradiance Ghdiff = AtmosphericData.y[7] "Diffuse irradiance on the façade";
  inner Modelica.SIunits.Irradiance Ghtot = AtmosphericData.y[6] "Total irradiance on the façade";
  inner Modelica.SIunits.Temperature T_ext = AtmosphericData.y[8]*0.1 + 273.15 "External temperature";
  inner Modelica.SIunits.Pressure Patm = AtmosphericData.y[9] "Atmospheric pressure";
  inner Modelica.SIunits.Angle incidence_angle = partialDSF.solarAnglesCalculation.incidence "Angle of incidence of sun in radians";
  inner Modelica.SIunits.Temperature T_int = DSFValidation.y[1] + 273.15 "Internal temperature";
  Modelica.SIunits.Temperature Tei_exp = DSFValidation.y[5] +273.15;
  Modelica.SIunits.Temperature Tie_exp = DSFValidation.y[6] +273.15;
  Modelica.SIunits.Temperature Tii_exp = DSFValidation.y[7] +273.15;
  Modelica.SIunits.Temperature Tei_M = DSFValidation.y[2] +273.15;
  Modelica.SIunits.Temperature Tie_M = DSFValidation.y[3] +273.15;
  Modelica.SIunits.Temperature Tii_M = DSFValidation.y[4] +273.15;

  PartialDSF partialDSF(
    directPropertiesFile=
        "S:/backup/DottoratoPolimi/Conferenze/2021_IBPSA2021/ArticoloAle/DSF/GlassProperties.txt",
    epsilon_glass=0.84,
    epsilon_glazing=0.037,
    H=5.45,
    A_in=0.11,
    A_out=0.09,
    s=0.58,
    L=3.55,
    alpha1_diff_ext=1*0.213,
    alpha2_diff_ext=1*0.213,
    alpha1_diff_int=1*0.071,
    alpha2_diff_int=0.75*0.077,
    tau_diff_ext=0.69,
    tilt=90,
    orientation=0,
    latitude=57.05,
    longitude=9.93,
    timezone=1,
    airCavity(ventilationType(cpValuesFile=
            "C:/Users/damaa/Documents/Ricerca/DSF/DSF_modelica/DSF_Library/cpValues.txt")))
    annotation (Placement(transformation(extent={{-40,-40},{40,40}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-112,-10},{-92,10}})));
  Modelica.Blocks.Sources.Clock clock
    annotation (Placement(transformation(extent={{-86,62},{-66,82}})));
  Modelica.Blocks.Tables.CombiTable1Ds AtmosphericData(
    tableOnFile=true,
    tableName="tab1",
    fileName=
        "S:/backup/DottoratoPolimi/Conferenze/2021_IBPSA2021/ArticoloAle/DSF/BoundaryConditions10min.txt",
    columns=2:10,
    smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments)
    annotation (Placement(transformation(extent={{-38,62},{-18,82}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature1
    annotation (Placement(transformation(extent={{90,-10},{70,10}})));
  Modelica.Blocks.Tables.CombiTable1Ds DSFValidation(
    tableOnFile=true,
    tableName="tab1",
    fileName=
        "S:/backup/DottoratoPolimi/Conferenze/2021_IBPSA2021/ArticoloAle/DSF/DSFValidation10min.txt",
    columns=2:8,
    smoothness=Modelica.Blocks.Types.Smoothness.MonotoneContinuousDerivative1)
    annotation (Placement(transformation(extent={{2,46},{22,66}})));
equation
    prescribedTemperature.T = T_ext;
    prescribedTemperature1.T = T_int;

  connect(prescribedTemperature.port, partialDSF.port_a)
    annotation (Line(points={{-92,0},{-40,0}}, color={191,0,0}));
  connect(AtmosphericData.u, clock.y)
    annotation (Line(points={{-40,72},{-65,72}}, color={0,0,127}));
  connect(prescribedTemperature1.port, partialDSF.port_b)
    annotation (Line(points={{70,0},{56,0},{56,0},{40,0}}, color={191,0,0}));
  connect(DSFValidation.u, clock.y) annotation (Line(points={{0,56},{-50,56},{-50,
          72},{-65,72}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=1294800,
      Interval=600,
      Tolerance=1e-06,
      __Dymola_Algorithm="Esdirk45a"));
end DSFValidation;

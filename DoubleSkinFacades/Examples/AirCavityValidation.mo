within DoubleSkinFacades.Examples;
model AirCavityValidation
  inner Real month = AtmosphericData.y[1];
  inner Real day = AtmosphericData.y[2];
  inner Real hour = AtmosphericData.y[3];
  inner Modelica.SIunits.Angle theta_wind = AtmosphericData.y[4] * (Modelica.Constants.pi/180.0) "Wind angle of attack";
  inner Modelica.SIunits.Velocity v_wind = AtmosphericData.y[5]*0.1 "Wind velocity";
  inner Modelica.SIunits.Irradiance Ghtot = AtmosphericData.y[6] "Diffuse irradiance on the façade";
  inner Modelica.SIunits.Irradiance Ghdiff = AtmosphericData.y[7] "Total irradiance on the façade";
  inner Modelica.SIunits.Temperature T_ext = AtmosphericData.y[8]*0.1 + 273.15 "External temperature";
  inner Modelica.SIunits.Pressure Patm = AtmosphericData.y[9] "Atmospheric pressure";
  inner Modelica.SIunits.Angle incidence_angle = solarAnglesCalculation.incidence "Angle of incidence of sun in radians";
  inner Modelica.SIunits.PerUnit tau_dir_ext = GlassPropertiesDirect.y[5] "Direct transmitivity of exterior skin";
  inner Modelica.SIunits.PerUnit tau_dir = GlassPropertiesDirect.y[4];
  inner Modelica.SIunits.Irradiance Gdiff;
  inner Modelica.SIunits.Irradiance Gdir;
  Modelica.SIunits.Irradiance Gdsf;
  inner Modelica.SIunits.Angle incidence_angle_in;
  Modelica.SIunits.Irradiance Ghdir;
  parameter Modelica.SIunits.PerUnit albedo = 0.2;
  parameter Real orientation(final unit="deg") = 0;

  Modelica.SIunits.Temperature Tei = ExperimentalData.y[1]*0.1+273.15;
  Modelica.SIunits.Temperature Tie = ExperimentalData.y[2]*0.1+273.15;
  Modelica.SIunits.Temperature ToutletM = ValidationData.y[1]+273.15;
  Modelica.SIunits.Temperature ToutletE = ValidationData.y[2]+273.15;
  Modelica.SIunits.Velocity v_mM = ValidationData.y[3];
  Modelica.SIunits.Velocity v_mANE = ValidationData.y[4];
  Modelica.SIunits.Velocity v_mCO2 = ValidationData.y[5];
  Modelica.SIunits.HeatFlowRate qvM = ValidationData.y[6];
  Modelica.SIunits.HeatFlowRate qv_ANE = ValidationData.y[7];
  Modelica.SIunits.HeatFlowRate qv_CO2= ValidationData.y[8];

  inner Modelica.SIunits.Angle theta_windmin = AtmosphericData1.y[4] * (Modelica.Constants.pi/180.0) "Wind angle of attack";
  inner Modelica.SIunits.Velocity v_windmin = AtmosphericData1.y[5]*0.1 "Wind velocity";
  inner Modelica.SIunits.Irradiance Ghtotmin = AtmosphericData1.y[6] "Diffuse irradiance on the façade";
  inner Modelica.SIunits.Irradiance Ghdiffmin = AtmosphericData1.y[7] "Total irradiance on the façade";
  inner Modelica.SIunits.Temperature T_extmin = AtmosphericData1.y[8]*0.1 + 273.15 "External temperature";
  inner Modelica.SIunits.Pressure Patmmin = AtmosphericData1.y[9] "Atmospheric pressure";
  inner Modelica.SIunits.Irradiance Gdiffmin;
  inner Modelica.SIunits.Irradiance Gdirmin;
  Modelica.SIunits.Irradiance Gdsfmin;
  Modelica.SIunits.Irradiance Ghdirmin;

  Components.AirChannel airCavity(
    s=0.58,
    H=5.45,
    L=3.55,
    A_in=0.11,
    A_out=0.09,
    epsilon_glass=0.84,
    tau_diff_ext=0.69,
    orientation=orientation)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature1
    annotation (Placement(transformation(extent={{60,-10},{40,10}})));

  Modelica.Blocks.Sources.Clock clock
    annotation (Placement(transformation(extent={{-32,48},{-12,68}})));
  Modelica.Blocks.Tables.CombiTable1Ds ExperimentalData(
    tableOnFile=true,
    tableName="tab1",
    columns=2:3,
    smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
    fileName="C:/Users/Jaime/OneDrive/TFM/BoundaryConditions_AirCavityValidation.txt")
    annotation (Placement(transformation(extent={{16,48},{36,68}})));

  Modelica.Blocks.Tables.CombiTable1Ds AtmosphericData(
    tableOnFile=true,
    tableName="tab1",
    fileName="C:/Users/Jaime/OneDrive/TFM/ExperimentalData10min.txt",
    columns=2:10,
    smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments)
    annotation (Placement(transformation(extent={{16,78},{36,98}})));
  Modelica.Blocks.Tables.CombiTable1Ds GlassPropertiesDirect(
    tableOnFile=true,
    tableName="tab1",
    fileName="C:/Users/Jaime/OneDrive/TFM/GlassProperties.txt",
    columns=2:6,
    smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative)
    annotation (Placement(transformation(extent={{-20,-64},{0,-44}})));
  AuxiliaryBlocks.SolarAnglesCalculation solarAnglesCalculation(
    tilt=90,
    orientation=0,
    latitude=57.05,
    longitude=9.93,
    timezone=1)
    annotation (Placement(transformation(extent={{46,-66},{66,-46}})));

  Modelica.Blocks.Tables.CombiTable1Ds ValidationData(
    tableOnFile=true,
    tableName="tab1",
    columns=2:10,
    smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
    fileName="C:/Users/Jaime/OneDrive/TFM/ValidationData_AirCavity.txt")
    annotation (Placement(transformation(extent={{62,32},{82,52}})));
  Modelica.Blocks.Tables.CombiTable1Ds AtmosphericData1(
    tableOnFile=true,
    tableName="tab1",
    fileName="C:/Users/Jaime/OneDrive/TFM/ExperimentalData10min.txt",
    columns=2:10,
    smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments)
    annotation (Placement(transformation(extent={{52,76},{72,96}})));
equation

    prescribedTemperature.T = Tei;
    prescribedTemperature1.T = Tie;

    GlassPropertiesDirect.u = incidence_angle_in*(180.0/Modelica.Constants.pi);

 Gdsf = Gdir + Gdiff;
 Ghdir = max(0, Ghtot-Ghdiff);

 Gdsfmin = Gdirmin + Gdiffmin;
 Ghdirmin = max(0, Ghtotmin-Ghdiffmin);

  if Modelica.Math.cos(solarAnglesCalculation.zenith)>0.01 then
      Gdir = max(0, Ghdir*cos(incidence_angle)/cos(solarAnglesCalculation.zenith));
      Gdiff = (0.5*Ghdiff + albedo*0.5*Ghtot);
      incidence_angle_in=incidence_angle;
      Gdirmin = max(0, Ghdirmin*cos(incidence_angle)/cos(solarAnglesCalculation.zenith));
      Gdiffmin = (0.5*Ghdiffmin + albedo*0.5*Ghtotmin);
    else
      Gdir = 0;
      Gdiff = Ghtot;
      incidence_angle_in=Modelica.Constants.pi/2;
      Gdirmin = 0;
      Gdiffmin = Ghtotmin;
  end if;

  connect(prescribedTemperature1.port, airCavity.port_b)
    annotation (Line(points={{40,0},{10,0}}, color={191,0,0}));
  connect(prescribedTemperature.port, airCavity.port_a)
    annotation (Line(points={{-40,0},{-10,0}}, color={191,0,0}));
  connect(ExperimentalData.u, clock.y) annotation (Line(points={{14,58},{-11,58}}, color={0,0,127}));
  connect(AtmosphericData.u, clock.y)
    annotation (Line(points={{14,88},{14,73},{-11,73},{-11,58}}, color={0,0,127}));
  connect(ValidationData.u, clock.y) annotation (Line(points={{60,42},{0,42},{0,
          58},{-11,58}}, color={0,0,127}));
  connect(AtmosphericData1.u, clock.y) annotation (Line(points={{50,86},{20,86},
          {20,58},{-11,58}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(
          preserveAspectRatio=false)));
end AirCavityValidation;

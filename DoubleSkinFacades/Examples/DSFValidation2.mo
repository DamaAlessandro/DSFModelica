within DoubleSkinFacades.Examples;
model DSFValidation2
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
  Modelica.SIunits.Temperature Tei_exp = DSFValidation.y[5] +273.15;
  Modelica.SIunits.Temperature Tie_exp = DSFValidation.y[6] +273.15;
  Modelica.SIunits.Temperature Tii_exp = DSFValidation.y[7] +273.15;
  Modelica.SIunits.Temperature Tei_M = DSFValidation.y[2] +273.15;
  Modelica.SIunits.Temperature Tie_M = DSFValidation.y[3] +273.15;
  Modelica.SIunits.Temperature Tii_M = DSFValidation.y[4] +273.15;

  PartialDSF partialDSF(
    epsilon_glass=0.84,
    epsilon_glazing=0.037,
    H=5.45,
    A_in=0.11,
    A_out=0.09,
    s=0.58,
    L=3.55,
    alpha1_diff_ext=0.213,
    alpha2_diff_ext=0.213,
    alpha1_diff_int=0.071,
    alpha2_diff_int=0.077,
    tau_diff_ext=0.69,
    tilt=90,
    orientation=0,
    latitude=57.05,
    longitude=9.93,
    timezone=1)
    annotation (Placement(transformation(extent={{-40,-40},{40,40}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-112,-10},{-92,10}})));
  Modelica.Blocks.Sources.Clock clock
    annotation (Placement(transformation(extent={{-86,62},{-66,82}})));
  Modelica.Blocks.Tables.CombiTable1Ds AtmosphericData(
    tableOnFile=true,
    tableName="tab1",
    fileName="C:/Users/Jaime/OneDrive/TFM/ExperimentalData10min.txt",
    columns=2:10,
    smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments)
    annotation (Placement(transformation(extent={{-38,62},{-18,82}})));
  Modelica.Blocks.Tables.CombiTable1Ds DSFValidation(
    tableOnFile=true,
    tableName="tab1",
    fileName="C:/Users/Jaime/OneDrive/TFM/DSFValidation10min.txt",
    columns=2:8,
    smoothness=Modelica.Blocks.Types.Smoothness.MonotoneContinuousDerivative1)
    annotation (Placement(transformation(extent={{2,46},{22,66}})));
  IDEAS.Fluid.MixingVolumes.MixingVolume vol(
    redeclare package Medium = IDEAS.Media.Air,
    m_flow_nominal=20,
    V=143.11)
    annotation (Placement(transformation(extent={{70,-6},{90,14}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor Conduction(R=0.1095)
    annotation (Placement(transformation(extent={{-34,-84},{-14,-64}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor Comvection(R=1/(2*108))
    annotation (Placement(transformation(extent={{20,-84},{40,-64}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C=2656465.4)
    annotation (Placement(transformation(extent={{-6,-92},{14,-112}})));
  IDEAS.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad(
    redeclare package Medium = IDEAS.Media.Water,
    Q_flow_nominal=1000,
    T_a_nominal=353.15,
    T_b_nominal=333.15)
    annotation (Placement(transformation(extent={{-34,-170},{-14,-150}})));
  IDEAS.Fluid.HeatExchangers.HeaterCooler_u hea(
    redeclare package Medium = IDEAS.Media.Water,
    m_flow_nominal=0.2,
    Q_flow_nominal=1000,
    dp_nominal=1000)
    annotation (Placement(transformation(extent={{-80,-170},{-60,-150}})));
  IDEAS.Fluid.Movers.FlowControlled_m_flow fan(redeclare package Medium =
        IDEAS.Media.Water, m_flow_nominal=2) annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={-52,-196})));
  IDEAS.Fluid.Sources.Boundary_pT bou(nPorts=1, redeclare package Medium =
        IDEAS.Media.Water)
    annotation (Placement(transformation(extent={{-120,-206},{-100,-186}})));
  Modelica.Blocks.Logical.Not not1
    annotation (Placement(transformation(extent={{78,-170},{98,-150}})));
  Modelica.Blocks.Logical.Hysteresis hysteresis(
    pre_y_start=true,
    uLow=294.15,
    uHigh=296.15)
    annotation (Placement(transformation(extent={{46,-170},{66,-150}})));
  Modelica.Blocks.Math.BooleanToReal booleanToReal
    annotation (Placement(transformation(extent={{106,-170},{126,-150}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
    annotation (Placement(transformation(extent={{14,-150},{34,-130}})));
  Modelica.Blocks.Sources.Constant const(k=1)
    annotation (Placement(transformation(extent={{6,-200},{26,-180}})));
equation
    prescribedTemperature.T = T_ext;

  connect(prescribedTemperature.port, partialDSF.port_a)
    annotation (Line(points={{-92,0},{-40,0}}, color={191,0,0}));
  connect(AtmosphericData.u, clock.y)
    annotation (Line(points={{-40,72},{-65,72}}, color={0,0,127}));
  connect(DSFValidation.u, clock.y) annotation (Line(points={{0,56},{-50,56},{-50,
          72},{-65,72}}, color={0,0,127}));
  connect(vol.heatPort, partialDSF.port_b)
    annotation (Line(points={{70,4},{56,4},{56,0},{40,0}}, color={191,0,0}));
  connect(Conduction.port_a, partialDSF.port_a) annotation (Line(points={{-34,-74},
          {-72,-74},{-72,0},{-40,0}}, color={191,0,0}));
  connect(Conduction.port_b, Comvection.port_a)
    annotation (Line(points={{-14,-74},{20,-74}}, color={191,0,0}));
  connect(heatCapacitor.port, Comvection.port_a)
    annotation (Line(points={{4,-92},{4,-74},{20,-74}}, color={191,0,0}));
  connect(Comvection.port_b, vol.heatPort) annotation (Line(points={{40,-74},{54,
          -74},{54,4},{70,4}}, color={191,0,0}));
  connect(hea.port_b, rad.port_a)
    annotation (Line(points={{-60,-160},{-34,-160}}, color={0,127,255}));
  connect(fan.port_b, bou.ports[1])
    annotation (Line(points={{-62,-196},{-100,-196}}, color={0,127,255}));
  connect(fan.port_b, hea.port_a) annotation (Line(points={{-62,-196},{-88,-196},
          {-88,-160},{-80,-160}}, color={0,127,255}));
  connect(rad.port_b, fan.port_a) annotation (Line(points={{-14,-160},{-4,-160},
          {-4,-196},{-42,-196}}, color={0,127,255}));
  connect(temperatureSensor.T, hysteresis.u) annotation (Line(points={{34,-140},
          {36,-140},{36,-160},{44,-160}}, color={0,0,127}));
  connect(hysteresis.y, not1.u)
    annotation (Line(points={{67,-160},{76,-160}}, color={255,0,255}));
  connect(not1.y, booleanToReal.u)
    annotation (Line(points={{99,-160},{104,-160}}, color={255,0,255}));
  connect(booleanToReal.y, hea.u) annotation (Line(points={{127,-160},{134,-160},
          {134,-224},{-130,-224},{-130,-154},{-82,-154}}, color={0,0,127}));
  connect(const.y, fan.m_flow_in) annotation (Line(points={{27,-190},{-12,-190},
          {-12,-184},{-52,-184}}, color={0,0,127}));
  connect(temperatureSensor.port, partialDSF.port_b) annotation (Line(points={{14,
          -140},{6,-140},{6,-118},{62,-118},{62,4},{56,4},{56,0},{40,0}}, color=
         {191,0,0}));
  connect(rad.heatPortCon, partialDSF.port_b) annotation (Line(points={{-26,-152.8},
          {-28,-152.8},{-28,-126},{6,-126},{6,-118},{62,-118},{62,4},{56,4},{56,
          0},{40,0}}, color={191,0,0}));
  connect(rad.heatPortRad, Comvection.port_a) annotation (Line(points={{-22,-152.8},
          {-20,-152.8},{-20,-106},{-8,-106},{-8,-74},{20,-74}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-160,-240},
            {160,100}})),                                        Diagram(
        coordinateSystem(preserveAspectRatio=false, extent={{-160,-240},{160,100}})));
end DSFValidation2;

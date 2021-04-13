within DoubleSkinFacades.Examples;
model PartialDSF
  import DoubleSkinFacades;
  parameter Modelica.SIunits.Height H = 5.45 "Height of the façade" annotation (Dialog(group="Façade specifications"));
  parameter Modelica.SIunits.Height L = 3.55 "Length of the façade" annotation (Dialog(group="Façade specifications"));
  parameter Modelica.SIunits.Height s = 0.58 "Thickness of the air cavity" annotation (Dialog(group="Façade specifications"));
  parameter Modelica.SIunits.Area A_in "Area of air cavity inlet per meter of width" annotation (Dialog(group="Façade specifications"));
  parameter Modelica.SIunits.Area A_out "Area of air cavity outlet per meter of width" annotation (Dialog(group="Façade specifications"));
  parameter String directPropertiesFile = "C:/Users/Jaime/OneDrive/TFM/GlassProperties.txt" "File where optical direct properties matrix is stored" annotation (Dialog(group="Optical properties"));
  parameter Modelica.SIunits.PerUnit epsilon_glass "Emissivity of glass" annotation (Dialog(group="Optical properties"));
  parameter Modelica.SIunits.PerUnit epsilon_glazing "Emissivity of glazing" annotation (Dialog(group="Optical properties"));
  parameter Modelica.SIunits.PerUnit alpha1_diff_ext "Diffuse absorptivity of external glass in exterior skin" annotation (Dialog(group="Optical properties"));
  parameter Modelica.SIunits.PerUnit alpha2_diff_ext "Diffuse absorptivity of internal glass in exterior skin" annotation (Dialog(group="Optical properties"));
  parameter Modelica.SIunits.PerUnit alpha1_diff_int "Diffuse absorptivity of external glass in interior skin" annotation (Dialog(group="Optical properties"));
  parameter Modelica.SIunits.PerUnit alpha2_diff_int "Diffuse absorptivity of internal glass in interior skin" annotation (Dialog(group="Optical properties"));
  parameter Modelica.SIunits.PerUnit tau_diff_ext "Diffuse transmitivity of exterior skin" annotation (Dialog(group="Optical properties"));
  parameter Real tilt(final unit="deg") "Façade tilt angle in degrees" annotation (Dialog(group="Façade specifications"));
  parameter Real orientation(final unit="deg") "Façade orientation angle in degrees" annotation (Dialog(group="Façade specifications"));
  parameter Real latitude(final unit="deg") "Latitude in degrees" annotation (Dialog(group="Location info"));
  parameter Real longitude(final unit="deg") "Longitude in degrees" annotation(Dialog(group="Location info"));
  parameter Integer timezone "Time zone + GMT" annotation(Dialog(group="Location info"));
  inner Modelica.SIunits.PerUnit alpha1_dir_ext = D*GlassPropertiesDirect.y[1] "Direct absorptivity of external glass in exterior skin";
  inner Modelica.SIunits.PerUnit alpha2_dir_ext = D*GlassPropertiesDirect.y[1] "Direct absorptivity of internal glass in exterior skin";
  inner Modelica.SIunits.PerUnit alpha1_dir_int = D*GlassPropertiesDirect.y[2] "Direct absorptivity of external glass in interior skin";
  inner Modelica.SIunits.PerUnit alpha2_dir_int = D*GlassPropertiesDirect.y[3] "Direct absorptivity of internal glass in interior skin";
  inner Modelica.SIunits.PerUnit tau_dir_ext = GlassPropertiesDirect.y[5] "Direct transmitivity of exterior skin";
  inner Modelica.SIunits.PerUnit tau_dir = GlassPropertiesDirect.y[4];
  inner Modelica.SIunits.Irradiance Gdiff;
  inner Modelica.SIunits.Irradiance Gdir;
  inner Modelica.SIunits.Angle zenith = solarAnglesCalculation.zenith;
  inner Modelica.SIunits.Angle elevation = solarAnglesCalculation.elevation;
  inner Modelica.SIunits.Angle HSA = solarAnglesCalculation.HSA;
  outer Modelica.SIunits.Angle incidence_angle;
  outer Modelica.SIunits.Irradiance Ghdiff;
  outer Modelica.SIunits.Irradiance Ghtot;
  inner Modelica.SIunits.Angle incidence_angle_in;
  Modelica.SIunits.Irradiance Ghdir;
  parameter Modelica.SIunits.PerUnit albedo = 0.2;
  parameter Real D = 1;

  DoubleSkinFacades.Components.ExternalWindow exteriorSkin(
    H=H,
    L=L,
    material_glass=Media.Materials.Glass(),
    material_gas=DoubleSkinFacades.Media.Gases.Air(),
    alpha1_diff=D*alpha1_diff_ext,
    alpha2_diff=D*alpha2_diff_ext,
    epsilon_glass=epsilon_glass,
    epsilon_glazing=epsilon_glazing,
    redeclare DoubleSkinFacades.BasicComponents.WindowTypes.SingleGlass
      glassType,
    s_glass=0.008,
    s_gap=0,
    orientation=orientation)
    annotation (Placement(transformation(extent={{-52,-10},{-32,10}})));
  DoubleSkinFacades.Components.AirChannel airCavity(
    L=L,
    H=H,
    s=s,
    epsilon_glass=epsilon_glass,
    tau_diff_ext=tau_diff_ext,
    redeclare
      DoubleSkinFacades.BasicComponents.AirChannelModels.VentilationType.NaturalVentilation
      ventilationType(cpValuesFile=
          "C:/Users/damaa/Documents/Ricerca/DSF/DSF_modelica/DSF_Library/cpValues.txt"),
    A_in=A_in,
    A_out=A_out,
    orientation=orientation)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  DoubleSkinFacades.Components.InternalWindow interiorSkin(
    H=H,
    L=L,
    s=s,
    material_glass=Media.Materials.Glass(),
    alpha1_diff=D*alpha1_diff_int,
    alpha2_diff=D*alpha2_diff_int,
    epsilon_glass=epsilon_glass,
    epsilon_glazing=epsilon_glazing,
    redeclare DoubleSkinFacades.BasicComponents.WindowTypes.DoubleGlass
      glassType,
    s_glass=0.004,
    s_gap=0.016,
    material_gas=DoubleSkinFacades.Media.Gases.Argon())
    annotation (Placement(transformation(extent={{30,-10},{50,10}})));

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation (
      Placement(transformation(extent={{-110,-10},{-90,10}}),
        iconTransformation(extent={{-110,-10},{-90,10}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b annotation (
      Placement(transformation(extent={{90,-10},{110,10}}), iconTransformation(
          extent={{90,-10},{110,10}})));
  Modelica.Blocks.Tables.CombiTable1Ds GlassPropertiesDirect(
    tableOnFile=true,
    tableName="tab1",
    fileName=
        "S:/backup/DottoratoPolimi/Conferenze/2021_IBPSA2021/ArticoloAle/DSF/GlassProperties.txt",
    columns=2:6,
    smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative)
    annotation (Placement(transformation(extent={{-10,56},{10,76}})));
  DoubleSkinFacades.AuxiliaryBlocks.SolarAnglesCalculation
    solarAnglesCalculation(
    tilt=tilt,
    orientation=orientation,
    latitude=latitude,
    longitude=longitude,
    timezone=timezone)
    annotation (Placement(transformation(extent={{-10,-76},{10,-56}})));
equation

  GlassPropertiesDirect.u = incidence_angle_in*(180.0/Modelica.Constants.pi);

  //Calculation of direct and diffuse irradiance during sunrise and sunset
  Ghdir = max(0, Ghtot-Ghdiff);

  if Modelica.Math.cos(zenith)>0.05 then
      Gdir = 1*max(0, Ghdir*cos(incidence_angle)/cos(zenith));
      Gdiff = 1*(0.5*Ghdiff + albedo*0.5*Ghtot);
      incidence_angle_in=incidence_angle;
    else
      Gdir = 0;
      Gdiff = Ghtot;
      incidence_angle_in=Modelica.Constants.pi/2;
  end if;

  connect(exteriorSkin.port_b, airCavity.port_a)
    annotation (Line(points={{-33.6667,0},{-10,0}}, color={191,0,0}));

  connect(airCavity.port_b, interiorSkin.port_a)
    annotation (Line(points={{10,0},{30,0}}, color={191,0,0}));
  connect(exteriorSkin.port_a, port_a)
    annotation (Line(points={{-50.3333,0},{-100,0}}, color={191,0,0}));
  connect(interiorSkin.port_b, port_b)
    annotation (Line(points={{50,0},{100,0}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          fillColor={170,255,255},
          fillPattern=FillPattern.None),
        Rectangle(
          extent={{-60,80},{-28,-80}},
          lineColor={28,108,200},
          fillColor={216,255,251},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{28,80},{60,-80}},
          lineColor={28,108,200},
          fillColor={216,255,251},
          fillPattern=FillPattern.Solid)}),                      Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end PartialDSF;

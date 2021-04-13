within DoubleSkinFacades.AuxiliaryBlocks;
block SolarAnglesCalculation "Block for solar angles calulation"
  parameter Real tilt(final unit="deg") "Façade tilt angle in degrees";
  parameter Real orientation(final unit="deg") "Façade orientation angle in degrees";
  parameter Real latitude(final unit="deg") "Latitude in degrees";
  parameter Real longitude(final unit="deg") "Longitude in degrees";
  parameter Integer timezone "Time zone + GMT";
  outer Real month;
  outer Real day;
  outer Integer hour;
  Modelica.SIunits.Angle delta "Declination angle";
  Modelica.SIunits.Angle fi "Latitude";
  Modelica.SIunits.Angle beta "Surface tilt angle";
  Modelica.SIunits.Angle gamma "Façade azimuth angle";
  Modelica.SIunits.Angle azimuth "Solar azimuth angle";
  Modelica.SIunits.Angle azimuthfit "Solar azimuth angle";
  Modelica.SIunits.Angle HSA "Relative azimuth angle / Horizontal Shading Angle";
  //Modelica.SIunits.Angle altitude "Solar altitude";
  Modelica.SIunits.Angle VSA "Profile angle / Vertical Shading Angle";
  Modelica.SIunits.Angle elevation;
  Modelica.SIunits.Angle h "Hour angle + pi/6";
  Modelica.SIunits.Angle wss;
  Modelica.SIunits.Angle incidence "Incidence angle";
  Modelica.SIunits.Angle zenith "Zenith angle";
  Modelica.SIunits.Time timelong "Correction factor for longitude";
  Real nn "Number of day";
  Modelica.SIunits.Angle w "Solar hour angle";
  Real LSTM;
  Real EoT;
  Real B;
  Real TC;
  Real LST;
  Real LT;

protected
  Modelica.SIunits.PerUnit T;
  Modelica.SIunits.PerUnit U;
  Modelica.SIunits.PerUnit V;
  Modelica.SIunits.PerUnit Tz;
  Modelica.SIunits.PerUnit Uz;
  Modelica.SIunits.PerUnit Vz;
  Integer C1;
  Integer C2;
  Integer C3;

equation

  LSTM = 15 * timezone;
  EoT = 9.87 * sin(2*B) - 7.53 * cos(B) + 1.5 *sin(B);
  B = 360/365 * (nn - 81)*Modelica.Constants.pi/180;
  TC = 4*(longitude - LSTM) + EoT;
  LST = LT + TC/60;
  LT = hour;
  w = (LST - 12) * 15 * Modelica.Constants.pi/180;

  nn = (month-1)*30.5 + (day -1);
  delta = 23.45*sin((284+nn)*2*Modelica.Constants.pi/365)*Modelica.Constants.pi/180;

  fi = latitude*Modelica.Constants.pi/180;
  timelong = (longitude - timezone*15)/360*24;

  h = (hour-1 + 30/60 + timelong -12)/24*2*Modelica.Constants.pi;

  wss =acos(-tan(delta)*tan(fi));

  gamma = orientation*Modelica.Constants.pi/180;

  beta = tilt*Modelica.Constants.pi/180;

  T= sin(delta)*(sin(fi)*cos(beta)-cos(fi)*sin(beta)*cos(gamma));
  U= cos(delta)*(cos(fi)*cos(beta)+sin(fi)*sin(beta)*cos(gamma));
  V= cos(delta)*sin(beta)*sin(gamma);

  Tz= sin(delta)*sin(fi);
  Uz= cos(delta)*cos(fi);
  Vz= 0;

  incidence = (acos(T+U*cos(h)+V*sin(h)));
  zenith = (acos(Tz+Uz*cos(h)+Vz*sin(h)));
  elevation = Modelica.Constants.pi/2-zenith;

 if  abs(h)<= wss then     C1=1;
  else
    C1=-1;
  end if;

  if (fi-delta) >= 0 then     C2=1;
  else
    C2=-1;
  end if;

  if h >= 0 then     C3=1;
  else
    C3=-1;
  end if;

  azimuthfit = asin(sin(h)*cos(delta)/sin(zenith));
  azimuth= C1*C2*azimuthfit+C3*(1-C1*C2)/2*Modelica.Constants.pi;

  HSA = azimuth - gamma;
  //altitude = acos(cos(incidence)/cos(HSA));
  VSA = atan(tan(elevation)/cos(HSA));

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Ellipse(
          extent={{-40,40},{40,-40}},
          lineColor={244,125,35},
          fillColor={255,227,69},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-4,70},{2,46}},
          lineColor={244,125,35},
          fillColor={255,227,69},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-4,-48},{2,-72}},
          lineColor={244,125,35},
          fillColor={255,227,69},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-3,12},{3,-12}},
          lineColor={244,125,35},
          fillColor={255,227,69},
          fillPattern=FillPattern.Solid,
          origin={-59,0},
          rotation=90),
        Rectangle(
          extent={{-3,12},{3,-12}},
          lineColor={244,125,35},
          fillColor={255,227,69},
          fillPattern=FillPattern.Solid,
          origin={61,0},
          rotation=90),
        Polygon(
          points={{-58,54},{-38,34},{-34,38},{-54,58},{-58,54}},
          lineColor={244,125,35},
          fillColor={255,227,69},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{30,-38},{50,-58},{54,-54},{34,-34},{30,-38}},
          lineColor={244,125,35},
          fillColor={255,227,69},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-12,8},{8,-12},{12,-8},{-8,12},{-12,8}},
          lineColor={244,125,35},
          fillColor={255,227,69},
          fillPattern=FillPattern.Solid,
          origin={46,46},
          rotation=90),
        Polygon(
          points={{-12,8},{8,-12},{12,-8},{-8,12},{-12,8}},
          lineColor={244,125,35},
          fillColor={255,227,69},
          fillPattern=FillPattern.Solid,
          origin={-46,-46},
          rotation=90)}),                                        Diagram(coordinateSystem(
          preserveAspectRatio=false)));
end SolarAnglesCalculation;

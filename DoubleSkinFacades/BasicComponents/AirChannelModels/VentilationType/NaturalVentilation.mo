within DoubleSkinFacades.BasicComponents.AirChannelModels.VentilationType;
model NaturalVentilation
  "This component models the fluid-dynamic behaviour of the naturally ventilated air cavity"
  extends VentilationType.PartialVentilation(
                                         final useInputs = true);
  outer Modelica.SIunits.Temperature T_ext;
  outer Modelica.SIunits.Pressure Patm;
  outer Modelica.SIunits.Temperature v_wind;
  outer Modelica.SIunits.Angle theta_wind;
  parameter Modelica.SIunits.Height H "Height of the double-skin façade";
  parameter Modelica.SIunits.Length s "Thickness of the air cavity";
  parameter Modelica.SIunits.Area A_in "Area of air cavity inlet";
  parameter Modelica.SIunits.Area A_out "Area of air cavity outlet";
  parameter Modelica.SIunits.DynamicViscosity mu = 1.8*10^(-5) "Dynamic viscosity of air";
  parameter Modelica.SIunits.Acceleration g = 9.806 "Value of gravity";
  parameter Modelica.SIunits.SpecificHeatCapacity Rstar = 8314/29 "Ideal gas R constant value of air";
  parameter Modelica.SIunits.Length Dh = 2*s "Hydraulic diameter of air cavity";
  parameter Modelica.SIunits.PerUnit epsilon = 0.01 "Superficial roughness";
  parameter String cpValuesFile = "C:/Users/damaa/Documents/Ricerca/DSF/DSF modelica/DoubleSkinFaçades Library/CpValues.txt" "File where pressure coefficient values matrix is stored" annotation (Dialog(group="Optical properties"));
  final parameter Modelica.SIunits.PerUnit Cd_top = 1 / (1/0.6 - A_out/s) "Discharge coefficient of top opening";
  final parameter Modelica.SIunits.PerUnit Cd_bottom = 1 / (1/0.6 - A_in/s) "Discharge coeffiecient of bottom opening";
  parameter Modelica.SIunits.PerUnit kappa = 1 "Factor that takes into account the discontinuity of wind";
  parameter Modelica.SIunits.PerUnit alpha = 0.2 "Wind obstacle coefficient alpha=0.15 rural area - alpha=0.4 highly built area";
  parameter Modelica.SIunits.Height h_ref= 10 "Reference height for wind measurements";
  parameter Real orientation(final unit="deg") "Façade orientation angle in degrees" annotation (Dialog(group="Façade specifications"));
  Modelica.SIunits.PerUnit Cp_top "Wind pressure drop coefficient at the top opening";
  Modelica.SIunits.PerUnit Cp_bottom "Wind pressure drop coefficient at the bottom opening";
  Modelica.SIunits.PerUnit Re(start=1) "Reynolds number of air inside the cavity";
  Modelica.SIunits.Pressure deltaP_wind "Pressure drop due to the incident wind";
  Modelica.SIunits.Pressure deltaP_thermal "Pressure drop due to buoyancy inside the cavity";
  //Modelica.SIunits.Pressure deltaP "Pressure drop";
  Modelica.SIunits.Pressure deltaP_tot(start=1*10^(-5)) "Total pressure drop inside the cavity";
  Modelica.SIunits.Density rho_m(start=1.225) "Mean density of air insdide the channel";
  Modelica.SIunits.Density rho_ext "Mean density of air outside";
  Modelica.SIunits.Velocity v_roof "Velocity of outside air at H roof";
  Modelica.SIunits.Temperature Tf_in;
  Real xi(start=1) "Colebrook and White parameter";
  Integer upwards;
  Real senseflow;
  Real senseflow1;

  Modelica.Blocks.Tables.CombiTable1Ds cpValuesTable(
    smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
    tableOnFile=true,
    tableName="tab1",
    fileName=
        "S:/backup/DottoratoPolimi/Conferenze/2021_IBPSA2021/ArticoloAle/DSF/CpValues.txt",
    columns=2:3)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
equation

  Tf_in = Tf;

  cpValuesTable.u = abs(theta_wind * 180 / (Modelica.Constants.pi) -orientation- 180);
  Cp_top = cpValuesTable.y[2];
  Cp_bottom = cpValuesTable.y[1];

  rho_m = Patm/(Rstar*Tf_in);
  rho_ext = Patm/(Rstar*T_ext);

  v_roof = v_wind * kappa * (H / h_ref)^alpha;

  deltaP_thermal = g*H*(rho_ext - rho_m);
  deltaP_wind = 0.5 * rho_ext * v_roof^2 * abs(Cp_bottom - Cp_top);
  deltaP_tot = deltaP_wind + deltaP_thermal;
  senseflow = g*H*(rho_ext - rho_m)+0.5 * rho_ext * v_roof^2 * (Cp_bottom - Cp_top);
  senseflow1 = 0.5 * rho_ext * v_roof^2 * (Cp_bottom - Cp_top);

  if senseflow > 0 then upwards =1;
  else upwards =0;
  end if;

  Re = rho_m*abs(v_m)*Dh/mu;

  xi=abs((-2*log((epsilon/Dh)/3.7065) + 2.5226/(Re*sqrt(abs(xi))))^(-2));

  v_m = sqrt(abs(2*abs(deltaP_tot)/rho_m/((s/(A_in*Cd_bottom))^2 + xi*H/Dh + (s/(A_out*Cd_top))^2)));

  annotation (Icon(graphics={Rectangle(extent={{-100,100},{100,-100}},
            lineColor={28,108,200}),
        Polygon(
          points={{-68,52},{-68,52}},
          lineColor={255,255,85},
          fillColor={255,255,85},
          fillPattern=FillPattern.Solid),
        Line(points={{-62,38},{22,38},{36,52},{24,66}},
            color={28,108,200},smooth=Smooth.Bezier),
        Line(points={{-68,-30},{16,-30},{26,-48},{14,-60}},
            color={28,108,200},smooth=Smooth.Bezier),
        Line(points={{-58,-10},{48,-10},{58,-40},{48,-58}},
            color={28,108,200},smooth=Smooth.Bezier),
        Line(points={{-56,18},{52,14},{72,40},{54,68}},
            color={28,108,200},smooth=Smooth.Bezier)}));
end NaturalVentilation;

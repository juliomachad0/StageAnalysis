from stage_analysis import EfficiencyStage
import numpy as np
import matplotlib.pyplot as plt

tomega = 1000000
bartopa = 100000
patobar = 1 / 100000

# ENGINE PROP AND TANK PARAMETERS
Pc = 60 * bartopa                       # Pa
tb = 1000                               # s, 4 min
ts = 10e-3                              # s
of = 1.64
PoGasPressureFed = 200 * bartopa        # Gas tank initial pressure of fed system
PoGasPumpFed = 200 * bartopa            # Gas tank initial pressure of pump systems

# PROPELLANTS PROPERTIES
ox_properties = [878, 'l']              # [density, phase]
fuel_properties = [1440, 'l']           # [density, phase]

# [Massa molar, To, specific ratio gases of camera]
gas_properties = [4.0026 / 1000, 288.15, 1.23]
gas_tank_properties = [1700, 3300 * tomega]
ox_tank_properties = [2800, 455 * tomega]
fuel_tank_properties = [2800, 455 * tomega]

# PUMP TURBINE AND GG PARAMETERS
T_itu = 900  # K temperatura da turbina
P_turb_ratio = 20  # relcao Pentrada/Psaida tubeira
turbina = [T_itu, P_turb_ratio]

# pump and turbine power density
p_pd_ox = 40e3  # pump power density ox W/kg
p_pd_fuel = 40e3  # pump power density fuel
t_powerd = 20e3  # turbine power density
pumpPowerDensity = [p_pd_ox, p_pd_fuel]

# material gas generator properties
gg_densitymaterial = 8890
gg_uts = 524e6  # pa
gg_mat = [gg_densitymaterial, gg_uts]
# gg gas properties
gg_densitygas = 0.8774  # kg/m³
M_gg = 14.2e-3  # molar mass of gg kg/mol
gamma_gg = 1.196  # gamma gg
gg_gas_prop = [gg_densitygas, M_gg, gamma_gg, Pc]
# efficiencys
n_pu_ox = 0.8
n_pu_fuel = 0.8
n_turbine = 0.8
efficiencys = [n_pu_ox, n_pu_fuel, n_turbine]
pump_eff = [n_pu_ox, n_pu_fuel]
# ELECTRIC PUMP, ELECTRIC ENGINE AND INVERTER PARAMETERS
delta_motor = 3.8e3 # W/kg
delta_inv = 60e3    # W/kg
delta_bap = 4.8e3   # W/kg
delta_baw = 130*3600     # Wh/kg
n_ee = 0.8          # rendimento motor eletrico
n_inv = 0.85

teste = EfficiencyStage.StageAnalysis(Pc=Pc, of=of, tb=tb,
                                     ox_properties=ox_properties,
                                     fuel_properties=fuel_properties,
                                     gas_properties=gas_properties,
                                     gas_tank_properties=gas_tank_properties,
                                     ox_tank_properties=ox_tank_properties,
                                     fuel_tank_properties=fuel_tank_properties)
teste.Set_PropellantTanksPressure()
teste.Set_Po_Gas(condition='calculated', PoGasPressureFed=PoGasPressureFed,
                 PoGasPumpFed=PoGasPumpFed,
                 PoGasPressureFed_to_Pc=10, PoGasPumpFed_to_Pc=1.5)
teste.Set_PumpParameters(PowerDensity=pumpPowerDensity, efficiencies=pump_eff)
teste.Set_ElectricCycle_parameters(motor_power_d=delta_motor, inverter_power_d=delta_inv,
                                   battery_power_d_bap=delta_bap,
                                   battery_power_d_baw=delta_baw,
                                   motor_n=n_ee, inverter_n=n_inv)
teste.Set_OpenCycle_parameters(ts=ts,
                               turbine_power_density=t_powerd,
                               turbine_efficiency=n_turbine,
                               turbine_pressure_ratio=P_turb_ratio,
                               turbine_temp=T_itu,
                               gg_mat=gg_mat, gg_gas_prop=gg_gas_prop)

teste.Update_mp(1000)
pressure_vector = np.arange(10 * bartopa, 1000 * bartopa, 0.5 * bartopa)
prop_vector = np.arange(10, 1000, 10)
po_gas_vector = np.arange(130 * bartopa, 300 * bartopa, 0.1 * bartopa)
time_vector = np.arange(10, 1000, 1)
xlegenda = 'pressure'
x_vector = pressure_vector
MR_fed = []
SE_fed = []
MR_pump = []
SE_pump = []
SE_electric = []
massa_estrutura = []
massa_gas = []
for i in x_vector:
    teste.Update_Pc(i)
    SE_fed.append(teste.StructuralEfficiency(fed_system='pressure_fed'))
    SE_pump.append(teste.StructuralEfficiency(fed_system='open_cycle'))
    SE_electric.append(teste.StructuralEfficiency(fed_system='electric_pump'))

if xlegenda == 'pressure':
    x_vector = [i*patobar for i in x_vector]

plt.plot(x_vector, SE_pump, color='red', label='open cycle')
plt.plot(x_vector, SE_fed, color='green', label='pressure fed')
plt.plot(x_vector, SE_electric, color='blue', label='electric pump')
plt.xlabel(xlegenda)
plt.ylabel('Structural Efficiency - mp/(mp+me)')
plt.legend()
#plt.savefig('C:/Users/julio.machado/Desktop/steste2.png', format='png', dpi=300)
plt.show()

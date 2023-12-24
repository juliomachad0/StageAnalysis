from stage_analysis.EficiencyStage import StageAnalysis
import numpy as np
import matplotlib.pyplot as plt

tomega = 1000000
bartopa = 100000
patobar = 1/100000
# ENGINE PROP AND TANK PARAMETERS
Pc = 30 * bartopa
tb = 240  # s
ts = 10e-3  # s
of = 1.64
PoGas = [200 * bartopa, 200 * bartopa]  # Gas tank initial pressure
ku = 1.2

# PROPELLANTS PROPERTIES
ox_properties = [878, 'l']  # [densidade, fase]
fuel_properties = [1440, 'l']
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
PowerDensity = [p_pd_ox, p_pd_fuel, t_powerd]

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
# ELECTRIC PUMP, ELECTRIC ENGINE AND INVERTER PARAMETERS


teste = StageAnalysis(Pc=Pc, of=of, tb=tb, PoGas=PoGas,
                      ku=ku, ox_properties=ox_properties,
                      fuel_properties=fuel_properties,
                      gas_properties=gas_properties,
                      gas_tank_properties=gas_tank_properties,
                      ox_tank_properties=ox_tank_properties,
                      fuel_tank_properties=fuel_tank_properties)
teste.Set_TPU_parameters(ts=ts, turbine=turbina, PowerDensity=PowerDensity,
                         gg_mat=gg_mat, gg_gas_prop=gg_gas_prop,
                         efficiencys=efficiencys)
teste.Update_mp(1000)

pressure_vector = np.arange(0.5 * tomega, 5 * tomega, 0.5 * tomega)
prop_vector = np.arange(10, 1000, 10)
po_gas_vector = np.arange(130 * bartopa, 300 * bartopa, 0.1 * bartopa)

x_vector = pressure_vector
MR_fed = []
SE_fed = []
MR_pump = []
SE_pump = []

for i in x_vector:
    teste.Update_Pc(i)
    SE_fed.append(teste.MassPropRatio())
    SE_pump.append(teste.MassPropRatio(feed_system='open_cycle'))

plt.plot(x_vector, SE_pump, color='red', label='open cycle')
plt.plot(x_vector, SE_fed, color='green', label='pressure fed')
plt.title('Eficiencia estrutural\n Pressao: 30 bar, PressaGas: 200 bar, Btime: 240 s')
plt.xlabel('propelente (kg)')
plt.legend()
plt.show()

#  MR_fed.append(teste.MassPropRatio())
#  MR_pump.append(teste.MassPropRatio(feed_system='open_cycle'))
'''
x_vector = pressure_vector
MR_fed = []
SE_fed = []
MR_pump = []
SE_pump = []

for i in x_vector:
    teste.Update_Pc(i)

    SE_fed.append(teste.StructuralEfficiency())
    SE_pump.append(teste.StructuralEfficiency(feed_system='open_cycle'))

plt.plot([i/100000 for i in x_vector], SE_pump, color='red', label='open cycle')
plt.plot([i/100000 for i in x_vector], SE_fed, color='green', label='pressure fed')
plt.title('Eficiencia estrutural\n Pressao: 30 bar, PressaGas: 200 bar, Btime: 240 s')
plt.xlabel('pressao (bar)')
plt.legend()
plt.show()

teste.Update_Pc(Pc=Pc)
teste.Update_mp(1000)
x_vector = po_gas_vector
MR_fed = []
SE_fed = []
MR_pump = []
SE_pump = []

for i in x_vector:
    teste.Update_Pc(i)

    SE_fed.append(teste.StructuralEfficiency())
    SE_pump.append(teste.StructuralEfficiency(feed_system='open_cycle'))

plt.plot([i/100000 for i in x_vector], SE_pump, color='red', label='open cycle')
plt.plot([i/100000 for i in x_vector], SE_fed, color='green', label='pressure fed')
plt.title('Eficiencia estrutural\n Pressao: 30 bar, PressaGas: 200 bar, Btime: 240 s')
plt.legend()
plt.show()
'''
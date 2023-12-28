import numpy as np


class StageAnalysisAndDesign:
    def __init__(self, Pc: float, of: float,
                 tb: float, r_ext_stage: float,
                 ox_properties: list,
                 fuel_properties: list, gas_properties: list,
                 gas_tank_properties: list, ox_tank_properties: list,
                 fuel_tank_properties: list):
        # constants
        self.Runi = 8.31    # J/molK
        self.g_o = 9.80665  # m/s2
        # PoGAS configuration
        self.PoGasPumpFed_to_Pc = None      # Po Gas/Pc ratio of the pressure fed system
        self.PoGasPressureFed_to_Pc = None  # Po Gas/Pc ratio of the pump fed systems
        self.PoGasPressureFed = None        # value of the Pressure fed system PoGas
        self.PoGasPumpFed = None            # value of the Pump fed system Po gas
        self.PoGas_condition = None         # condition to estimate Po Gas
        self.PoGas_to_Ptank = 1.5          # in case of PoGas lower than Ptank
        # safety constants
        self.ku = 1.05      # safe constant, volume V_tank = ku*V_prop
        self.kg = 1.3       # safe constant, gas mass, mass the left in the gas tank
        self.ktg = 2.5      # safe constant, mass of gas tank
        self.ktp = 1.25     # safe constant, mass of propellant tanks
        self.k_gg = 2.5     # safe constant, gas generator
        self.kb = 1.2       # safe constant, battery mass safety
        self.FsTanks = 1.2  # safety factor Fs = Te/Ta
        # ENGINE PARAMETERS
        self.Pc = Pc         # chamber pressure, Pascal
        self.of = of         # mixture ratio
        self.burn_time = tb  # burning time
        self.ts = None       # stay time
        # STAGE PARAMETERS
        self.r_ext_stage = r_ext_stage
        # ISP parameters
        self.v_exit_tub = None
        self.DP_exit_tub = None
        self.A_exit_tub = np.pi * (self.r_ext_stage ** 2)
        self.I_sp = None
        # system total mass
        self.m_tp = None  # massa total do sistema de turbo bomba
        self.m_ep = None  # massa total do sistema de bomba eletrica
        self.m_pg = None  # massa total do sistema de pressurização
        # Massa dos componentes
        self.m_g = None   # massa do gas pressurizante
        self.m_tg = None  # massa do tanque de gas
        self.m_to = None  # massa do tanque de oxidante
        self.m_tf = None  # massa do tanque de combustivel
        # pumps, gg, turbine
        self.m_pu = None   # massa das bombas
        self.m_tu = None   # massa das turbinas
        self.m_gg = None   # massa do gerador de gas
        self.m_ptu = None  # turbine driven propellant mass
        # electric
        self.m_ee = None   # massa do motor eletrico
        self.m_inv = None  # massa do inversor
        self.m_bat = None  # massa das baterias
        # massa do propelente
        self.mp = None     # massa de propelente
        # propellants tanks
        self.PropellantTankPressure_condition = 'calculated'
        self.OxPropellantTankPressure = [0, 0, 0]    # list: [PressureFed, OpenCycle, Electric_Pump]
        self.FuelPropellantTankPressure = [0, 0, 0]  # list: [PressureFed, OpenCycle, Electric_Pump]
        self.PropellantTankPressure_method = 'ratios'
        self.Ratio_PressureFedPtankOx_Pc = 1.8
        self.Ratio_PressureFedPtankFuel_Pc = 1.8
        self.Ratio_PumpFedPtankOx_Pc = 0.3
        self.Ratio_PumpFedPtankFuel_Pc = 0.3
        # propellant and gas properties
        self.alpha = None
        self.alpha_ox = None
        self.alpha_fuel = None
        self.ox_properties = ox_properties  # [densidade]
        self.fuel_properties = fuel_properties
        self.pre_gas_molar_mass = gas_properties[0]  # [Massa molar, To, specific ratio]
        self.pre_gas_temp = gas_properties[1]
        self.pre_gas_gamma = gas_properties[2]
        # propellant and gas tanks properties
        self.gas_tank_properties = gas_tank_properties    # [density, tensile strength]
        self.ox_tank_properties = ox_tank_properties      # density, tensile strength
        self.fuel_tank_properties = fuel_tank_properties  # density, tensile strength
        # PUMPS, TURBINES, GG PARAMETERS
        self.gg_densitymaterial = None
        self.gg_uts = None
        self.pump_power_density_ox = None
        self.pump_power_density_fuel = None
        self.turbine_power_densi = None
        self.gg_densitygases = None
        self.M_gg = None                      # massa molar do gas gg
        self.gamma_gg = None                  # gamma gg
        self.P_gg = None                      # pressao na camara do gg
        self.T_in_tu = None                   # temperatura da turbina
        self.turbine_pressure_ratio = None
        self.n_pu_ox = None                   # eficiencia da bomba de oxidante
        self.n_pu_fuel = None                 # eficiencia da bombda de combustivel
        self.n_turbine = None                 # eficiencia da turbina
        # ELECTRIC PUMP SYSTEM FED
        self.motor_power_density = None
        self.inverter_power_density = None
        self.battery_power_density_bap = None
        self.battery_power_density_baw = None
        self.motor_efficiency = None
        self.inverter_efficiency = None

    def Set_Po_Gas(self, condition="calculated", PoGasPressureFed=None, PoGasPumpFed=None,
                   PoGasPressureFed_to_Pc=None, PoGasPumpFed_to_Pc=None,
                   PoGas_to_Ptank=None):
        if PoGas_to_Ptank is not None:
            self.PoGas_to_Ptank = PoGas_to_Ptank
        if condition.lower() == "informed":
            self.PoGasPressureFed = PoGasPressureFed
            self.PoGasPumpFed = PoGasPumpFed
            self.PoGas_condition = 'informed'
            return
        elif condition.lower() == "calculated":
            if PoGasPressureFed_to_Pc is None or PoGasPumpFed_to_Pc is None:
                print("Stage Analysis:Set_Po_Gas: A least one of the PoGas/Pc ratios not provided")
                if PoGasPressureFed_to_Pc is None:
                    self.PoGasPressureFed_to_Pc = 10
                if PoGasPumpFed_to_Pc is None:
                    self.PoGasPumpFed_to_Pc = 1.5
                print("Set_Po_Gas: PoGas/Pc Ratios used: PoGasPressureFed/Pc = {}, PoGasPumpFedToPc = {}".format(
                    self.PoGasPressureFed_to_Pc,
                    self.PoGasPumpFed_to_Pc))
                self.PoGasPressureFed = self.PoGasPressureFed_to_Pc * self.Pc
                self.PoGasPumpFed = self.PoGasPumpFed_to_Pc * self.Pc
                self.PoGas_condition = 'calculated'
                return
            self.PoGasPressureFed_to_Pc = PoGasPressureFed_to_Pc
            self.PoGasPumpFed_to_Pc = PoGasPumpFed_to_Pc
            self.PoGasPressureFed = self.PoGasPumpFed_to_Pc * self.Pc
            self.PoGasPumpFed = self.PoGasPumpFed_to_Pc * self.Pc
            self.PoGas_condition = 'calculated'
            return
        else:
            print("Stage Analysis:Set_Po_Gas: wrong input parameter informed on "
                  "'condition', options: 'informed','calculated'")
            print("Stage Analysis:Set_Po_Gas: PoGas_to_Pc ratios used: Pressure fed: 10, Pump fed: 1.5")
            self.PoGasPressureFed_to_Pc = 10
            self.PoGasPumpFed_to_Pc = 1.5
            self.PoGasPressureFed = self.PoGasPressureFed_to_Pc*self.Pc
            self.PoGasPumpFed = self.PoGasPumpFed_to_Pc*self.Pc
            self.PoGas_condition = 'calculated'

    def Set_PropellantTanksPressure(self, OxPropellantTankPressure=None, FuelPropellantTankPressure=None,
                                    condition='calculated',
                                    method='ratios',
                                    Ratio_PressureFedPtankOx_Pc=1.8,
                                    Ratio_PressureFedPtankFuel_Pc=1.8,
                                    Ratio_PumpFedPtankOx_Pc=0.3,
                                    Ratio_PumpFedPtankFuel_Pc=0.3):
        if OxPropellantTankPressure is None:
            self.OxPropellantTankPressure = [0, 0, 0]
        else:
            self.OxPropellantTankPressure = OxPropellantTankPressure
        if FuelPropellantTankPressure is None:
            self.FuelPropellantTankPressure = [0, 0, 0]
        else:
            self.FuelPropellantTankPressure = FuelPropellantTankPressure

        self.PropellantTankPressure_condition = condition
        self.PropellantTankPressure_method = method
        self.Ratio_PressureFedPtankOx_Pc = Ratio_PressureFedPtankOx_Pc
        self.Ratio_PressureFedPtankFuel_Pc = Ratio_PressureFedPtankFuel_Pc
        self.Ratio_PumpFedPtankOx_Pc = Ratio_PumpFedPtankOx_Pc
        self.Ratio_PumpFedPtankFuel_Pc = Ratio_PumpFedPtankFuel_Pc

    def Set_PumpParameters(self, PowerDensity: list, efficiencies: list):
        self.n_pu_ox = efficiencies[0]  # oxidizer pump efficiency
        self.n_pu_fuel = efficiencies[1]  # fuel pump efficiency
        self.pump_power_density_ox = PowerDensity[0]  # propellant pump power density  - ox
        self.pump_power_density_fuel = PowerDensity[1]  # propellant pump power density - fuel

    def Set_OpenCycle_parameters(self, ts,
                                 turbine_temp: float,
                                 turbine_pressure_ratio: float,
                                 turbine_power_density: float,
                                 turbine_efficiency: float,
                                 gg_mat: list, gg_gas_prop: list):
        self.ts = ts  # tempo de permanencia
        # TURBINE PARAMETERS
        self.n_turbine = turbine_efficiency  # eficiencia da turbina
        self.T_in_tu = turbine_temp  # temperatura da turbina
        self.turbine_pressure_ratio = turbine_pressure_ratio  # taxa de P_entrada/P_saida turbina
        self.turbine_power_densi = turbine_power_density  # turbine power densit
        # GAS GENERATOR PARAMETERS
        # material
        self.gg_densitymaterial = gg_mat[0]  # densidade do material
        self.gg_uts = gg_mat[1]  # tensao de escoamento
        # gas product of combustion
        self.gg_densitygases = gg_gas_prop[0]  # densidade do gas gerador
        self.M_gg = gg_gas_prop[1]  # massa molar do gerador de gas
        self.gamma_gg = gg_gas_prop[2]  # gamma do gerador de gas
        self.P_gg = gg_gas_prop[3]  # gas generator pressure

    def Set_ElectricCycle_parameters(self, motor_power_d, inverter_power_d,
                                     battery_power_d_bap, battery_power_d_baw,
                                     motor_n, inverter_n):
        self.motor_power_density = motor_power_d
        self.inverter_power_density = inverter_power_d
        self.battery_power_density_bap = battery_power_d_bap
        self.battery_power_density_baw = battery_power_d_baw
        self.motor_efficiency = motor_n
        self.inverter_efficiency = inverter_n

    def Set_Isp_Parameters(self, v_exit_tub, DP_exit_tub=0, r_exit_tub=None):
        if r_exit_tub is not None:
            self.A_exit_tub = np.pi*r_exit_tub**2
        self.v_exit_tub = v_exit_tub
        self.DP_exit_tub = DP_exit_tub

    def Isp(self, condition='sea_level'):
        if condition == 'sea_level':
            self.I_sp = self.v_exit_tub / self.g_o
            return self.I_sp
        elif condition == 'vacuum':
            m_dot = self.mp/self.burn_time
            self.I_sp = (m_dot * self.v_exit_tub + self.DP_exit_tub * self.A_exit_tub) / (m_dot * self.g_o)
            return self.I_sp

    def AlphaFunc(self):  # calcula os alfs citados no doc com base na taxa de mistura
        self.alpha_ox = (self.of / self.ox_properties[0]) * (1 / (1 + self.of))
        self.alpha_fuel = (1 / self.fuel_properties[0]) * (1 / (1 + self.of))
        self.alpha = self.alpha_fuel + self.alpha_ox

    # PROPELLANT AND GAS MASSES AND PRESSURE
    def PropellantTankPressure(self, fed_system, propellant_type=None):
        if self.PropellantTankPressure_condition == 'calculated':
            # based on ratios
            if self.PropellantTankPressure_method == 'ratios':
                if propellant_type == 'ox':
                    if fed_system == 'pressure_fed':
                        self.OxPropellantTankPressure[0] = self.Ratio_PressureFedPtankOx_Pc*self.Pc
                        return self.OxPropellantTankPressure[0]
                    elif fed_system == 'open_cycle':
                        self.OxPropellantTankPressure[1] = self.Ratio_PumpFedPtankOx_Pc*self.Pc
                        return self.OxPropellantTankPressure[1]
                    elif fed_system == 'electric_pump':
                        self.OxPropellantTankPressure[2] = self.Ratio_PumpFedPtankOx_Pc*self.Pc
                        return self.OxPropellantTankPressure[2]
                elif propellant_type == 'fuel':
                    if fed_system == 'pressure_fed':
                        self.FuelPropellantTankPressure[0] = self.Ratio_PressureFedPtankFuel_Pc * self.Pc
                        return self.FuelPropellantTankPressure[0]
                    elif fed_system == 'open_cycle':
                        self.FuelPropellantTankPressure[1] = self.Ratio_PumpFedPtankFuel_Pc * self.Pc
                        return self.FuelPropellantTankPressure[1]
                    elif fed_system == 'electric_pump':
                        self.FuelPropellantTankPressure[2] = self.Ratio_PumpFedPtankFuel_Pc * self.Pc
                        return self.FuelPropellantTankPressure[2]
            # base on kessaev aproximation
            elif self.PropellantTankPressure_method == 'kessaev':
                DeltaPcOx = self.DeltaP(self.Pc, self.ox_properties[1])
                DeltaPcFuel = self.DeltaP(self.Pc, self.fuel_properties[1])
                if propellant_type == 'ox':
                    self.OxPropellantTankPressure = [self.Pc + DeltaPcOx, self.Pc + DeltaPcOx, self.Pc + DeltaPcOx]
                    return self.OxPropellantTankPressure[0]

                elif propellant_type == 'fuel':
                    self.FuelPropellantTankPressure = [self.Pc + DeltaPcFuel, self.Pc + DeltaPcFuel,
                                                       self.Pc + DeltaPcFuel]
                    return self.FuelPropellantTankPressure[0]

        elif self.PropellantTankPressure_condition == 'informed':
            if fed_system == 'pressure_fed':
                if propellant_type == 'ox':
                    return self.OxPropellantTankPressure[0]
                if propellant_type == 'fuel':
                    return self.FuelPropellantTankPressure[0]
            elif fed_system == 'open_cycle':
                if propellant_type == 'ox':
                    return self.OxPropellantTankPressure[1]
                if propellant_type == 'fuel':
                    return self.FuelPropellantTankPressure[1]
            elif fed_system == 'electric_pump':
                if propellant_type == 'ox':
                    return self.OxPropellantTankPressure[2]
                if propellant_type == 'fuel':
                    return self.FuelPropellantTankPressure[2]

    def PoGas(self, condition, fed_system):
        # open cycle and electric pump system - turbo pump systems
        if fed_system == 'open_cycle' or fed_system == 'electric_pump':
            if condition == 'informed':
                return self.PoGasPumpFed
            elif condition == 'calculated':
                self.PoGasPumpFed = self.PoGasPumpFed_to_Pc*self.Pc
                prop_tank_ox_pressure = self.PropellantTankPressure(fed_system=fed_system, propellant_type='ox')
                prop_tank_fuel_pressure = self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel')
                tank_pressure = self.Compare('max', prop_tank_ox_pressure, prop_tank_fuel_pressure)
                if self.PoGasPumpFed <= tank_pressure:
                    self.PoGasPumpFed = self.PoGas_to_Ptank*tank_pressure
                    print("PoGasPressureFed lower than Ptank-> PoGas/Ptank = {}".format(self.PoGas_to_Ptank))
                return self.PoGasPumpFed
        # pressure fed system
        elif fed_system == 'pressure_fed':
            if condition == 'informed':
                return self.PoGasPressureFed
            elif condition == 'calculated':
                self.PoGasPressureFed = self.PoGasPressureFed_to_Pc * self.Pc
                prop_tank_ox_pressure = self.PropellantTankPressure(fed_system=fed_system, propellant_type='ox')
                prop_tank_fuel_pressure = self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel')
                tank_pressure = self.Compare('max', prop_tank_ox_pressure, prop_tank_fuel_pressure)
                if self.PoGasPressureFed <= tank_pressure:
                    self.PoGasPressureFed = self.PoGas_to_Ptank*tank_pressure
                    print("PoGasPressureFed lower than Ptank-> PoGas/Ptank = {}".format(self.PoGas_to_Ptank))
                return self.PoGasPressureFed

    def GasMass(self, kp, PoGas):
        self.AlphaFunc()
        termo1 = (self.kg * kp * self.Pc * self.pre_gas_molar_mass)/(self.pre_gas_temp*self.Runi)
        termo2 = (self.ku * self.alpha * self.mp)**self.pre_gas_gamma
        termo3 = 1 - kp * (self.Pc / PoGas)
        termo4 = 1/termo3
        self.m_g = termo1 * termo2 * termo4
        return self.m_g

    def GasTankMass(self, kp, PoGas):
        termo1 = (3 * self.gas_tank_properties[0]) / (2 * self.gas_tank_properties[1])
        self.AlphaFunc()
        termo2 = self.ktg * self.kg * kp * self.ku * self.pre_gas_gamma * self.alpha
        termo3 = 1 - kp * (self.Pc / PoGas)
        termo4 = (self.mp * self.Pc) / termo3
        self.m_tg = termo1 * termo2 * termo4
        return self.m_tg

    def OxTankMass(self, kp):
        termo1 = (3 * self.ox_tank_properties[0]) / (2 * self.ox_tank_properties[1])
        termo2 = self.ktp * kp * self.ku * self.alpha_ox * self.mp * self.Pc
        self.m_to = termo1 * termo2
        return self.m_to

    def FuelTankMass(self, kp):
        termo1 = (3 * self.fuel_tank_properties[0]) / (2 * self.fuel_tank_properties[1])
        termo2 = self.ktp * kp * self.ku * self.alpha_fuel * self.mp * self.Pc
        self.m_to = termo1 * termo2
        return self.m_to

    def PumpsMass(self, fed_system):
        Dpox = self.DeltaP(self.Pc, self.ox_properties[1])
        Dpfuel = self.DeltaP(self.Pc, self.fuel_properties[1])

        Dpuox = self.Pc + Dpox - self.PropellantTankPressure(fed_system=fed_system, propellant_type='ox')
        Dpufuel = self.Pc + Dpfuel - self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel')

        kpiox = Dpuox / self.Pc
        kpifuel = Dpufuel / self.Pc

        kpox = self.PropellantTankPressure(fed_system=fed_system, propellant_type='ox') / self.Pc
        kpfuel = self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel') / self.Pc

        termo1ox = (1 + kpiox - kpox) * (self.alpha_ox * self.Pc * self.mp)
        pump_mass_ox = termo1ox / (self.pump_power_density_ox * self.burn_time)

        termo1fuel = (1 + kpifuel - kpfuel) * (self.alpha_fuel * self.Pc * self.mp)
        pump_mass_fuel = termo1fuel / (self.pump_power_density_fuel * self.burn_time)

        self.m_pu = pump_mass_ox + pump_mass_fuel
        return self.m_pu

    def TurbineMass(self, fed_system):
        kpox = self.PropellantTankPressure(fed_system=fed_system, propellant_type='ox') / self.Pc
        kpfuel = self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel') / self.Pc
        kp = self.Compare(condition='max', a=kpox, b=kpfuel)
        Dpfuel = self.DeltaP(self.Pc, self.fuel_properties[1])
        Dpufuel = self.Pc + Dpfuel - self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel')
        kpifuel = Dpufuel / self.Pc
        termo1 = (1 + kpifuel - kp)
        termo2 = (self.Pc * self.mp) / self.burn_time
        termo3 = self.alpha_fuel / self.n_pu_fuel + self.alpha_ox / self.n_pu_ox
        termo4 = 1 / self.turbine_power_densi
        self.m_tu = termo1 * termo2 * termo3 * termo4
        return self.m_tu

    def PropMassbyTurbine(self, fed_system):
        kpox = self.PropellantTankPressure(fed_system=fed_system, propellant_type='ox') / self.Pc
        kpfuel = self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel') / self.Pc
        kp = self.Compare(condition='max', a=kpox, b=kpfuel)
        Dpfuel = self.DeltaP(self.Pc, self.fuel_properties[1])
        Dpufuel = self.Pc + Dpfuel - self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel')
        kpifuel = Dpufuel / self.Pc
        termo1 = (1 + kpifuel - kp)
        termo2 = (self.Pc * self.mp) / self.n_turbine
        termo3 = self.alpha_fuel / self.n_pu_fuel + self.alpha_ox / self.n_pu_ox
        termo4 = self.M_gg * (self.gamma_gg - 1) / (self.Runi * self.gamma_gg)
        termo51 = (self.gamma_gg - 1) / self.gamma_gg
        termo52 = 1 - (1 / self.turbine_pressure_ratio) ** termo51
        termo5 = 1 / (self.T_in_tu * termo52)
        self.m_ptu = termo1 * termo2 * termo3 * termo4 * termo5
        return termo5

    def GGMas(self, Pgg=None):
        if Pgg is None:
            Pgg = self.Pc
        if self.m_ptu is None:
            return "Massa consumida pela turbina deve ser calculada primeiro"
        vgg = (self.ts * self.m_ptu) / (self.burn_time * self.gg_densitygases)
        termo1 = 3 * self.gg_densitymaterial / (2 * self.gg_uts)
        termo2 = self.k_gg * Pgg
        self.m_gg = termo1 * termo2 * vgg
        return self.m_gg

    def Electric_EngineMass(self, kp, kpi):
        termo1 = (1+kpi-kp)*((self.Pc*self.mp)/self.burn_time)
        termo2 = self.alpha_ox/self.n_pu_ox + self.alpha_fuel/self.n_pu_fuel
        termo3 = 1/self.motor_power_density
        self.m_ee = termo1*termo2*termo3
        return self.m_ee

    def Inverter_Mass(self, kp, kpi):
        termo1 = (1 + kpi - kp) * ((self.Pc * self.mp) / self.burn_time)
        termo2 = self.alpha_ox / self.n_pu_ox + self.alpha_fuel / self.n_pu_fuel
        termo3 = 1 / (self.inverter_power_density*self.motor_efficiency)
        self.m_inv = termo1 * termo2 * termo3
        return self.m_inv

    def Battery_Mass(self, kp, kpi):
        termo1 = 1 + kpi - kp
        termo2 = self.Pc * self.mp
        termo3 = self.alpha_ox / self.n_pu_ox + self.alpha_fuel / self.n_pu_fuel
        termo4 = self.kb/(self.inverter_efficiency*self.motor_efficiency)
        mbap = termo1*termo2*(1/self.burn_time)*termo3*termo4*(1/self.battery_power_density_bap)
        mbaw = termo1*termo2*termo3*termo4*(1/self.battery_power_density_baw)
        mba = self.Compare(condition='max', a=mbap, b=mbaw)
        self.m_bat = mba
        return self.m_bat

    def PressureFeedSystemMass(self, PoGas, kpox, kpfuel):
        kpmax = self.Compare(condition='max', a=kpox, b=kpfuel)
        gas_system_mass = self.GasMass(kpmax, PoGas) + self.GasTankMass(kpmax, PoGas)
        prop_tank_mass = self.OxTankMass(kpox) + self.FuelTankMass(kpox)
        self.m_pg = gas_system_mass + prop_tank_mass
        return self.m_pg

    def OpenCycleGGSystemMass(self):
        system_type = 'open_cycle'
        kpox = self.PropellantTankPressure(fed_system=system_type, propellant_type='ox') / self.Pc
        kpfuel = self.PropellantTankPressure(fed_system=system_type, propellant_type='fuel') / self.Pc
        PoGas = self.PoGas(condition=self.PoGas_condition, fed_system=system_type)
        self.m_pg = self.PressureFeedSystemMass(PoGas=PoGas, kpox=kpox, kpfuel=kpfuel)
        self.m_pu = self.PumpsMass(fed_system=system_type)  # massa das bombas
        self.m_tu = self.TurbineMass(fed_system=system_type)  # massa das turbinas
        self.m_ptu = self.PropMassbyTurbine(fed_system=system_type)  # massa do gerador de gas
        self.m_gg = self.GGMas(Pgg=self.P_gg)  # turbine driven propellant mass
        # turbo pump system mass
        self.m_tp = self.m_pg + self.m_pu + self.m_tu + self.m_ptu + self.m_gg
        return self.m_tp

    def ElectricPumpSystemMass(self, fed_system):
        system_type = 'open_cycle'
        kpox = self.PropellantTankPressure(fed_system=system_type, propellant_type='ox') / self.Pc
        kpfuel = self.PropellantTankPressure(fed_system=system_type, propellant_type='fuel') / self.Pc
        PoGas = self.PoGas(condition=self.PoGas_condition, fed_system=system_type)
        self.m_pg = self.PressureFeedSystemMass(PoGas=PoGas, kpox=kpox, kpfuel=kpfuel)
        self.m_pu = self.PumpsMass(fed_system=system_type)  # massa das bombas

        Dpox = self.DeltaP(self.Pc, self.ox_properties[1])
        Dpfuel = self.DeltaP(self.Pc, self.fuel_properties[1])
        Dpuox = self.Pc + Dpox - self.PropellantTankPressure(fed_system=fed_system, propellant_type='ox')
        Dpufuel = self.Pc + Dpfuel - self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel')
        kpiox = Dpuox / self.Pc
        kpifuel = Dpufuel / self.Pc

        kpi = self.Compare(condition='max', a=kpiox, b=kpifuel)
        kp = self.Compare(condition='max', a=kpox, b=kpfuel)

        self.m_ee = self.Electric_EngineMass(kp=kp, kpi=kpi)
        self.m_inv = self.Inverter_Mass(kp=kp, kpi=kpi)
        self.m_bat = self.Battery_Mass(kp=kp, kpi=kpi)
        self.m_ep = self.m_pg + self.m_ee + self.m_inv + self.m_bat + self.m_pu
        return self.m_ep

    def StructuralEfficiency(self, fed_system):
        match fed_system:
            case 'pressure_fed':
                kpox = self.PropellantTankPressure(fed_system=fed_system, propellant_type='ox') / self.Pc
                kpfuel = self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel') / self.Pc
                PoGas = self.PoGas(condition=self.PoGas_condition, fed_system=fed_system)
                if kpox*self.Pc >= PoGas or kpfuel*self.Pc >= PoGas:
                    return 0
                self.m_pg = self.PressureFeedSystemMass(PoGas=PoGas, kpox=kpox, kpfuel=kpfuel)
                return self.mp / (self.mp + self.m_pg)
            case 'electric_pump':
                kpox = self.PropellantTankPressure(fed_system=fed_system, propellant_type='ox') / self.Pc
                kpfuel = self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel') / self.Pc
                PoGas = self.PoGas(condition=self.PoGas_condition, fed_system=fed_system)
                if kpox*self.Pc >= PoGas or kpfuel*self.Pc >= PoGas:
                    return 0
                self.m_ep = self.ElectricPumpSystemMass(fed_system=fed_system)
                return self.mp / (self.mp + self.m_ep)

            case 'open_cycle':
                kpox = self.PropellantTankPressure(fed_system=fed_system, propellant_type='ox') / self.Pc
                kpfuel = self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel') / self.Pc
                PoGas = self.PoGas(condition=self.PoGas_condition, fed_system=fed_system)
                if kpox * self.Pc >= PoGas or kpfuel * self.Pc >= PoGas:
                    return 0
                self.m_tp = self.OpenCycleGGSystemMass()
                return self.mp / (self.mp + self.m_tp)
            case _:
                return "Wrong input parameter. options: 'pressure_feed', 'electric_pump' 'open_cycle'"

    def MassRatio(self, fed_system):
        match fed_system:
            case 'pressure_fed':
                kpox = self.PropellantTankPressure(fed_system=fed_system, propellant_type='ox') / self.Pc
                kpfuel = self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel') / self.Pc
                PoGas = self.PoGas(condition=self.PoGas_condition, fed_system=fed_system)
                if kpox * self.Pc >= PoGas or kpfuel * self.Pc >= PoGas:
                    return 0
                self.m_pg = self.PressureFeedSystemMass(PoGas=PoGas, kpox=kpox, kpfuel=kpfuel)
                return self.m_pg / self.mp
            case 'electric_pump':
                kpox = self.PropellantTankPressure(fed_system=fed_system, propellant_type='ox') / self.Pc
                kpfuel = self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel') / self.Pc
                PoGas = self.PoGas(condition=self.PoGas_condition, fed_system=fed_system)
                if kpox * self.Pc >= PoGas or kpfuel * self.Pc >= PoGas:
                    return 0
                self.m_ep = self.ElectricPumpSystemMass(fed_system=fed_system)
                return self.m_ep / self.mp

            case 'open_cycle':
                kpox = self.PropellantTankPressure(fed_system=fed_system, propellant_type='ox') / self.Pc
                kpfuel = self.PropellantTankPressure(fed_system=fed_system, propellant_type='fuel') / self.Pc
                PoGas = self.PoGas(condition=self.PoGas_condition, fed_system=fed_system)
                if kpox * self.Pc >= PoGas or kpfuel * self.Pc >= PoGas:
                    return 0
                self.m_tp = self.OpenCycleGGSystemMass()
                return self.m_tp / self.mp
            case _:
                return "Wrong input parameter. options: 'pressure_feed', 'electric_pump' 'open_cycle'"

    def Update_mp(self, mp):
        self.mp = mp

    def Update_Pc(self, Pc):
        self.Pc = Pc

    def Update_Tb(self, Tb):
        self.burn_time = Tb

    def Update_PoGas(self, PoGasPressureFed, PoGasPumpFed):
        self.PoGasPressureFed = PoGasPressureFed
        self.PoGasPumpFed = PoGasPumpFed

    def Infos(self):
        print("Pc: {:.4f} bar".format(self.Pc / 100000))
        print("Massa de Propelente: {:.4f} kg".format(self.mp))
        print("Massa pressure Fed: {:.4f} kg".format(self.PressureFeedSystemMass()))
        print("Structural Efficiency: {:.4f} ".format(self.StructuralEfficiency()))
        print("Massa sistema/ Massa de Propelente: {:.4f} ".format(
            self.PressureFeedSystemMass() / self.mp
        ))

    @staticmethod
    def DeltaP(Pc, phase='l'):
        if phase == 'l':
            dp = 80 * np.sqrt(10 * Pc)
            if dp < 400000:
                return 400000
            else:
                return dp
        elif phase == 'g':
            dp = 40 * np.sqrt(10 * Pc)
            if dp < 400000:
                return 400000
            else:
                return dp
        elif phase == 'h':
            dp = 60 * np.sqrt(10 * Pc)
            if dp < 400000:
                return 400000
            else:
                return dp

    @staticmethod
    def Compare(condition, a, b):
        v = [a, b]
        if condition == 'max':
            return max(v)
        if condition == 'min':
            return min(v)

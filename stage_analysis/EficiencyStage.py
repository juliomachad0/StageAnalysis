import numpy as np


class StageAnalysis:
    def __init__(self, Pc, of, tb, PoGas, ku, ox_properties: list,
                 fuel_properties: list, gas_properties: list,
                 gas_tank_properties: list, ox_tank_properties: list,
                 fuel_tank_properties: list):
        # constants
        self.Runi = 8.31  # J/molK

        # safety constants
        self.ku = 1.05  # safe constant, volume V_tank = ku*V_prop
        self.kg = 1.3  # safe constant, gas mass, mass the left in the gas tank
        self.ktg = 2.5  # safe constant, mass of gas tank
        self.ktp = 1.25  # safe constant, mass of propellant tanks
        self.k_gg = 2.5  # safe constant, gas generator

        self.Pc = Pc
        self.of = of
        self.burn_time = tb
        self.ts = None  # tempo de permanencia

        # massa total dos sistemas
        self.m_tp = None  # massa total do sistema de turbo bomba
        self.m_ep = None  # massa total do sistema de bomba eletrica
        self.m_pg = None  # massa total do sistema de pressurização
        # Massa dos componentes
        self.m_g = None  # massa do gas pressurizante
        self.m_tg = None  # massa do tanque de gas
        self.m_to = None  # massa do tanque de oxidante
        self.m_tf = None  # massa do tanque de combustivel
        # pumps, gg, turbine
        self.m_pu = None  # massa das bombas
        self.m_tu = None  # massa das turbinas
        self.m_gg = None  # massa do gerador de gas
        self.m_ptu = None  # turbine driven propellant mass
        # electric
        self.m_ee = None  # massa do motor eletrico
        self.m_inv = None  # massa do inversor

        self.m_bat = None  # massa das baterias

        # massa do propelente
        self.mp = None  # massa de propelente
        # propellant and gas properties
        self.PoGas = PoGas
        self.alpha = None
        self.alpha_ox = None
        self.alpha_fuel = None
        self.ox_properties = ox_properties  # [densidade]
        self.fuel_properties = fuel_properties
        self.gas_properties = gas_properties  # [Massa molar, To, specific ratio]
        # propellant and gas tanks properties
        self.gas_tank_properties = gas_tank_properties
        self.ox_tank_properties = ox_tank_properties
        self.fuel_tank_properties = fuel_tank_properties
        self.ku = ku  # extra volume, v_tank/v_propellant

        # PUMPS, TURBINES, GG PARAMETERS
        self.gg_densitymaterial = None
        self.gg_uts = None

        self.pump_powerd_ox = None
        self.pump_powerd_fuel = None
        self.turbine_power_densi = None

        self.gg_densitygases = None
        self.M_gg = None  # massa molar do gas gg
        self.gamma_gg = None  # gamma gg
        self.P_gg = None  # pressao na camara do gg

        self.T_in_tu = None  # temperatura da turbina
        self.Pratio_turb = None

        self.n_pu_ox = None  # eficiencia da bomba de oxidante
        self.n_pu_fuel = None  # eficiencia da bombda de combustivel
        self.n_turbine = None  # eficiencia da turbina

    def AlphaFunc(self): # calcula os alfs citados no doc com base na taxa de mistura
        self.alpha_ox = (self.of / self.ox_properties[0]) * (1 / (1 + self.of))
        self.alpha_fuel = (1 / self.fuel_properties[0]) * (1 / (1 + self.of))
        self.alpha = self.alpha_fuel + self.alpha_ox

    # PROPELLANT AND GAS MASSES AND PRESSURE
    def PropellantTankPressure(self, system_type):
        if system_type == 'pressure_fed':
            return self.Pc * 1.8
        elif system_type == 'pump':
            return self.Pc * 0.3

    def GasMass(self, kp, PoGas):
        self.AlphaFunc()
        termo1 = kp * self.ku * self.gas_properties[2] * self.alpha
        termo2 = (self.gas_properties[0] / (self.Runi * self.gas_properties[1]))
        termo3 = (self.mp * self.Pc) / (1 - kp * (self.Pc / PoGas))
        self.m_g = self.kg * termo1 * termo2 * termo3
        return self.m_g

    def GasTankMass(self, kp, PoGas):
        termo1 = (3 * self.gas_tank_properties[0]) / (2 * self.gas_tank_properties[1])
        self.AlphaFunc()
        termo2 = self.ktg * self.kg * kp * self.ku * self.gas_properties[2] * self.alpha
        termo3 = (self.mp * self.Pc) / (1 - kp * (self.Pc / PoGas))
        self.m_tg = termo1 * termo2 * termo3
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

    # OPEN CYCLE FUNCTIONS, TURBINE, PUMP, GG, TURBINEDRIVE

    def Set_TPU_parameters(self, ts, turbine, PowerDensity,
                           gg_mat: list, gg_gas_prop: list,
                           efficiencys: list, ):
        self.ts = ts  # tempo de permanencia

        self.T_in_tu = turbine[0]  # temperatura da turbina
        self.Pratio_turb = turbine[1]  # taxa de P_entrada/P_saida turbina

        self.pump_powerd_ox = PowerDensity[0]  # propellant pump power density  - ox
        self.pump_powerd_fuel = PowerDensity[1]  # propellant pump power density - fuel
        self.turbine_power_densi = PowerDensity[2]  # turbine power density

        self.gg_densitymaterial = gg_mat[0]  # densidade do material
        self.gg_uts = gg_mat[1]  # tensao de escoamento

        self.gg_densitygases = gg_gas_prop[0]  # densidade do gas gerador
        self.M_gg = gg_gas_prop[1]  # massa molar do gerador de gas
        self.gamma_gg = gg_gas_prop[2]  # gamma do gerador de gas
        self.P_gg = gg_gas_prop[3]
        self.n_pu_ox = efficiencys[0]  # eficiencia da bomba de oxidante
        self.n_pu_fuel = efficiencys[1]  # eficiencia da bombda de combustivel
        self.n_turbine = efficiencys[2]  # eficiencia da turbina

    def PumpsMass(self):
        Dpox = self.DeltaP(self.Pc, self.ox_properties[1])
        Dpfuel = self.DeltaP(self.Pc, self.fuel_properties[1])

        Dpuox = self.Pc + Dpox - self.PropellantTankPressure(system_type='pump')
        Dpufuel = self.Pc + Dpfuel - self.PropellantTankPressure(system_type='pump')

        kpiox = Dpuox / self.Pc
        kpifuel = Dpufuel / self.Pc

        kp = self.PropellantTankPressure(system_type='pump') / self.Pc
        termo1ox = (1 + kpiox - kp) * (self.alpha_ox * self.Pc * self.mp)
        pump_mass_ox = termo1ox / (self.pump_powerd_ox * self.burn_time)

        termo1fuel = (1 + kpifuel - kp) * (self.alpha_fuel * self.Pc * self.mp)
        pump_mass_fuel = termo1fuel / (self.pump_powerd_fuel * self.burn_time)

        self.m_pu = pump_mass_ox + pump_mass_fuel
        return self.m_pu

    def TurbineMass(self):
        kp = self.PropellantTankPressure(system_type='pump') / self.Pc
        Dpfuel = self.DeltaP(self.Pc, self.fuel_properties[1])
        Dpufuel = self.Pc + Dpfuel - self.PropellantTankPressure(system_type='pump')
        kpifuel = Dpufuel / self.Pc
        termo1 = (1 + kpifuel - kp)
        termo2 = (self.Pc * self.mp) / self.burn_time
        termo3 = self.alpha_fuel / self.n_pu_fuel + self.alpha_ox / self.n_pu_ox
        termo4 = 1 / self.turbine_power_densi
        self.m_tu = termo1 * termo2 * termo3 * termo4
        return self.m_tu

    def PropMassbyTurbine(self):
        kp = self.PropellantTankPressure(system_type='pump') / self.Pc
        Dpfuel = self.DeltaP(self.Pc, self.fuel_properties[1])
        Dpufuel = self.Pc + Dpfuel - self.PropellantTankPressure(system_type='pump')
        kpifuel = Dpufuel / self.Pc
        termo1 = (1 + kpifuel - kp)
        termo2 = (self.Pc * self.mp) / self.n_turbine
        termo3 = self.alpha_fuel / self.n_pu_fuel + self.alpha_ox / self.n_pu_ox
        termo4 = self.M_gg * (self.gamma_gg - 1) / (self.Runi * self.gamma_gg)
        termo51 = (self.gamma_gg - 1) / self.gamma_gg
        termo52 = 1 - (1 / self.Pratio_turb) ** termo51
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

    def PressureFeedSystemMass(self):
        system_type = 'pressure_fed'
        kp = self.PropellantTankPressure(system_type=system_type) / self.Pc
        gas_system_mass = self.GasMass(kp, self.PoGas[0]) + self.GasTankMass(kp, self.PoGas[0])
        prop_tank_mass = self.OxTankMass(kp) + self.FuelTankMass(kp)
        self.m_pg = gas_system_mass + prop_tank_mass
        return self.m_pg

    def ElectricPumpSystemMass(self):
        pass

    def OpenCycleGGSystemMass(self):
        system_type = 'pump'
        kp = self.PropellantTankPressure(system_type=system_type) / self.Pc
        gas_system_mass = self.GasMass(kp, self.PoGas[1]) + self.GasTankMass(kp, self.PoGas[1])
        prop_tank_mass = self.OxTankMass(kp) + self.FuelTankMass(kp)
        self.m_pg = gas_system_mass + prop_tank_mass
        self.m_pu = self.PumpsMass()  # massa das bombas
        self.m_tu = self.TurbineMass()  # massa das turbinas
        self.m_ptu = self.PropMassbyTurbine()  # massa do gerador de gas
        self.m_gg = self.GGMas(Pgg=self.P_gg)  # turbine driven propellant mass
        # turbo pump system mass
        self.m_tp = self.m_pg + self.m_pu + self.m_tu + self.m_ptu + self.m_gg
        return self.m_tp

    def StructuralEfficiency(self, feed_system='pressure_feed'):
        match feed_system:
            case 'pressure_feed':
                self.m_pg = self.PressureFeedSystemMass()
                return self.mp / (self.mp + self.m_pg)
            case 'electric_pump':
                return self.mp / (self.mp + self.m_ep)
            case 'open_cycle':
                self.m_tp = self.OpenCycleGGSystemMass()
                return self.mp / (self.mp + self.m_tp)
            case _:
                return "Wrong input parameter. options: 'pressure_feed', 'electric_pump' 'open_cycle'"

    def MassPropRatio(self, feed_system='pressure_feed'):
        match feed_system:
            case 'pressure_feed':
                self.m_pg = self.PressureFeedSystemMass()
                return self.m_pg / self.mp
            case 'electric_pump':
                return self.m_ep / self.mp
            case 'open_cycle':
                self.m_tp = self.OpenCycleGGSystemMass()
                return self.m_tp / self.mp
            case _:
                return "Wrong input parameter. options: 'pressure_feed', 'electric_pump' 'open_cycle'"

    def Update_mp(self, mp):
        self.mp = mp

    def Update_Pc(self, Pc):
        self.Pc = Pc

    def Update_PoGas(self, PoGas):
        self.PoGas = PoGas

    def Infos(self):
        print("Pc: {:.4f} bar".format(self.Pc / 100000))
        print("Massa de Propelente: {:.4f} kg".format(self.mp))
        print("Massa pressure Fed: {:.4f} kg".format(self.PressureFeedSystemMass()))
        print("Eficiencia estrutural: {:.4f} ".format(self.StructuralEfficiency()))
        print("Massa sistema/ Massa de Propelente: {:.4f} ".format(
            self.PressureFeedSystemMass() / self.mp
        ))

    @staticmethod
    def DeltaP(Pc, phase='l'):
        if phase == 'l':
            return 1.3 * 80 * np.sqrt(10 * Pc)
        elif phase == 'g':
            return 1.3 * 40 * np.sqrt(10 * Pc)
        elif phase == 'h':
            return 1.3 * 60 * np.sqrt(10 * Pc)

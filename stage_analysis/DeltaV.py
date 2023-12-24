import numpy as np
import pandas as pd


class AnalysisOfDeltaV:
    def __init__(self, stages: list, jettison_mass: list = None,
                 payload_mass: float = 0.0, name: str = 'Rocket'):
        self.payload_mass = payload_mass  # kg
        # stage configuration: [mp, me, isp] #kg, s
        self.stages = stages
        # [number of stage to jettison (nsj), mass] #kg
        # if not to be jettisoned, put 'n' in nsj parameter
        self.jettison_mass = jettison_mass
        self.g = 9.80665  # m/s²
        self.analysis_name = name

    @staticmethod
    def me_function(mp, efic):  # structural mass from eficiency
        return mp / efic - mp

    @staticmethod
    def efic_m(mp, me):  # eficiency from structural mass (me-massa estrutural)
        return mp/(mp+me)

    def deltaV(self, mpay=None):  # calculo automatic de deltaV dos estagios
        stages = self.stages
        jettison_mass = self.jettison_mass
        if mpay is None:
            mpay = self.payload_mass
        else:
            self.payload_mass = mpay

        dv = 0
        total_mass = 0
        me_i = 0
        me_f = 0
        if jettison_mass is None:
            pass
        else:
            for i in range(len(jettison_mass)):
                total_mass += jettison_mass[i][1]

        for i in range(len(stages)):  # variando estágios
            if jettison_mass is None:
                pass
            else:
                for k in range(len(jettison_mass)):  # avaliacao alijamento de massa
                    if isinstance(jettison_mass[k][0], int):
                        if i == jettison_mass[k][0] - 1:
                            me_i = total_mass
                            total_mass = total_mass - jettison_mass[k][1]
                            me_f = total_mass
                            break
                        else:
                            me_i = total_mass
                            me_f = total_mass

            mo = 0
            for j in range(i, len(stages)):  # calculando massa inicial de cada passo
                mo = mo + stages[j][0] + stages[j][1]

            mf = 0
            for j in range(i, len(stages)):  # calculando massa final de cada passo
                if j == i:
                    mf = mf + stages[j][1]
                else:
                    mf = mf + stages[j][0] + stages[j][1]

            dvl = self.g * stages[i][2] * np.log((mo + me_i + mpay) / (mf + me_f + mpay))  # dV
            dv = dv + dvl
        return dv

    def add_stage(self, isp, position='last', mp=0, me=0, efic=0):  # adicionao estagio
        if (efic == 0) and (me == 0):
            print("me, efic: a least one must be informed. Stage not inserted")
            return
        elif me == 0:
            me = self.me_function(mp, efic)
            new_stage = [mp, me, isp]
            if position == 'last':
                self.stages.append(new_stage)
            elif isinstance(position, int):
                self.stages.insert(position, new_stage)
            else:
                print("Position parameter must be 'last' or int (number of stage, down to up)")
        elif efic == 0:
            new_stage = [mp, me, isp]
            if position == 'last':
                self.stages.append(new_stage)
                return
            elif isinstance(position, int):
                self.stages.insert(position-1, new_stage)
                return
            else:
                print("Position parameter must be 'last' or int (number of stage, down to up)")
                return

    def remove_stage(self, position='last'):
        if position == 'last':
            self.stages.pop(-1)
            return
        elif isinstance(position, int):
            self.stages.pop(position-1)
            return
        else:
            print("Position parameter must be 'last' or int (number of stage, down to up)")
            return

    def update_stage(self, position, isp, mp, me=0, efic=0):
        self.remove_stage(position=position)
        self.add_stage(position=position, isp=isp, mp=mp, me=me, efic=efic)

    def update_jettison_mass(self, jettison_mass: list):
        self.jettison_mass = jettison_mass

    def update_payload_mass(self, payload_mass: float):
        self.payload_mass = payload_mass

    def show_up_stages(self, delivery='no'):
        df = pd.DataFrame(self.stages, columns=['mp', 'me', 'isp'])
        print(30*'-*')
        print("Stages Configuration - {}".format(self.analysis_name))
        print(df)
        print(30*'-*')
        if delivery.lower() == 'yes':
            return df
        elif delivery.lower() == 'no':
            pass
        else:
            print("delivery option must be 'yes' or 'no'")

    def get_minimal_prop_qtd(self, deltaV_target: float, stage_position='last', efic: float = 0,
                             mpay: float = None, max_prop: float = 10000, step=1):
        if mpay is None:
            pass
        else:
            self.payload_mass = mpay
        if stage_position.lower() == 'last':
            position = -1
        elif isinstance(stage_position, int):
            position = stage_position-1
        else:
            print("stage position must be 'last' or the number of stage.")
            return
        if efic == 0:
            efic = self.efic_m(self.stages[position][0], self.stages[position][1])

        self.stages[position][0] = 0
        self.stages[position][1] = self.me_function(self.stages[position][0], efic)
        dv = self.deltaV(mpay=self.payload_mass)
        x_mprop = [self.stages[position][0]]
        y_dv = [dv]
        y_dv_target = [deltaV_target]
        while dv < deltaV_target:
            self.stages[position][0] += step
            self.stages[position][1] = self.me_function(self.stages[position][0], efic)
            dv = self.deltaV()
            x_mprop.append(self.stages[position][0])
            y_dv.append(dv)
            y_dv_target.append(deltaV_target)
            if self.stages[position][0] > max_prop:
                break
        if dv > deltaV_target:
            print("DeltaV reached")
        else:
            print("DeltaV not reached")
        return self.stages[position][0], x_mprop, y_dv, y_dv_target

    def get_max_payload_mass(self, deltaV_reference: float, step=1):
        self.payload_mass = 0
        dv = self.deltaV()
        x_mpay = [self.payload_mass]
        y_dv = [dv]
        while dv > deltaV_reference:
            self.payload_mass += step
            dv = self.deltaV()
            x_mpay.append(self.payload_mass)
            y_dv.append(dv)
        return self.payload_mass, x_mpay, y_dv

#A single GM and GM set class
import numpy as np
from scipy.interpolate import interp1d

class GM:
    
    def __init__(self, accel, dt, unit = 9.81, label='A single GroundMotion'):
        self.label = label
        self.accel = np.append(0.0,accel)
        self.dt = dt
        self.steps = len(self.accel)
        self.time = np.arange(0,self.steps) * self.dt
        self.dur = self.time[-1]
        self.PGA = max(abs(self.accel))
        self.unit = unit
        self.scale_factor = 1.0
        self.is_spectrum = False
        self.accel_scaled = self.accel

    #NewMark version 1
    def Newmark(self, dt, AccelG, Ts, damping=0.05, beta=0.25, gamma=0.5):
        length = len(Ts)
        m = 1.0
        k = 4.0 * np.pi**2 * m / Ts**2
        cc = 2.0 * (k * m)**0.5
        c = cc * damping
        a_s = np.zeros(length)
        v_s = np.zeros(length)
        d_s = np.zeros(length)
        a_s_abs = np.zeros(length)
        meq = m + c * gamma * dt + k * beta * dt**2
        for i in range(1,len(AccelG)):
            dag = AccelG[i] - AccelG[i-1]
            a1 = a_s[-1]
            v1 = v_s[-1]
            d1 = d_s[-1]
            dd = v1 * dt + a1 * dt**2 / 2.0
            dv = a1 * dt
            da = (-1.0 * m * dag - c * dv - k * dd) / meq
            a2 = a1 + da
            v2 = v1 + dv + da * gamma * dt
            d2 = d1 + dd + da * beta * dt**2
            a_s_abs = np.vstack((a_s_abs, a2 +  AccelG[i] * np.ones(length)))
            a_s = np.vstack((a_s, a2))
            v_s = np.vstack((v_s, v2))
            d_s = np.vstack((d_s, d2))

        return v_s, d_s, a_s, a_s_abs

    def response(self, T, damping=0.05, beta=0.25, gamma=0.5):
        v_s, d_s, a_s, a_s_abs = self.Newmark(self.dt, self.accel * self.unit, np.array([T]), damping=damping, beta=beta, gamma=gamma)
        return a_s/self.unit, v_s, d_s, a_s_abs/self.unit

    def Spectrum(self, Ts, damping=0.05, beta=0.25, gamma=0.5):
        self.Ts = Ts
        v_s, d_s, a_s, a_s_abs = self.Newmark(self.dt, self.accel * self.unit, Ts, damping=damping, beta=beta, gamma=gamma)
        A = a_s/self.unit
        A_abs = a_s_abs/self.unit
        V = v_s
        D = d_s
        self.Sa = np.max(abs(A_abs),axis=0)
        self.Sd = np.max(abs(D), axis=0)
        self.Sv = np.max(abs(V), axis=0)
        self.is_spectrum = True
        #self.Sa_f = interp1d(self.Ts, self.Sa)
        #self.Sd_f = interp1d(self.Ts, self.Sd)
        #self.Sv_f = interp1d(self.Ts, self.Sv)
        self.scale(self.scale_factor)

    def Sa_f(self, T, is_scale=False):
        if is_scale:
            temp_f = interp1d(self.Ts, self.Sa_scaled)
        else:
            temp_f = interp1d(self.Ts, self.Sa)
        return temp_f(T)

    def Sv_f(self, T, is_scale=False):
        if is_scale:
            temp_f = interp1d(self.Ts, self.Sv_scaled)
        else:
            temp_f = interp1d(self.Ts, self.Sv)
        return temp_f(T)[0]

        
    def Sd_f(self, T, is_scale=False):
        if is_scale:
            temp_f = interp1d(self.Ts, self.Sd_scaled)
        else:
            temp_f = interp1d(self.Ts, self.Sd)
        return temp_f(T)

    def scale(self, scale_factor):
        self.scale_factor = scale_factor
        self.accel_scaled = self.accel * self.scale_factor
        self.PGA_scaled = self.PGA * self.scale_factor
        if self.is_spectrum:
            self.Sa_scaled = self.Sa * self.scale_factor
            self.Sv_scaled = self.Sv * self.scale_factor
            self.Sd_scaled = self.Sd * self.scale_factor

    def scale_PGA(self, PGA_target):
        factor = PGA_target / self.PGA
        self.scale(factor)

    def savetxt_scaled(self, path):
        np.savetxt(path, np.array([self.accel_scaled[1:]]).T)





#GM set class
class GM_Set():

    def __init__(self):
        self.GMs = list([]);
        self.GM_num = 0
        self.dts = np.array([])
        self.labels = list([])
        self.scale_factors = np.array([])

    def append(self, accel, dt, unit = 9.81, label='A single GroundMotion'):
        cur_GM = GM(accel, dt, unit = unit, label=label)
        self.append_GM(cur_GM)

    def append_GM(self, GM):
        self.GMs.append(GM)
        self.dts = np.append(self.dts, GM.dt)
        self.labels.append(GM.label)
        self.scale_factors = np.append(self.scale_factors, GM.scale_factor)
        self.GM_num += 1

    def Spectrum(self, Ts, damping=0.05, beta=0.25, gamma=0.5):
        self.Ts = Ts
        for i in range(self.GM_num):
            cur_GM = self.GMs[i]
            print cur_GM.label
            if not cur_GM.is_spectrum:
                cur_GM.Spectrum(Ts, damping=damping, beta=beta, gamma=gamma)
            if i == 0:
                self.Sa_list = cur_GM.Sa_scaled
                self.Sv_list = cur_GM.Sv_scaled
                self.Sd_list = cur_GM.Sd_scaled
            else:
                self.Sa_list = np.vstack((self.Sa_list, cur_GM.Sa_scaled))
                self.Sv_list = np.vstack((self.Sv_list, cur_GM.Sv_scaled))
                self.Sd_list = np.vstack((self.Sd_list, cur_GM.Sd_scaled))
        self.Sa_mean = np.mean(self.Sa_list, axis=0)
        self.Sv_mean = np.mean(self.Sv_list, axis=0)
        self.Sd_mean = np.mean(self.Sd_list, axis=0)

    def savetxt(self, path):
        output_Sa = np.vstack((self.Ts, self.Sa_mean, self.Sa_list))
        output_Sv = np.vstack((self.Ts, self.Sv_mean, self.Sv_list))
        output_Sd = np.vstack((self.Ts, self.Sd_mean, self.Sd_list))
        np.savetxt(path + '_Sa.out', output_Sa.T)
        np.savetxt(path + '_Sv.out', output_Sv.T)
        np.savetxt(path + '_Sd.out', output_Sd.T)

    def loadtxt(self,path):
        input_Sa = np.loadtxt(path + '_Sa.out').T
        input_Sv = np.loadtxt(path + '_Sv.out').T
        input_Sd = np.loadtxt(path + '_Sd.out').T
        self.Ts = input_Sa[0]
        Sa_mean = input_Sa[1]
        Sa_list = input_Sa[2:]
        Sv_mean = input_Sv[1]
        Sv_list = input_Sv[2:]
        Sd_mean = input_Sd[1]
        Sd_list = input_Sd[2:]
        for i in range(self.GM_num):
            GM = self.GMs[i]
            GM.Ts = self.Ts
            GM.Sa = Sa_list[i]
            GM.Sv = Sv_list[i]
            GM.Sd = Sd_list[i]
            GM.is_spectrum = True
            GM.scale(GM.scale_factor)
        self.Spectrum(self.Ts)

    def scale_Sa_mean(self, T1, Sa_target):
        temp_f = interp1d(self.Ts, self.Sa_mean)
        Sa_T1 = temp_f(T1)
        scale_factor = Sa_target / Sa_T1
        self.scale_factors = self.scale_factors * scale_factor
        self.re_scale()

    def scale_Sa(self, T1, Sa_target):
        for i in range(self.GM_num):
            cur_GM = self.GMs[i]
            temp_Sa = cur_GM.Sa_f(T1, is_scale=False)
            self.scale_factors[i] = Sa_target / temp_Sa
        self.re_scale()

    def scale_PGA(self, PGA_target):
        for GM in self.GMs:
            GM.scale_PGA(PGA_target)
        self.Spectrum(self.Ts)

            

    def re_scale(self):
        for scale_factor, GM in zip(self.scale_factors, self.GMs):
            GM.scale(scale_factor)
        self.Spectrum(self.Ts)

    def reset_scale(self, factor=1.0):
        self.scale_factors = np.ones(self.GM_num) * factor
        self.re_scale()

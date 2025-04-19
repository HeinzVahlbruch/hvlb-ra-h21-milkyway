import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.signal import find_peaks, peak_prominences, chirp, peak_widths

class Observation:
    c_kms = 299792458
    def __init__(self, obs_id, file_path):
        self._obs_id = obs_id
        self._file_path= file_path
        self._instrument = ""
        self._origin = ""
        self._observer = ""
        self._telescope = ""
        self._longitude = -1
        self._latitude = -1
        self._obs_date = ""
    @property
    def obs_id(self):
        return self._obs_id
    @property
    def instrument(self):
        return self._instrument
    @property
    def origin(self):
        return self._origin
    @property
    def observer(self):
        return self._observer
    @property
    def telescope(self):
        return self._telescope
    @property
    def longitude(self):
        return self._longitude
    @property
    def latitude(self):
        return self._latitude
    @property
    def obs_date(self):
        return self._obs_date
    @property
    def values(self):
        return self._dfValues
    def read_data(self):
        self._dfValues = pd.DataFrame()
    def get_peaks(self, height, min_prominence, width, min_velo, max_velo):
        peaksAll, _ = find_peaks(self._dfValues['t'], height=height, width=width, prominence=(min_prominence, None))
        peakList = []
        for x in peaksAll:
            if (self._dfValues['v'][x] > min_velo) & (self._dfValues['v'][x] < max_velo):
                peakList.append(x)
        peaks = np.array(peakList)
        return peaks
    def get_mins(self, height):
        mins, _ = find_peaks(self._dfValues['t']*-1, height=height)
        return mins
    def get_max_velocity(self, mins, peaks):
        if (self._dfValues.empty) | (peaks.size == 0) | (mins.size == 0):
            maxv = 0
        else:
            maxv = max(mins[self._dfValues.v[mins] > max(self._dfValues.v[peaks])])
        return maxv 

class SalsaFITS(Observation):
    def __init__(self, obs_id, file_path):
        super().__init__(obs_id, file_path)
        self.read_data()
    def read_data(self):
        hdu_list = fits.open(self._file_path)
        hdr = hdu_list[0].header
        
        naxis = hdr['NAXIS1']
        crval1 = hdr['CRVAL1']
        crpix = hdr['CRPIX1']
        cdelt1 = hdr['CDELT1']
        vLSR = hdr['VELO-LSR'] * 1000
        self._longitude = int(round(float(hdr['CRVAL2'])))
        self._latitude = int(round(float(hdr['CRVAL3'])))
        self._obs_date = hdr['DATE-OBS']
        self._instrument = hdr['INSTRUME']
        self._origin = hdr['ORIGIN']
        self._observer = hdr['OBSERVER']
        self._telescope = hdr['TELESCOP']

        t = hdu_list[0].data[0,0,:]
        hdu_list.close()
        
        l = np.full((naxis), self._longitude)
        f = np.array([crval1 + (i - crpix) * cdelt1 for i in range(naxis)])
        v = (- self.c_kms * (f - crval1) / crval1 - vLSR) / 1000
                
        spec = {'v': v, 'f':f, 't': t, 'l': l}
        self._dfValues = pd.DataFrame(spec)

def simple_plot(obs_dic, steps):
    for k, s in obs_dic.items():
        if k % steps == 0:
            df_spec = s.values
            plt.plot(df_spec.v, df_spec.t)

    plt.xlabel('V')  
    plt.ylabel('T')  
    plt.show()

def detail_plot(spec, peak_finder_params):
    df_spec = spec.values

    peak_height = peak_finder_params['peak_height']
    min_prominence = peak_finder_params['min_prominence']
    peak_width = peak_finder_params['peak_width']
    min_velo = peak_finder_params['min_velo']
    max_velo = peak_finder_params['max_velo']
    minimum_height = peak_finder_params['minimum_height']

    peaks = spec.get_peaks(peak_height, min_prominence, peak_width, min_velo, max_velo)
    mins = spec.get_mins(minimum_height)
    maxv_min = spec.get_max_velocity(mins, peaks)

    plt.plot(df_spec.v, df_spec.t, color="Orange")
    plt.plot(df_spec.v[peaks], df_spec.t[peaks], "x")
    plt.plot(df_spec.v[maxv_min], df_spec.t[maxv_min], "o")
    plt.xlabel('V')  
    plt.ylabel('T')  
    plt.show()

def peaks_plots_analysis(obs_dic, columns, peak_finder_params):
    long, velocities, max_velocities, instrument, telescope = [],[],[],[],[]
    pir = pic = 0
    
    peak_height = peak_finder_params[0]['peak_height']
    min_prominence = peak_finder_params[0]['min_prominence']
    peak_width = peak_finder_params[0]['peak_width']
    min_velo = peak_finder_params[0]['min_velo']
    max_velo = peak_finder_params[0]['max_velo']
    minimum_height = peak_finder_params[0]['minimum_height']
    
    obs_count = len(obs_dic.items())
    pr = int(obs_count / columns)

    fig, ax = plt.subplots(pr, columns, figsize=(10*columns,120))

    for k, s in obs_dic.items():
        df_spec = s.values

        if (k in peak_finder_params):
            tmp_peak_height = peak_finder_params[k]['peak_height']
            tmp_min_prominence = peak_finder_params[k]['min_prominence']
            tmp_peak_width = peak_finder_params[k]['peak_width']
            tmp_min_velo = peak_finder_params[k]['min_velo']
            tmp_max_velo = peak_finder_params[k]['max_velo']
            tmp_minimum_height = peak_finder_params[k]['minimum_height']

            peaks = s.get_peaks(tmp_peak_height, tmp_min_prominence, tmp_peak_width, tmp_min_velo, tmp_max_velo)
            mins = s.get_mins(tmp_minimum_height)
        else:
            peaks = s.get_peaks(peak_height, min_prominence, peak_width, min_velo, max_velo)
            mins = s.get_mins(minimum_height)
        maxv_min = s.get_max_velocity(mins, peaks)

        for p in peaks:
            long.append(k)
            velocities.append(df_spec.v[p])
            max_v = df_spec.v[maxv_min]
            max_velocities.append(max_v)
            instrument.append(s.instrument)
            telescope.append(s.telescope)
            
        title = 'LSR Velocity for Galactic Coordinates Long ' + str(s.longitude) + ' - Lat ' + str(s.latitude)

        ax[pir, pic].plot(df_spec.v, df_spec.t, c='orange')
        ax[pir, pic].plot(df_spec.v[peaks], df_spec.t[peaks], "x", c='blue')
        ax[pir, pic].plot(df_spec.v[maxv_min], df_spec.t[maxv_min], "o", c='green')
        ax[pir, pic].set_title(title)
        ax[pir, pic].set_xlabel('LSR Velocity (km/s)')
        ax[pir, pic].set_ylabel('Relative Intensity', labelpad=150)
        ax[pir, pic].spines['left'].set_position(('data', 0))
        ax[pir, pic].spines['bottom'].set_position(('data', 0))
        ax[pir, pic].spines['top'].set_visible(False)
        ax[pir, pic].spines['right'].set_visible(False)
        ax[pir, pic].grid(axis='y')
        ax[pir, pic].set_ylim(0,170)
        ax[pir, pic].set_xlim(-120,172)

        pic+=1
        if pic == columns:
            pir+=1
            pic = 0

    plt.show()

    peak_results = pd.DataFrame({'longitude':long, 'velocity':velocities, 'max_velocity':max_velocities, 'telescope':telescope, 'instrument':instrument})
    return (peak_results)

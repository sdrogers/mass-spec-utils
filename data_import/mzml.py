import numpy as np
import pymzml
import bisect

class MZMLScan(object):
    def __init__(self,scan_no,source_file,ms_level,peaks,rt_in_minutes,precursor_mz = None):
        self.scan_no = scan_no
        self.source_file = source_file
        self.ms_level = ms_level
        self.peaks = peaks
        self.rt_in_minutes = rt_in_minutes
        self.rt_in_seconds = 60.0*rt_in_minutes
        self.precursor_mz = precursor_mz
        self.previous_ms1 = None
        self.next_ms1 = None
        self.purity = None 
        self.boxes = []

        self.peaks = sorted(self.peaks,key = lambda x: x[0]) # just to be safe
    
    def compute_purity(self,ms1_scan,half_isolation_width = 0.75):
        pre = ms1_scan.get_precursor(self.precursor_mz)
        if pre:
            tot = ms1_scan.get_total_int(self.precursor_mz,half_isolation_width = half_isolation_width)
            return pre[1]/tot
        else:
            return None
    
    def get_total_int(self,precursor_mz,half_isolation_width=0.75):
        mz,intensity = zip(*self.peaks)
        mz = np.array(mz)
        intensity = np.array(intensity)
        pos1 = bisect.bisect_right(mz,precursor_mz - half_isolation_width)
        pos2 = bisect.bisect_right(mz,precursor_mz + half_isolation_width)
        return sum(intensity[pos1:pos2])

    def get_peaks_in_window(self,window_center,half_isolation_width=0.75):
        """
            A method to return all the peaks within a window

            Arguments:
            window_center -- the center of the window
            half_isolation_width -- half of the window width (i.e. all peaks
            are extracted that are with += half_isolation_width of window_center)
        """
        mz,intensity = zip(*self.peaks)
        mz = np.array(mz)
        intensity = np.array(intensity)
        pos1 = bisect.bisect_right(mz,window_center - half_isolation_width)
        pos2 = bisect.bisect_right(mz,window_center + half_isolation_width)
        peaks = []
        for pos in range(pos1,pos2):
            peaks.append((mz[pos],intensity[pos]))
        return peaks
    
    def get_precursor(self,precursor_mz,max_error = 10):
        mz,intensity = zip(*self.peaks)
        mz = np.array(mz)
        err = abs(mz - precursor_mz)
        pos = np.argmin(err)
        ppm_err = 1e6*err[pos]/precursor_mz
        if ppm_err <= max_error:
            return (mz[pos],intensity[pos])
        else:
            print("Can't find peak within {} ppm, closest = {}".format(max_error,ppm_err))
            return None

    def get_intensity(self,query_mz,dp = 2):
        # checks for exact mz match
        inte = 0.0
        b,_ = safe_round(query_mz,dp)
        for mz,i in self.peaks:
            a,_ = safe_round(mz,dp)
            if a==b:
                inte += i

        return inte

def filter_spectrum(spectrum,min_intensity = 5e3):
    new_peaks = list(filter(lambda x: x[1]>=min_intensity,spectrum.peaks))
    spectrum.peaks = new_peaks


def safe_round(mz,dp):
    mz*=10**dp
    i_mz = int(mz)
    rem = mz - i_mz
    if rem >= 0.5:
        i_mz += 1
    return i_mz,i_mz/10**dp

class MZMLFile(object):
    def __init__(self,file_name):
        self.file_name = file_name
        self.scans = []
        self._load_file()
        self._link_scans()

    def _link_scans(self):
        # populates the previous ms1 and next ms1 attributes
        # find the first ms1
        for i,scan in enumerate(self.scans):
            # look back for previous ms1
            scan.previous_ms1 = self._look_back(i)
            scan.next_ms1 = self._look_forward(i)

    def _look_back(self,pos):
        # find the first ms1 scan looking back from pos
        temp_pos = pos - 1
        while temp_pos >= 0 and self.scans[temp_pos].ms_level > 1:
            temp_pos -=1
        if temp_pos >= 0:
            return self.scans[temp_pos]
        else:
            return None

    def _look_forward(self,pos):
        temp_pos = pos + 1
        while temp_pos < len(self.scans) and self.scans[temp_pos].ms_level > 1:
            temp_pos += 1
        if temp_pos < len(self.scans):
            return self.scans[temp_pos]
        else:
            return None
    
    def _load_file(self):
        reader = pymzml.run.Reader(self.file_name,obo_version = '4.1.12')
        scan_no = 0
        for scan in reader:
            rt = scan.scan_time_in_minutes()
            ms_level = scan.ms_level
            peaks = scan.peaks('centroided')
            if ms_level == 2:
                precursor_mz = scan.selected_precursors[0]['mz']
            else:
                precursor_mz = None
            self.scans.append(MZMLScan(scan_no,self.file_name,
                              ms_level,peaks,rt,precursor_mz))
            scan_no += 1
        print("Loaded {} scans".format(scan_no))

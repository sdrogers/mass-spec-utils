# code for ms2 processing
import pymzml
import csv
import numpy as np
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

class PickedBox(object):
    def __init__(self,peak_id,mz,rt,mz_min,mz_max,rt_min,rt_max,area = None,height = None):
        self.peak_id = peak_id
        self.rt = rt
        self.rt_in_seconds = rt*60.0
        self.mz = mz
        self.rt_in_minutes = rt
        self.mz_range = [mz_min,mz_max]
        self.rt_range = [rt_min,rt_max]
        self.ms2_scans = []
        self.area = area
        self.height = height
        self.rt_range_in_seconds = [60.0*r for r in self.rt_range]
    
    def __str__(self):
        return str(self.peak_id) + ": " + str(self.mz_range) + " " + str(self.rt_range)
    
    def find_ms2_scans(self,mzml_file):
        mz_range = self.mz_range
        rt_range = self.rt_range
        scans = list(filter(lambda x: 
                    x.ms_level == 2 and
                            x.precursor_mz >= mz_range[0] and
                            x.precursor_mz <= mz_range[1] and
                            x.rt_in_minutes >= rt_range[0] and
                            x.rt_in_minutes <= rt_range[1],
                    mzml_file.scans))
        self.ms2_scans = scans
        

def load_picked_boxes(csv_name):
    boxes = []
    with open(csv_name,'r') as f:
        reader = csv.reader(f)
        heads = next(reader)
        id_pos = 0
        mz_pos = 1
        rt_pos = 2
        rt_start = [i for i,h in enumerate(heads) if 'RT start' in h][0]
        rt_end = [i for i,h in enumerate(heads) if 'RT end' in h][0]
        mz_min = [i for i,h in enumerate(heads) if 'm/z min' in h][0]
        mz_max = [i for i,h in enumerate(heads) if 'm/z max' in h][0]
        try:
            area_pos = [i for i,h in enumerate(heads) if 'Peak area' in h][0]
            height_pos = [i for i,h in enumerate(heads) if 'Peak height' in h][0]
        except:
            area_pos = None
            height_pos = None
        for line in reader:
            if area_pos:
                boxes.append(PickedBox(int(line[id_pos]),
                                    float(line[mz_pos]),
                                    float(line[rt_pos]),
                                    float(line[mz_min]),
                                    float(line[mz_max]),
                                    float(line[rt_start]),
                                    float(line[rt_end]),
                                    area = float(line[area_pos]),
                                    height = float(line[height_pos])))
            else:
                boxes.append(PickedBox(int(line[id_pos]),
                                    float(line[mz_pos]),
                                    float(line[rt_pos]),
                                    float(line[mz_min]),
                                    float(line[mz_max]),
                                    float(line[rt_start]),
                                    float(line[rt_end])))

        
    return boxes


def map_boxes_to_scans(mzml_file,boxes,half_isolation_window = 0.75,allow_last_overlap = False):
    scans2boxes = {}
    boxes2scans = {}
    for scan in mzml_file.scans:
        if scan.ms_level == 1:
            continue
        rt = scan.rt_in_minutes
        if allow_last_overlap:
            previous_ms1 = scan.previous_ms1
            if previous_ms1:
                rt = previous_ms1.rt_in_minutes
        min_mz = scan.precursor_mz - half_isolation_window
        max_mz = scan.precursor_mz + half_isolation_window
        sub_boxes = list(filter(lambda x: 
                        min_mz <= x.mz_range[1] and
                        max_mz >= x.mz_range[0] and
                        rt >= x.rt_range[0] and 
                        rt <= x.rt_range[1], boxes))
        if len(sub_boxes) > 0:
            scans2boxes[scan] = sub_boxes
            for box in sub_boxes:
                if not box in boxes2scans:
                    boxes2scans[box] = []
                boxes2scans[box].append(scan)
    return scans2boxes,boxes2scans

def traverse_boxes_scans(boxes2scans,scans2boxes,start_scan):
    subset_scans = set()
    subset_boxes = set()
    scans_to_explore = set()
    scans_to_explore.add(start_scan)
    while len(scans_to_explore) > 0:
        this_scan = scans_to_explore.pop()
        subset_scans.add(this_scan)
        sub_boxes = scans2boxes[this_scan]
        for b in sub_boxes:
            if not b in subset_boxes:
                subset_boxes.add(b)
                for sc in boxes2scans[b]:
                    if not sc in subset_scans:
                        scans_to_explore.add(sc)
    return subset_scans,subset_boxes
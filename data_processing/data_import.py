# code for ms2 processing)
import pymzml
import csv
import numpy as np
import bisect
from scipy.optimize import nnls
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

def interpolate(precursor_mz,ms1_1,ms1_2,forced_rt):
    pre = ms1_1.get_precursor(precursor_mz)
    post = ms1_2.get_precursor(precursor_mz)
    if pre == None or post == None:
        return 0.0
    rt1 = ms1_1.rt_in_minutes
    rt2 = ms1_2.rt_in_minutes
    rt_fac = (forced_rt - rt1)/(rt2-rt1)
    inte_diff = post[1] - pre[1]
    return pre[1]+inte_diff*rt_fac

def create_deconvolution_matrices(subset_scans,subset_boxes,min_ms2_intensity = 5000):
    n_scans = len(subset_scans)
    n_ms1 = len(subset_boxes)
    scan_list = list(subset_scans)
    box_list = list(subset_boxes)
    ms1_intensity = np.zeros((n_scans,n_ms1))
    for i,s in enumerate(scan_list):
        for j,b in enumerate(box_list):
            inte = interpolate(b.mz,s.previous_ms1,s.next_ms1,s.rt_in_minutes)
            ms1_intensity[i,j] = inte
    

    frags_idx = {}
    frag_pos = 0
    all_specs = []
    for i,s in enumerate(scan_list):
        rspec = []
        for mz,intensity in s.peaks:
            rmz = round(mz,2)
            if intensity < min_ms2_intensity:
                continue
            if not rmz in frags_idx:
                frags_idx[rmz] = frag_pos
                frag_pos += 1
            rspec.append((rmz,intensity))
        all_specs.append(rspec)
        
    frag_mat = np.zeros((n_scans,len(frags_idx)))
    for i,spec in enumerate(all_specs):
        for rmz,intensity in spec:
            pos = frags_idx[rmz]
            frag_mat[i,pos] += intensity


    box_idx = {box:i for i,box in enumerate(box_list)}
    return ms1_intensity,frag_mat,frags_idx,box_idx

def deconvolve(ms1_intensity,frag_mat,frags_idx):
    from scipy.stats import pearsonr
    n_scans,n_boxes = ms1_intensity.shape
    if n_scans < n_boxes:
        print("NO. CAN. DO")
        return None
    
    # check the correlation
    corr_thresh = 0.75
    n_corr_greater_thresh = 0
    for i in range(n_scans-1):
        for j in range(i+1,n_scans):
            rho,h = pearsonr(ms1_intensity[i,:],ms1_intensity[j,:])
            if rho > corr_thresh:
                n_corr_greater_thresh += 1
    if n_corr_greater_thresh > n_scans - n_boxes:

        print("Looks nasty, scared")
        print("{} pairs of observations (of {}) have correlation > {}".format(n_corr_greater_thresh,n_scans*(n_scans-1)/2,corr_thresh))
        print(ms1_intensity)
        # return
        

    B = frag_mat
    A = ms1_intensity
    all_o = []
    out_mat = []
    for i in range(len(frags_idx)):
        o = nnls(np.matrix(A),B[:,i],maxiter=10000)
        all_o.append(o)
        out_mat.append(o[0])
    out_mat = np.array(out_mat)
    out_mat = out_mat.T
    # out_mat/=out_mat.sum(axis=1)[:,None]
    return out_mat

def interpolate_ms1_intensities(ms2scan,ms1scan_1,ms1scan_2,half_isolation_width = 0.75,matching_ppm = 10,min_ms1_intensity = 0):
    """
    Takes a particular ms2 scan and:
        - finds all the peaks within the isolation window in both ms1 scans
        - matches them across both scans (setting intensity = 0 if something
        is not found)
        - Interpolates them to get their intensity at the rt of the ms2 scan
        - Returns a list of mz,intensity tuples
    Arguments:
    ms2scan -- the ms2scan of interest
    ms1scan_1 -- the first (in rt) ms1scan
    ms1scan_2 -- the second (in rt) ms1scan

    Keyword arguments:
    half_isolation_width -- half of the isolation window width (default = 0.75)
    """
    
    # get the ms1 peaks of interest
    peaks1 = ms1scan_1.get_peaks_in_window(ms2scan.precursor_mz,half_isolation_width = half_isolation_width)
    peaks2 = ms1scan_2.get_peaks_in_window(ms2scan.precursor_mz,half_isolation_width = half_isolation_width)
    # match them
    matched_peaks = match_peaks(peaks1,peaks2,ppm=10)


    ms1scan_1.rt_in_minutes
    ms1scan_2.rt_in_minutes

    # interpolate
    interpolated_intensities = []
    rt_fac = (ms2scan.rt_in_minutes - ms1scan_1.rt_in_minutes)/(ms1scan_2.rt_in_minutes - ms1scan_1.rt_in_minutes)
    for mz,i1,i2 in matched_peaks:
        i_diff = i2-i1
        interpolated_i = i1 + rt_fac*i_diff
        interpolated_intensities.append((mz,interpolated_i))
    
    return list(filter(lambda x: x[1] >= min_ms1_intensity,interpolated_intensities))

def match_peaks(peaks1,peaks2,ppm=10):
    """
    Matches two lists of peaks to find the overlaps (within tolerance)
    """
    # sort by descending intensity
    peaks1.sort(key = lambda x: x[1],reverse = True)
    peaks2.sort(key = lambda x: x[1],reverse = True)
    used_from_2 = set()
    remaining_in_2 = list(range(len(peaks2)))
    matches = []
    for i,(mz,intensity) in enumerate(peaks1):
        min_err = 1e6
        min_pos = -1
        for pos,j in enumerate(remaining_in_2):
            err = abs(1e6*(mz - peaks2[j][0])/mz)
            if err < min_err:
                min_err = err
                min_pos = pos
        if min_err <= ppm:
            matches.append((i,remaining_in_2[min_pos]))
            del remaining_in_2[min_pos]
            found_match = True
        if not found_match:
            matches.append((i,None))
    for j in remaining_in_2:
        matches.append((None,j))
    
    # degubbing
    # for p1,p2 in matches:
    #     line = ""
    #     if not p1 is None:
    #         line += "{:.4f}".format(peaks1[p1][0])
    #     else:
    #         line += "None"
    #     if not p2 is None:
    #         line += "\t{:.4f}".format(peaks2[p2][0])
    #     else:
    #         line += "\tNone"
    #     print(line)

    # print(matches)

    output_matches = []
    for p1,p2 in matches:
        if not p1 is None and p2 is not None:
            mean_mz = (peaks1[p1][0] + peaks2[p2][0])/2
            output_matches.append((mean_mz,peaks1[p1][1],peaks2[p2][1]))
        elif p2 is None:
            output_matches.append((peaks1[p1][0],peaks1[p1][1],0.0))
        elif p1 is None:
            output_matches.append((peaks2[p2][0],0.0,peaks2[p2][1]))
    return output_matches

def simple_align(peak_set,peaks,sample_id,ppm = 10):
    import copy
    for this_mz,intensity in peaks:
        best_err = 1e6
        best_mz = None
        for mz in peak_set:
            err = 1e6*abs(mz - this_mz)/mz
            if err < best_err:
                best_err = err
                best_mz = mz
        if best_err <= ppm:
            peak_set[best_mz].append((sample_id,this_mz,intensity))
            new_mz = sum([m[1] for m in peak_set[best_mz]])/len(peak_set[best_mz])
            if not new_mz == best_mz:
                peak_set[new_mz] = peak_set[best_mz]
                del peak_set[best_mz]
        else:
            peak_set[this_mz] = [(sample_id,this_mz,intensity)]
    


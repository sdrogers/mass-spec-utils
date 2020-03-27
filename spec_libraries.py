# mb_parser
import os
import glob
import math
import random
from massbank import load_massbank
from hmdb import load_hmdb_msms_records


def modified_cosine_similarity(spectrum1,spectrum2,tol,min_match):
    if spectrum1.n_peaks == 0 or spectrum2.n_peaks == 0:
        return 0.0,[]

    spec1 = spectrum1.normalised_peaks
    spec2 = spectrum2.normalised_peaks

    zero_pairs = find_pairs(spec1,spec2,tol,shift=0.0)

    shift = spectrum1.parent_mz - spectrum2.parent_mz

    nonzero_pairs = find_pairs(spec1,spec2,tol,shift = shift)

    matching_pairs = zero_pairs + nonzero_pairs

    matching_pairs = sorted(matching_pairs,key = lambda x: x[2], reverse = True)

    used1 = set()
    used2 = set()
    score = 0.0
    used_matches = []
    for m in matching_pairs:
        if not m[0] in used1 and not m[1] in used2:
            score += m[2]
            used1.add(m[0])
            used2.add(m[1])
            used_matches.append(m)
    if len(used_matches) < min_match:
        score = 0.0
    return score,used_matches


def find_pairs(spec1,spec2,tol,shift=0):
    matching_pairs = []
    spec2lowpos = 0
    spec2length = len(spec2)
    
    for idx,(mz,intensity) in enumerate(spec1):
        # do we need to increase the lower idx?
        while spec2lowpos < spec2length and spec2[spec2lowpos][0] + shift < mz - tol:
            spec2lowpos += 1
        if spec2lowpos == spec2length:
            break
        spec2pos = spec2lowpos
        while(spec2pos < spec2length and spec2[spec2pos][0] + shift < mz + tol):
            matching_pairs.append((idx,spec2pos,intensity*spec2[spec2pos][1]))
            spec2pos += 1
        
    return matching_pairs    

def cosine_similarity(spectrum1,spectrum2,tol,min_match):
    # spec 1 and spec 2 have to be sorted by mz
    if spectrum1.n_peaks == 0 or spectrum2.n_peaks == 0:
        return 0.0,[]
    # find all the matching pairs
    
    spec1 = spectrum1.normalised_peaks
    spec2 = spectrum2.normalised_peaks
    
    matching_pairs = find_pairs(spec1,spec2,tol,shift = 0.0)
    
        
        
    matching_pairs = sorted(matching_pairs,key = lambda x:x[2],reverse = True)
    used1 = set()
    used2 = set()
    score = 0.0
    used_matches = []
    for m in matching_pairs:
        if not m[0] in used1 and not m[1] in used2:
            score += m[2]
            used1.add(m[0])
            used2.add(m[1])
            used_matches.append(m)
    if len(used_matches) < min_match:
        score = 0.0
    return score,used_matches












class SpectralLibrary(object):
    def __init__(self):
        self.records = {}
        self.sorted_record_list = []

    def _candidates(self,query_mz,ms1_tol):
        from sortedcontainers import SortedList
        pmz_list = SortedList([m.precursor_mz for m in self.sorted_record_list])
        lower = query_mz - ms1_tol
        upper = query_mz + ms1_tol
        start = pmz_list.bisect(lower)
        end = pmz_list.bisect(upper)
        return self.sorted_record_list[start:end]
    
    def _dic2list(self):
        # - makes a list of the values in records
        # - removes any without a precursor mz
        # - sorts
        record_list = list(self.records.values())
        filtered_list = list(filter(lambda x: not x.precursor_mz is None,record_list))
        return sorted(filtered_list)

    def spectral_match(self,query,
            scoring_function = 'cosine',
            ms2_tol = 0.2,
            min_match_peaks = 1,
            ms1_tol = 0.2,
            score_thresh = 0.7):
        
        candidates = self._candidates(query.precursor_mz,ms1_tol)
        hits = []
        for c in candidates:
            if scoring_function == 'cosine':
                sc,_ = cosine_similarity(query,c,ms2_tol,min_match_peaks)
            elif scoring_function == 'modified cosine':
                sc,_ = modified_cosine_similarity(query,c,ms2_tol,min_match_peaks)
            else:
                print("Unrecognised scoring function: ",scoring_function)
                return None
            if sc >= score_thresh:
                hits.append((c.spectrum_id,sc,c))
        return hits


class MassBankLibrary(SpectralLibrary):
    def __init__(self,mb_dir = '/Users/simon/git/MassBank-data/'):
        super().__init__()
        self.records = load_massbank(mb_dir = mb_dir)
        self.sorted_record_list = self._dic2list() # sorted by precursor mz

class HMDBLibrary(SpectralLibrary):
    def __init__(self,csv_file = '/Users/simon/hmdb/hmdb_metabolites.csv',
                 msms_file_folder = '/Users/simon/hmdb/hmdb_experimental_msms_spectra',
                 mode = 'positive',
                 transformations = ['[M+H]+']):
        super().__init__()
        self.records = load_hmdb_msms_records(msms_file_folder,csv_file,target_mode = 'positive',transformations = ['[M+H]+'])
        self.sorted_record_list = self._dic2list() # sorted by precursor mz



if __name__ == '__main__':

    # example of massbank library
    mbl = MassBankLibrary()

    record_list = list(mbl.records.values())
    a = random.randrange(len(record_list))
    example_spectrum = record_list[a]

    hits = mbl.spectral_match(example_spectrum,score_thresh = 0.5)
    print(example_spectrum)
    for hit in hits:
        print(hit)

    # example of HMDB
    hmdbl = HMDBLibrary()
    
    hits = hmdbl.spectral_match(example_spectrum,score_thresh = 0.5)
    for hit in hits:
        print(hit)

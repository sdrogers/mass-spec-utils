# mb_parser
import os
import glob
import math
import random


def fast_cosine_shift(spectrum1,spectrum2,tol,min_match):
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

def fast_cosine(spectrum1,spectrum2,tol,min_match):
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





class MBRecord(object):
    def __init__(self,precursor_mz,peaks,metadata,original_file,spectrum_id):
        self.peaks = peaks
        self.metadata = metadata
        self.original_file = original_file
        self.precursor_mz = precursor_mz
        self._sqrt_normalise()
        self.n_peaks = len(self.peaks)
        self.spectrum_id = spectrum_id
        

    def __str__(self):
        return ", ".join([self.spectrum_id,self.metadata['names'][0],str(self.precursor_mz),self.original_file])

    def __repr__(self):
        return self.__str__()

    def __lt__(self,other):
        return(self.precursor_mz < other.precursor_mz)

    def _sqrt_normalise(self):
        temp = []
        total = 0.0
        for mz,intensity in self.peaks:
            temp.append((mz,math.sqrt(intensity)))
            total += intensity
        norm_facc = math.sqrt(total)
        self.normalised_peaks = []
        for mz,intensity in temp:
            self.normalised_peaks.append((mz,intensity/norm_facc))


def parse_mb_file(file_name):
    with open(file_name,'r') as f:
        in_peaks = False
        peaks = []
        metadata = {}
        names = []
        precursor_mz = None
        precursor_type = None
        smiles = None
        iupac = None
        for line in f:
            if len(line) == 0:
                continue
            line = line.rstrip()
            if line.startswith('PK$PEAK'):
                in_peaks = True
                continue
            if line.startswith('MS$FOCUSED_ION: PRECURSOR_M/Z'):
                precursor_mz_string = line.split()[-1]
                try:
                    precursor_mz = float(precursor_mz_string)
                except:
                    if '/' in precursor_mz_string:
                        precursor_mz = float(precursor_mz_string.split('/')[0])
                    else:
                        precursor_mz = None
                continue
            if line.startswith('MS$FOCUSED_ION: PRECURSOR_TYPE'):
                precursor_type_string = line.split()[-1]
                if '/' in precursor_type_string:
                    precursor_type = precursor_type_string.split('/')[0]
                else:
                    precursor_type = precursor_type_string

                continue
            if line.startswith('CH$NAME'):
                names.append(' '.join(line.split(':')[1:]))
                continue
            if line.startswith('CH$SMILES'):
                smiles = line.split()[-1]
                continue
            if line.startswith('CH$IUPAC'):
                iupac = line.split()[-1]
                continue
            if in_peaks:
                if line.startswith('//'):
                    in_peaks = False
                else:
                    mz,inte,rel_inte = line.split()
                    peaks.append((float(mz),float(inte)))
        metadata['precursor_mz'] = precursor_mz
        metadata['precursor_type'] = precursor_type
        metadata['smiles'] = smiles
        metadata['iupac'] = iupac
        metadata['names'] = names
        spectrum_id = file_name.split(os.sep)[-1].split('.')[0]
        record = MBRecord(precursor_mz,peaks,metadata,file_name,spectrum_id)

        return record
        

def load_folder(folder_name,records = {},verbose = True):
    file_list = glob.glob(os.path.join(folder_name,'*.txt'))
    for file_name in file_list:
        spec_id = file_name.split(os.sep)[-1].split('.')[0]
        if verbose:
            print("Loading ",spec_id)
        assert not spec_id in records
        records[spec_id] = parse_mb_file(file_name)
        records[spec_id].metadata['spec_id'] = spec_id


def load_massbank(mb_dir = '/Users/simon/git/MassBank-data/'):
    folders = glob.glob(mb_dir + '*/')
    records = {}
    for folder in folders[:5]:
        n_records = len(records)
        print("Loading records from ",folder)
        load_folder(folder,records = records,verbose = False)
        print("\t Loaded {} new records".format(len(records)-n_records))
    return records


def dic2list(records):
    record_list = list(records.values())
    filtered_list = list(filter(lambda x: not x.precursor_mz is None,record_list))
    return sorted(filtered_list)


class MassBankLibrary(object):
    def __init__(self,mb_dir = '/Users/simon/git/MassBank-data/'):
        self.records = load_massbank(mb_dir = mb_dir)
        self.sorted_record_list = dic2list(self.records) # sorted by precursor mz

    def _candidates(self,query_mz,ms1_tol):
        from sortedcontainers import SortedList
        pmz_list = SortedList([m.precursor_mz for m in self.sorted_record_list])
        lower = query_mz - ms1_tol
        upper = query_mz + ms1_tol
        start = pmz_list.bisect(lower)
        end = pmz_list.bisect(upper)
        return self.sorted_record_list[start:end]
    
    def spectral_match(self,query,
            scoring_function = fast_cosine,
            ms2_tol = 0.2,
            min_match_peaks = 1,
            ms1_tol = 0.2,
            score_thresh = 0.7):
        
        candidates = self._candidates(query.precursor_mz,ms1_tol)
        hits = []
        for c in candidates:
            sc,_ = scoring_function(query,c,ms2_tol,min_match_peaks)
            if sc >= score_thresh:
                hits.append((c.spectrum_id,sc,c))
        return hits


if __name__ == '__main__':
    mb_dir = '/Users/simon/git/MassBank-data/Chubu_Univ'
    # files = glob.glob(os.path.join(mb_dir,'*.txt'))
    
    # record = parse_mb_file(files[0])
    # print(record)
    # records = {}
    # load_folder(mb_dir,records,verbose = True)
    # # print(len(records))

    # records = load_massbank()
    # sorted_records = dic2list(records)
    # pmz = [r.precursor_mz for r in sorted_records]
    # print(pmz)

    mbl = MassBankLibrary()

    # ca = mbl._candidates(123.345,10)
    # for c in ca:
    #     print(c)

    record_list = list(mbl.records.values())
    a = random.randrange(len(record_list))
    example_spectrum = record_list[a]

    hits = mbl.spectral_match(example_spectrum)
    print(example_spectrum)
    for hit in hits:
        print(hit)
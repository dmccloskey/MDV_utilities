from molmass.molmass import Formula

import numpy
import re

class mass_isotopomer_distributions():
    def __init__(self):
        return;

    
    def build_precursorSpectrumFromMRMs(self,peakSpectrum_I,blankSpectrum_I):
        '''extract maximum intensity peak'''
        # Input:
        #   peakSpectrum_I = {fragment:{(precursor_mass,product_mass):intensity}} 
        #   peakSpectrum_I = {fragment:{(precursor_mass,product_mass):intensity}} 
        # Output: 
        #   peakSpectrum_theoretical = {fragment:{mass:intensity}}
        #   peakSpectrum_measured = {fragment:{mass:[measuredMass,intensity]}} 
        #   peakSpectrum_corrected = {fragment:{mass:[measuredMass,intensity]}} 
        #   peakSpectrum_normalized = {fragment:{mass:[measuredMass,intensity]}} 

        fragments_I = list(peakSpectrum_I.keys());

        # round all precursor/product masses in input for comparison:
        peakSpectrum_copy_I = {};
        for frag,spec in peakSpectrum_I.items():
            peakSpectrum_tmp = {};
            for masses,intensity in spec.items():
                peakSpectrum_tmp[(numpy.around(masses[0]),numpy.around(masses[1]))] = intensity;
            peakSpectrum_copy_I[frag] = peakSpectrum_tmp;
        blankSpectrum_copy_I = {};
        for frag,spec in blankSpectrum_I.items():
            blankSpectrum_tmp = {};
            for masses,intensity in spec.items():
                blankSpectrum_tmp[(numpy.around(masses[0]),numpy.around(masses[1]))] = intensity;
            blankSpectrum_copy_I[frag] = blankSpectrum_tmp;


        peakSpectrum_theoretical = self.report_fragmentSpectrum_normMax(fragments_I,True);
        # determine masses from fragments
        masses = [];
        peakSpectrum_measured = {};
        peakSpectrum_normalized = {};
        peakSpectrum_corrected = {};
        for frag,spec in peakSpectrum_theoretical.items():
            peakSpectrum_measured[frag] = None;
            peakSpectrum_corrected[frag] = None;
            peakSpectrum_normalized[frag] = None;

            if not spec: continue; #check if a carbon is even contained in the fragment

            masses = list(spec.keys());
            masses.sort(); # sort mass in massList
            masses_rounded = numpy.around(masses); # round masses to nearest digit for comparison


            # 1. copy data from peakSpectrum_I to peakSpectrum_measured based on theoretical fragments
            # 2. generate corrected spectrum
            intensityList = [];
            if frag in peakSpectrum_I:
                precursor_masses = [k[0] for k in peakSpectrum_copy_I[frag].keys()];
                measured_spec = {};
                corrected_spec = {};
                for i,mass in enumerate(masses_rounded): #iterate through theoretical precursor masses
                    measured = 0.0;
                    corrected = 0.0;
                    if mass in precursor_masses:
                        product_masses = [k[1] for k in peakSpectrum_copy_I[frag].keys() if k[0]==mass];
                        for product in product_masses: #iterate through measured product masses
                            if frag in blankSpectrum_copy_I:
                                blank_precursor_masses = [k[0] for k in blankSpectrum_copy_I[frag].keys()];
                                if mass in blank_precursor_masses:
                                    blank_product_masses = [k[1] for k in blankSpectrum_copy_I[frag].keys() if k[0]==mass];
                                    if product in blank_product_masses:
                                        if blankSpectrum_copy_I[frag][(mass,product)]<0.5*peakSpectrum_copy_I[frag][(mass,product)]:
                                            corrected += peakSpectrum_copy_I[frag][(mass,product)]-blankSpectrum_copy_I[frag][(mass,product)];
                                            measured += peakSpectrum_copy_I[frag][(mass,product)]
                                        else:
                                            corrected += 0.0;
                                            measured += peakSpectrum_copy_I[frag][(mass,product)]
                                    else:
                                        corrected += peakSpectrum_copy_I[frag][(mass,product)];
                                        measured += peakSpectrum_copy_I[frag][(mass,product)]
                                else:
                                    corrected += peakSpectrum_copy_I[frag][(mass,product)];
                                    measured += peakSpectrum_copy_I[frag][(mass,product)]
                            else:
                                corrected += peakSpectrum_copy_I[frag][(mass,product)];
                                measured += peakSpectrum_copy_I[frag][(mass,product)];
                    measured_spec[masses[i]] = measured;
                    corrected_spec[masses[i]] = corrected;
                    intensityList.append(corrected);
                peakSpectrum_measured[frag] = measured_spec;
                peakSpectrum_corrected[frag] = corrected_spec;

            # normalize each spectrum:
            #NOTE: normalization by max to allow for later conversion to normalization by sum
            normalized = {};
            intensityListMax = max(intensityList);
            for k,v in peakSpectrum_corrected[frag].items():
                if intensityListMax != 0: normalized[k] = v/intensityListMax;
                else: normalized[k] = None;
            peakSpectrum_normalized[frag] = normalized;

        return peakSpectrum_measured, peakSpectrum_corrected, peakSpectrum_normalized;
    def build_productSpectrumFromMRMs(self,peakSpectrum_I,blankSpectrum_I):
        '''extract maximum intensity peak'''
        # Input:
        #   peakSpectrum_I = {fragment:{(product_mass,product_mass):intensity}} 
        #   peakSpectrum_I = {fragment:{(product_mass,product_mass):intensity}} 
        # Output: 
        #   peakSpectrum_theoretical = {fragment:{mass:intensity}}
        #   peakSpectrum_measured = {fragment:{mass:intensity}} 
        #   peakSpectrum_corrected = {fragment:{mass:intensity}} 
        #   peakSpectrum_normalized = {fragment:{mass:intensity}} 

        fragments_I = list(peakSpectrum_I.keys());

        # round all precursor/product masses in input for comparison:
        peakSpectrum_copy_I = {};
        for frag,spec in peakSpectrum_I.items():
            peakSpectrum_tmp = {};
            for masses,intensity in spec.items():
                peakSpectrum_tmp[(numpy.around(masses[0]),numpy.around(masses[1]))] = intensity;
            peakSpectrum_copy_I[frag] = peakSpectrum_tmp;
        blankSpectrum_copy_I = {};
        for frag,spec in blankSpectrum_I.items():
            blankSpectrum_tmp = {};
            for masses,intensity in spec.items():
                blankSpectrum_tmp[(numpy.around(masses[0]),numpy.around(masses[1]))] = intensity;
            blankSpectrum_copy_I[frag] = blankSpectrum_tmp;


        peakSpectrum_theoretical = self.report_fragmentSpectrum_normMax(fragments_I,True);
        # determine masses from fragments
        masses = [];
        peakSpectrum_measured = {};
        peakSpectrum_normalized = {};
        peakSpectrum_corrected = {};
        for frag,spec in peakSpectrum_theoretical.items():
            peakSpectrum_measured[frag] = None;
            peakSpectrum_corrected[frag] = None;
            peakSpectrum_normalized[frag] = None;

            if not spec: continue; #check if a carbon is even contained in the fragment

            masses = list(spec.keys());
            masses.sort(); # sort mass in massList
            masses_rounded = numpy.around(masses); # round masses to nearest digit for comparison


            # 1. copy data from peakSpectrum_I to peakSpectrum_measured based on theoretical fragments
            # 2. generate corrected spectrum
            intensityList = [];
            if frag in peakSpectrum_I:
                product_masses = [k[1] for k in peakSpectrum_copy_I[frag].keys()];
                measured_spec = {};
                corrected_spec = {};
                for i,mass in enumerate(masses_rounded): #iterate through theoretical product masses
                    measured = 0.0;
                    corrected = 0.0;
                    if mass in product_masses:
                        precursor_masses = [k[0] for k in peakSpectrum_copy_I[frag].keys() if k[1]==mass];
                        for precursor in precursor_masses: #iterate through measured precursor masses
                            if frag in blankSpectrum_copy_I:
                                blank_product_masses = [k[1] for k in blankSpectrum_copy_I[frag].keys()];
                                if mass in blank_product_masses:
                                    blank_precursor_masses = [k[0] for k in blankSpectrum_copy_I[frag].keys() if k[1]==mass];
                                    if precursor in blank_precursor_masses:
                                        if blankSpectrum_copy_I[frag][(precursor,mass)]<0.5*peakSpectrum_copy_I[frag][(precursor,mass)]:
                                            corrected += peakSpectrum_copy_I[frag][(precursor,mass)]-blankSpectrum_copy_I[frag][(precursor,mass)];
                                            measured += peakSpectrum_copy_I[frag][(precursor,mass)]
                                        else:
                                            corrected += 0.0;
                                            measured += peakSpectrum_copy_I[frag][(precursor,mass)]
                                    else:
                                        corrected += peakSpectrum_copy_I[frag][(precursor,mass)];
                                        measured += peakSpectrum_copy_I[frag][(precursor,mass)]
                                else:
                                    corrected += peakSpectrum_copy_I[frag][(precursor,mass)];
                                    measured += peakSpectrum_copy_I[frag][(precursor,mass)]
                            else:
                                corrected += peakSpectrum_copy_I[frag][(precursor,mass)];
                                measured += peakSpectrum_copy_I[frag][(precursor,mass)];
                    measured_spec[masses[i]] = measured;
                    corrected_spec[masses[i]] = corrected;
                    intensityList.append(corrected);
                peakSpectrum_measured[frag] = measured_spec;
                peakSpectrum_corrected[frag] = corrected_spec;

            # normalize each spectrum:
            #NOTE: normalization by max to allow for later conversion to normalization by sum
            normalized = {};
            intensityListMax = max(intensityList);
            for k,v in peakSpectrum_corrected[frag].items():
                if intensityListMax != 0: normalized[k] = v/intensityListMax;
                else: normalized[k] = None;
            peakSpectrum_normalized[frag] = normalized;

        return peakSpectrum_measured, peakSpectrum_corrected, peakSpectrum_normalized;
    def compare_peakSpectrum_normMax(self,peakSpectrum_normalized_list_I,return_theoretical = False):
        # Input:
        #   peakSpectrum_normalized_list_I = [{fragment:{mass:intensity}}]
        # Output:
        #   peakSpectrum_stats_O = {fragment:{mass:{'n':integer,
        #                                   'mean':fraction,
        #                                   'stdDev':fraction,
        #                                   'absDev':fraction}}

        fragments_all = [];
        for row in peakSpectrum_normalized_list_I:
            fragments_all.extend(list(row.keys()));
        fragments_I = list(set(fragments_all));
        #fragments_I = peakSpectrum_normalized_list_I[0].keys(); 
        peakSpectrum_theoretical = self.report_fragmentSpectrum_normMax(fragments_I,True);

        peakSpectrum_stats_O = {};
        for frag in fragments_I:
            peakSpectrum_stats_O[frag] = {'n':None,
                               'mean':None,
                               'stdDev':None,
                               'absDev':None};

            if not peakSpectrum_theoretical[frag]: continue; # no carbons in fragment

            intensityList = [];
            masses = [];
            stats = {};
            for peakSpectrum in peakSpectrum_normalized_list_I:
                intensityDict = {};
                peakSpectrumMasses = list(peakSpectrum_theoretical[frag].keys());
                for mass in peakSpectrumMasses:
                    if frag in peakSpectrum and mass in peakSpectrum[frag] and peakSpectrum[frag][mass] and peakSpectrum[frag][mass] > 0.0: 
                        intensityDict[mass] = peakSpectrum[frag][mass];
                    else: 
                        intensityDict[mass] = 0.0;
                    if not mass in masses: masses.append(mass);
                intensityList.append(intensityDict);
                ## uncomment to only compare measured masses
                #intensityDict = {};
                #peakSpectrumMasses = peakSpectrum[frag].keys();
                #for mass in peakSpectrumMasses:
                #    if peakSpectrum[frag][mass] > 0.0: 
                #        intensityDict[mass] = peakSpectrum[frag][mass];
                #        if not mass in masses: masses.append(mass);
                #intensityList.append(intensityDict);
            for mass in masses:
                stats[mass] = None;
                data = [];
                for intensity in intensityList:
                    if intensity[mass]>0.0:data.append(intensity[mass]);
                if data:
                    intensity_array = numpy.array(data);
                    if peakSpectrum_theoretical[frag][mass]:abs_dev = abs(intensity_array.mean() - peakSpectrum_theoretical[frag][mass]);
                    else: abs_dev = None;
                    stats[mass] = {'n':len(intensity_array),
                                   'mean':intensity_array.mean(),
                                   'stdDev':intensity_array.std(),
                                   'absDev':abs_dev};
                else:
                    stats[mass] = {'n':0.0,
                                   'mean':0.0,
                                   'stdDev':0.0,
                                   'absDev':None};
            if stats: peakSpectrum_stats_O[frag] = stats;

        if return_theoretical:
            return peakSpectrum_stats_O,peakSpectrum_theoretical;
        else:
            return peakSpectrum_stats_O;
    def compare_peakSpectrum_normSum(self,peakSpectrum_normalized_list_I,return_theoretical = False):
        # Input:
        #   peakSpectrum_normalized_list_I = [{fragment:{mass:[measuredMass,intensity]}}]
        # Output:
        #   peakSpectrum_stats_O = {fragment:{mass:{'n':integer,
        #                                   'mean':fraction,
        #                                   'stdDev':fraction,
        #                                   'absDev':fraction}}

        fragments_all = [];
        for row in peakSpectrum_normalized_list_I:
            fragments_all.extend(list(row.keys()));
        fragments_I = list(set(fragments_all));
        #fragments_I = peakSpectrum_normalized_list_I[0].keys(); 
        peakSpectrum_theoretical = self.report_fragmentSpectrum_normSum(fragments_I,True);

        peakSpectrum_stats_O = {};
        for frag in fragments_I:
            peakSpectrum_stats_O[frag] = {'n':None,
                               'mean':None,
                               'stdDev':None,
                               'absDev':None};

            if not peakSpectrum_theoretical[frag]: continue; # no carbons in fragment

            intensityList = [];
            masses = [];
            stats = {};
            for peakSpectrum in peakSpectrum_normalized_list_I:
                intensityDict = {};
                peakSpectrumMasses = list(peakSpectrum_theoretical[frag].keys());
                for mass in peakSpectrumMasses:
                    if frag in peakSpectrum and frag in peakSpectrum and mass in peakSpectrum[frag] and peakSpectrum[frag][mass] and peakSpectrum[frag][mass] > 0.0: 
                        intensityDict[mass] = peakSpectrum[frag][mass];
                    else: 
                        intensityDict[mass] = 0.0;
                    if not mass in masses: masses.append(mass);
                intensityList.append(intensityDict);
                ## uncomment to only compare measured masses
                #intensityDict = {};
                #peakSpectrumMasses = peakSpectrum[frag].keys();
                #for mass in peakSpectrumMasses:
                #    if peakSpectrum[frag][mass] > 0.0: 
                #        intensityDict[mass] = peakSpectrum[frag][mass];
                #        if not mass in masses: masses.append(mass);
                #intensityList.append(intensityDict);
            for mass in masses:
                stats[mass] = None;
                data = [];
                for intensity in intensityList:
                    if intensity[mass]>0.0:data.append(intensity[mass]);
                if data:
                    intensity_array = numpy.array(data);
                    if peakSpectrum_theoretical[frag][mass]:abs_dev = abs(intensity_array.mean() - peakSpectrum_theoretical[frag][mass]);
                    else: abs_dev = None;
                    stats[mass] = {'n':len(intensity_array),
                                   'mean':intensity_array.mean(),
                                   'stdDev':intensity_array.std(),
                                   'absDev':abs_dev};
                else:
                    stats[mass] = {'n':0.0,
                                   'mean':0.0,
                                   'stdDev':0.0,
                                   'absDev':None};
            if stats: peakSpectrum_stats_O[frag] = stats;
            
        if return_theoretical:
            return peakSpectrum_stats_O,peakSpectrum_theoretical;
        else:
            return peakSpectrum_stats_O;
    def report_fragmentSpectrum_normMax(self,fragments_I,round_mass=False):
        '''calculate the format spectrum as a list'''
        # Input: formula_str_I
        # Output: spectrum_lst_O

        fragmentSpectrum_tmp = {};
        fragmentSpectrum_O = {};

        for formula_str_I in fragments_I:
            fragmentSpectrum_tmp[formula_str_I] = None;
            fragmentSpectrum_O[formula_str_I] = None;
            formula_str = re.sub('[+-]', '', formula_str_I);
            n12C = 0
            n13C = 0
            if 'C' not in Formula(formula_str)._elements: continue; #check if a carbon is even contained in the formula
            if 0 in Formula(formula_str)._elements['C']:
                n12C += Formula(formula_str)._elements['C'][0]; #get the # of Carbons
            if 13 in Formula(formula_str)._elements['C']:
                n13C += Formula(formula_str)._elements['C'][13]
            mnumber = Formula(formula_str).isotope.massnumber #get the nominal mass number
            spectrum = Formula(formula_str).spectrum() #get the spectrum
            fragmentSpectrum = {}
            intensityList = [];
            for c in range(-n13C, n12C + 1):
                if c<0:
                    fragmentSpectrum[Formula(formula_str).isotope.mass-1]=0.0;
                    intensityList.append(0.0);
                else:
                    if mnumber+c in spectrum:
                        fragmentSpectrum[spectrum[mnumber+c][0]]=spectrum[mnumber+c][1];
                        intensityList.append(spectrum[mnumber+c][1]);
                    else:
                        fragmentSpectrum[Formula(formula_str).isotope.mass + c]=0.0;
                        intensityList.append(0.0);
            fragmentSpectrum_tmp[formula_str_I] = fragmentSpectrum;

            # by default, the spectrum is normalized to the sum of all intensities measured
            # convert sum-normalized spectrum to max-normalized spectrum
            intensityListMax = max(intensityList);
            fragmentSpectrum = {};
            for k,v in fragmentSpectrum_tmp[formula_str_I].items():
                if round_mass:
                    fragmentSpectrum[int(numpy.round(k))] = v/intensityListMax;
                else:
                    fragmentSpectrum[k] = v/intensityListMax;
            fragmentSpectrum_O[formula_str_I] = fragmentSpectrum;

        return fragmentSpectrum_O;
    def report_fragmentSpectrum_normSum(self,fragments_I,round_mass=False):
        '''calculate the fragment spectrum'''
        # Input: formula_str_I
        # Output: spectrum_lst_O

        fragmentSpectrum_O = {};

        for formula_str_I in fragments_I:
            fragmentSpectrum_O[formula_str_I] = None;
            formula_str = re.sub('[+-]', '', formula_str_I);
            n12C = 0
            n13C = 0
            if 'C' not in Formula(formula_str)._elements: break; #check if a carbon is even contained in the formula
            if 0 in Formula(formula_str)._elements['C']:
                n12C += Formula(formula_str)._elements['C'][0]; #get the # of Carbons
            if 13 in Formula(formula_str)._elements['C']:
                n13C += Formula(formula_str)._elements['C'][13]
            mnumber = Formula(formula_str).isotope.massnumber #get the nominal mass number
            spectrum = Formula(formula_str).spectrum() #get the spectrum
            fragmentSpectrum = {}
            for c in range(-n13C, n12C + 1):
                if c<0:
                    exact_mass = Formula(formula_str).isotope.mass+c;
                    if round_mass:
                        fragmentSpectrum[int(numpy.round(exact_mass))]=0.0;
                    else:
                        fragmentSpectrum[exact_mass]=0.0;
                else:
                    if mnumber+c in spectrum:
                        exact_mass = spectrum[mnumber+c][0];
                        if round_mass:
                            fragmentSpectrum[int(numpy.round(exact_mass))]=spectrum[mnumber+c][1];
                        else:
                            fragmentSpectrum[exact_mass]=spectrum[mnumber+c][1];
                    else:
                        exact_mass = Formula(formula_str).isotope.mass + c
                        if round_mass:
                            fragmentSpectrum[int(numpy.round(exact_mass))]=0.0;
                        else:
                            fragmentSpectrum[exact_mass]=0.0;
            fragmentSpectrum_O[formula_str_I] = fragmentSpectrum;

        return fragmentSpectrum_O;
    def extract_peakData_normMax(self, peakData_I, fragments_I, res_I=0.3, round_mass=False):
        '''extract maximum intensity peak'''
        # Input: peakData_I = mass:intensity
        #        res_I = mass window/resolution (default = 0.3);
        # Output: 
        #   peakSpectrum_theoretical = {fragment:{mass:intensity}}
        #   peakSpectrum_measured = {fragment:{mass:intensity}} 
        #   peakSpectrum_corrected = {fragment:{mass:intensity}} 
        #   peakSpectrum_normalized = {fragment:{mass:intensity}} 

        '''The algorithm implement below does not track the peak width for calculation of peak area,
        nor for calculate of resolution using FWHM.  However, compared to peak-picking algorithm
        implemented in analyst(r) and peakView(r), the intensities for most compounds match
        the intensities calculated as peaks (compare 140228_MRM_EPI/..._EPI to ..._EPI_peakList
        or 140228_ER_EPI/...I to ..._ER).'''

        # min peak height
        detectionThreshold = 2500.0

        # pre-sort for efficiency
        # sort masses in peakData
        keys = list(peakData_I.keys());
        keys.sort();

        # determine baseline intensity
        # based on the most occuring intensity (background threshold);
        values = numpy.array(list(peakData_I.values()));
        values_median = mode(values)[0];
        if len(values_median) > 1:
            baseline = float(max(values_median)); # min returned too much junk
        else:
            baseline = float(values_median);
            
        if round_mass:
            peakSpectrum_theoretical = self.report_fragmentSpectrum_normMax(fragments_I,True);
        else:
            peakSpectrum_theoretical = self.report_fragmentSpectrum_normMax(fragments_I);
        # determine masses from fragments
        masses = [];
        peakSpectrum_measured_qcqa = {};
        peakSpectrum_normalized_qcqa = {};
        peakSpectrum_corrected_qcqa = {};
        peakSpectrum_measured = {};
        peakSpectrum_normalized = {};
        peakSpectrum_corrected = {};
        for frag,spec in peakSpectrum_theoretical.items():
            peakSpectrum_measured_qcqa[frag] = None;
            peakSpectrum_corrected_qcqa[frag] = None;
            peakSpectrum_normalized_qcqa[frag] = None;
            peakSpectrum_measured[frag] = None;
            peakSpectrum_corrected[frag] = None;
            peakSpectrum_normalized[frag] = None;

            if not spec: continue; #check if a carbon is even contained in the fragment

            masses = list(spec.keys());
            masses.sort(); # sort mass in massList

            keyIndex = 0;
            keyMax = len(keys);
            
            measured_qcqa = {};
            measured = {};
            for mass in masses: # iterate through each mass
                maxPeak = 0.0;
                keyMaxPeak = None;
                measured_qcqa[mass] = [keyMaxPeak,maxPeak];
                measured[mass] = maxPeak;
                while keyIndex<keyMax:
                    if keys[keyIndex] >= mass - res_I and keys[keyIndex] < mass + res_I:
                        peak = peakData_I[keys[keyIndex]];
                        if peak > maxPeak: 
                            maxPeak = peak;
                            keyMaxPeak = keys[keyIndex];
                        keyIndex += 1;
                    elif keys[keyIndex] < mass - res_I:
                        keyIndex += 1;
                        continue;
                    elif keys[keyIndex] >= mass + res_I:
                        measured_qcqa[mass] = [keyMaxPeak,maxPeak];
                        measured[mass] = maxPeak;
                        break;
            if measured: 
                peakSpectrum_measured_qcqa[frag] = measured_qcqa;
                peakSpectrum_measured[frag] = measured;
            else: break #no peaks were found for the fragment

            # correct intensity for background:
            corrected_qcqa = {};
            #intensityList = [];
            for k,v in peakSpectrum_measured_qcqa[frag].items():
                if v[1] > detectionThreshold:
                    if v[1] - baseline > 0.0: 
                        corrected_qcqa[k] = [v[0],v[1] - baseline];
                    else:
                        corrected_qcqa[k] = [v[0],0.0];
                else:
                    corrected_qcqa[k] = [v[0],0.0];
                #intensityList.append(corrected_qcqa[k][1]);
            peakSpectrum_corrected_qcqa[frag] = corrected_qcqa

            corrected = {};
            intensityList = [];
            for k,v in peakSpectrum_measured[frag].items():
                if v > detectionThreshold:
                    if v - baseline > 0.0: 
                        corrected[k] = v - baseline;
                    else:
                        corrected[k] = 0.0;
                    intensityList.append(corrected[k]);
                else:
                    corrected[k] = 0.0;
                intensityList.append(corrected[k]);
            peakSpectrum_corrected[frag] = corrected;

            # normalize each spectrum:
            normalized_qcqa = {};
            intensityListMax_qcqa = max(intensityList);
            for k,v in peakSpectrum_corrected_qcqa[frag].items():
                if intensityListMax_qcqa != 0: normalized_qcqa[k] = [v[0],v[1]/intensityListMax_qcqa];
                else: normalized_qcqa[k] = [v[0], None];
            peakSpectrum_normalized_qcqa[frag] = normalized_qcqa;
            
            normalized = {};
            intensityListMax = max(intensityList);
            for k,v in peakSpectrum_corrected[frag].items():
                if intensityListMax != 0: normalized[k] = v/intensityListMax;
                else: normalized[k] = None;
            peakSpectrum_normalized[frag] = normalized;

        return peakSpectrum_measured, peakSpectrum_corrected, peakSpectrum_normalized;
    def extract_peakData_normSum(self, peakData_I, fragments_I, res_I=0.3,round_mass=False):
        '''extract maximum intensity peak'''
        # Input: peakData_I = mass:intensity
        #        res_I = mass window/resolution (default = 0.3);
        # Output: 
        #   peakSpectrum_theoretical = {fragment:{mass:intensity}}
        #   peakSpectrum_measured = {fragment:{mass:intensity}} 
        #   peakSpectrum_corrected = {fragment:{mass:intensity}} 
        #   peakSpectrum_normalized = {fragment:{mass:intensity}} 

        # min peak height
        detectionThreshold = 1000.0

        # pre-sort for efficiency
        # sort masses in peakData
        keys = list(peakData_I.keys());
        keys.sort();

        # determine baseline intensity
        # based on the most occuring intensity (background threshold);
        values = numpy.array(list(peakData_I.values()));
        values_median = mode(values)[0];
        if len(values_median) > 1:
            baseline = float(max(values_median)); # min returned too much junk
        else:
            baseline = float(values_median);
            

        if round_mass:
            peakSpectrum_theoretical = self.report_fragmentSpectrum_normMax(fragments_I,True);
        else:
            peakSpectrum_theoretical = self.report_fragmentSpectrum_normMax(fragments_I);
        # determine masses from fragments
        masses = [];
        peakSpectrum_measured_qcqa = {};
        peakSpectrum_normalized_qcqa = {};
        peakSpectrum_corrected_qcqa = {};
        peakSpectrum_measured = {};
        peakSpectrum_normalized = {};
        peakSpectrum_corrected = {};
        for frag,spec in peakSpectrum_theoretical.items():
            peakSpectrum_measured_qcqa[frag] = None;
            peakSpectrum_corrected_qcqa[frag] = None;
            peakSpectrum_normalized_qcqa[frag] = None;
            peakSpectrum_measured[frag] = None;
            peakSpectrum_corrected[frag] = None;
            peakSpectrum_normalized[frag] = None;

            if not spec: continue; #check if a carbon is even contained in the fragment

            masses = list(spec.keys());
            masses.sort(); # sort mass in massList

            keyIndex = 0;
            keyMax = len(keys);
            
            measured_qcqa = {};
            measured = {};
            for mass in masses: # iterate through each mass
                maxPeak = 0.0;
                keyMaxPeak = None;
                measured_qcqa[mass] = [keyMaxPeak,maxPeak];
                measured[mass] = maxPeak;
                while keyIndex<keyMax:
                    if keys[keyIndex] >= mass - res_I and keys[keyIndex] < mass + res_I:
                        peak = peakData_I[keys[keyIndex]];
                        if peak > maxPeak: 
                            maxPeak = peak;
                            keyMaxPeak = keys[keyIndex];
                        keyIndex += 1;
                    elif keys[keyIndex] < mass - res_I:
                        keyIndex += 1;
                        continue;
                    elif keys[keyIndex] >= mass + res_I:
                        measured_qcqa[mass] = [keyMaxPeak,maxPeak];
                        measured[mass] = maxPeak;
                        break;
            if measured: 
                peakSpectrum_measured_qcqa[frag] = measured_qcqa;
                peakSpectrum_measured[frag] = measured;
            else: break #no peaks were found for the fragment

            # correct intensity for background:
            corrected_qcqa = {};
            #intensityList = [];
            for k,v in peakSpectrum_measured_qcqa[frag].items():
                if v[1] > detectionThreshold:
                    if v[1] - baseline > 0.0: 
                        corrected_qcqa[k] = [v[0],v[1] - baseline];
                    else:
                        corrected_qcqa[k] = [v[0],0.0];
                else:
                    corrected_qcqa[k] = [v[0],0.0];
                #intensityList.append(corrected_qcqa[k][1]);
            peakSpectrum_corrected_qcqa[frag] = corrected_qcqa

            corrected = {};
            intensityList = [];
            for k,v in peakSpectrum_measured[frag].items():
                if v > detectionThreshold:
                    if v - baseline > 0.0: 
                        corrected[k] = v - baseline;
                    else:
                        corrected[k] = 0.0;
                    intensityList.append(corrected[k]);
                else:
                    corrected[k] = 0.0;
                intensityList.append(corrected[k]);
            peakSpectrum_corrected[frag] = corrected;

            # normalize each spectrum:
            normalized_qcqa = {};
            intensityListSum_qcqa = sum(intensityList);
            for k,v in peakSpectrum_corrected_qcqa[frag].items():
                if intensityListSum_qcqa != 0: normalized_qcqa[k] = [v[0],v[1]/intensityListSum_qcqa];
                else: normalized_qcqa[k] = [v[0], None];
            peakSpectrum_normalized_qcqa[frag] = normalized_qcqa;
            
            normalized = {};
            intensityListSum = sum(intensityList);
            for k,v in peakSpectrum_corrected[frag].items():
                if intensityListSum != 0: normalized[k] = v/intensityListSum;
                else: normalized[k] = None;
            peakSpectrum_normalized[frag] = normalized;

        return peakSpectrum_measured, peakSpectrum_corrected, peakSpectrum_normalized;
    def extract_peakList_normMax(self, peakSpectrum_I, fragments_I, round_mass=False):
        '''extract peak spectrum from peak list'''
        # Input:
        #   peakSpectrum_I = {fragment:{(precursor_mass,product_mass):intensity}} 
        #   fragments_I = [fragments] 
        # Output: 
        #   peakSpectrum_corrected = {fragment:{mass:intensity}} 
        #   peakSpectrum_normalized = {fragment:{mass:intensity}} 
        
        # round all precursor/product masses in input for comparison:
        peakSpectrum_copy_I = {};
        for frag,spec in peakSpectrum_I.items():
            peakSpectrum_tmp = {};
            for masses,intensity in spec.items():
                peakSpectrum_tmp[numpy.around(masses)] = intensity;
            peakSpectrum_copy_I[frag] = peakSpectrum_tmp;

        if round_mass:
            peakSpectrum_theoretical = self.report_fragmentSpectrum_normMax(fragments_I,True);
        else:
            peakSpectrum_theoretical = self.report_fragmentSpectrum_normMax(fragments_I);
        # determine masses from fragments
        masses = [];
        peakSpectrum_normalized = {};
        peakSpectrum_corrected = {};
        for frag,spec in peakSpectrum_theoretical.items():
            peakSpectrum_corrected[frag] = None;
            peakSpectrum_normalized[frag] = None;

            if not spec: continue; #check if a carbon is even contained in the fragment

            masses = list(spec.keys());
            masses.sort(); # sort mass in massList
            masses_rounded = numpy.around(masses); # round masses to nearest digit for comparison


            # 1. copy data from peakSpectrum_I to peakSpectrum_corrected based on theoretical fragments
            intensityList = [];
            if frag in peakSpectrum_I:
                fragment_masses = [k for k in peakSpectrum_copy_I[frag].keys()];
                corrected_spec = {};
                for i,mass in enumerate(masses_rounded):
                    corrected = 0.0;
                    if mass in fragment_masses:
                        corrected = peakSpectrum_copy_I[frag][mass];
                    corrected_spec[masses[i]] = corrected;
                    intensityList.append(corrected);
                peakSpectrum_corrected[frag] = corrected_spec;
            else: 
                corrected_spec = {};
                for i,mass in enumerate(masses_rounded):
                    corrected = 0.0;
                    corrected_spec[masses[i]] = corrected;
                    intensityList.append(corrected);
                peakSpectrum_corrected[frag] = corrected_spec;


            # normalize each spectrum:
            #NOTE: normalization by max to allow for later conversion to normalization by sum
            normalized = {};
            intensityListMax = max(intensityList);
            for k,v in peakSpectrum_corrected[frag].items():
                if v:
                    if intensityListMax != 0: normalized[k] = v/intensityListMax;
                    else: normalized[k] = None;
                else: normalized[k] = None;
            peakSpectrum_normalized[frag] = normalized;

        return peakSpectrum_corrected, peakSpectrum_normalized;
    def extract_peakList_normSum(self, peakSpectrum_I, fragments_I, round_mass=False):
        '''extract peak spectrum from peak list'''
        # Input:
        #   peakSpectrum_I = {fragment:{mass:intensity}} 
        #   fragments_I = [fragments] 
        # Output: 
        #   peakSpectrum_corrected = {fragment:{mass:intensity}} 
        #   peakSpectrum_normalized = {fragment:{mass:intensity}} 
        
        # round all precursor/product masses in input for comparison:
        peakSpectrum_copy_I = {};
        for frag,spec in peakSpectrum_I.items():
            peakSpectrum_tmp = {};
            for masses,intensity in spec.items():
                peakSpectrum_tmp[numpy.around(masses)] = intensity;
            peakSpectrum_copy_I[frag] = peakSpectrum_tmp;

        if round_mass:
            peakSpectrum_theoretical = self.report_fragmentSpectrum_normSum(fragments_I,True);
        else:
            peakSpectrum_theoretical = self.report_fragmentSpectrum_normSum(fragments_I);
        # determine masses from fragments
        masses = [];
        peakSpectrum_normalized = {};
        peakSpectrum_corrected = {};
        for frag,spec in peakSpectrum_theoretical.items():
            peakSpectrum_corrected[frag] = None;
            peakSpectrum_normalized[frag] = None;

            if not spec: continue; #check if a carbon is even contained in the fragment

            masses = list(spec.keys());
            masses.sort(); # sort mass in massList
            masses_rounded = numpy.around(masses); # round masses to nearest digit for comparison


            # 1. copy data from peakSpectrum_I to peakSpectrum_corrected based on theoretical fragments
            intensityList = [];
            if frag in peakSpectrum_I:
                fragment_masses = [k for k in peakSpectrum_copy_I[frag].keys()];
                corrected_spec = {};
                for i,mass in enumerate(masses_rounded):
                    corrected = 0.0;
                    if mass in fragment_masses and peakSpectrum_copy_I[frag][mass]:
                        corrected = peakSpectrum_copy_I[frag][mass];
                    corrected_spec[masses[i]] = corrected;
                    intensityList.append(corrected);
                peakSpectrum_corrected[frag] = corrected_spec;
            else: 
                corrected_spec = {};
                for i,mass in enumerate(masses_rounded):
                    corrected = 0.0;
                    corrected_spec[masses[i]] = corrected;
                    intensityList.append(corrected);
                peakSpectrum_corrected[frag] = corrected_spec;


            # normalize each spectrum:
            normalized = {};
            intensityListSum = sum(intensityList);
            for k,v in peakSpectrum_corrected[frag].items():
                if v>0.0:
                    if intensityListSum != 0: normalized[k] = v/intensityListSum;
                    else: normalized[k] = None;
                else: normalized[k] = None;
            peakSpectrum_normalized[frag] = normalized;

        return peakSpectrum_corrected, peakSpectrum_normalized;
    def recombine_dilutionsMRMs(self,peakData_I):
        '''Method to "recombine" MRMs from one dilution to the next'''

        # input: peakData_I = {frag:[mass:{'intensity':intensity,
        #							    'dilution':dilution,
        #							    'used_':used_,
        #                               'comment_':comment_}]}
        # e.g.: {frag:[100:{'dilution':'high',...}],
        #             [101:{'dilution':'low','comment_':'Recombine',...}],
        #             [101:{'dilution':'high','comment_':'Recombine',...}],
        #             [102:{'dilution':'low','comment_':'Recombine',...}],
        #             [103:{'dilution':'low',...}],...}
        # NOTE: dictionary > List of dictionaries
        # NOTE: input list of masses must be sorted in ascending order
        #                followed by 'dilutions' in descending order as shown below!
        # output: peakData_O = {frag:{mass:{'intensity':intensity,
        #							    'dilution':dilution,
        #							    'used_':used_,
        #                               'comment_':comment_}}}
        #         peakData_O_false = {frag:{mass:{'intensity':intensity,
        #							    'dilution':dilution,
        #							    'used_':used_,
        #                               'comment_':comment_}}}
        # Note: second output structure needed to update rows that are changed to false

        '''Algorithm:
        start:
	    dilution	m	comment	used
	    'low'	0	''	false
	    'high'	0	''	true
	    'low'	1	'Recombine'	true
	    'high'	1	'Recombine'	true
	    'low'	2	'Recombine'	true
	    'high'	2	''	false
	    'low'	3	''	true
	    'high'	3	''	false
	    recombine...
	    end:
	    dilution	m	comment	used
	    'low'	0	''	false
	    'high'	0	''	true
	    'low'	1	'Recombine'	false
	    'high'	1	'Recombine'	true
	    'low'	2	'Recombine'	true
	    'high'	2	''	false
	    'low'	3	''	true
	    'high'	3	''	false
	    ...
	    done prior: set normalized intensity to diluion 'low', m 1 to 1;
	                recalculate the rest of the normalized intensities for the dilutions 'low', m 2,3,4,...;
	    calculate the percent change from dilution 'low', m 1 to dilution 'low', m 2; from dilution 'low', m 2 to dilution 'low', m 3; ...;
	    replace dilution 'high', m 2 with the normalized intensity for dilution 'low', m 1 - the percent change from dilution 'low', m 1 to dilution 'low', m 2;
		    replace dilution 'low', m 3 with the new normalized intensity for m 2 - the percent change from dilution 'low', m 2 to dilution 'low', m 3;
		    ...;'''
        
        peakData_O = {};
        peakData_O_false = {};
        #iterate through each fragment
        for frag,spec in peakData_I.items():
            peakData_O[frag] = None;
            peakData_O_false[frag] = None;
            spec_O = {};
            spec_O_false = {};
            if not spec: continue; #check if there is data for the fragment
            # extract out dilutions
            dilutions = [];
            for d in spec:
                values = list(d.values())[0];
                dilutions.append(values['dilution']);
            dilutions = list(set(dilutions));
            dilutions.sort();
            dilutions_dict = dict(list(zip(dilutions,['low','high'])));
            #iterate through each spectrum
            intensity_prev = 0.0
            intensity_new = 0.0;
            intensity_difference = 0.0;
            recombine_cnt = 0;
            for spec_dict in spec:
                mass = list(spec_dict.keys())[0];
                data = list(spec_dict.values())[0];
                spec_O[mass] = None;
                data_O = {};
                if not data['intensity']: 
                    data_O['dilution'] = None;
                    data_O['intensity'] = None;
                    data_O['comment_'] = None;
                    data_O['used_'] = None;
                    spec_O[mass] = data_O;
                    continue;
                if data['comment_'] == 'Recombine':
                    if recombine_cnt == 0: # 1st recombination event
                        if dilutions_dict[data['dilution']] != 'low': print('bad input');
                        intensity_prev = data['intensity'];
                        data['used_'] = False;
                        # copy the data
                        data_O['dilution'] = data['dilution'];
                        data_O['intensity'] = data['intensity'];
                        data_O['comment_'] = data['comment_'];
                        data_O['used_'] = data['used_'];
                        spec_O_false[mass] = data_O;
                        recombine_cnt += 1;
                        continue
                    elif recombine_cnt == 1: # 2nd recombination event
                        if dilutions_dict[data['dilution']] != 'high': print('bad input');
                        intensity_new = data['intensity'];
                        recombine_cnt += 1;
                    elif recombine_cnt == 2: # 3rd recombination event
                        if dilutions_dict[data['dilution']] != 'low': print('bad input');
                        intensity_difference = data['intensity']/intensity_prev;
                        intensity_prev = data['intensity'];
                        intensity_new = intensity_new*intensity_difference;
                        data['intensity'] = intensity_new;
                        recombine_cnt += 1;
                elif recombine_cnt >= 3:
                    if dilutions_dict[data['dilution']] != 'low': print('bad input');
                    intensity_difference = data['intensity']/intensity_prev;
                    intensity_prev = data['intensity'];
                    intensity_new = intensity_new*intensity_difference;
                    data['intensity'] = intensity_new;
                    recombine_cnt += 1;
                # copy data
                data_O['dilution'] = data['dilution'];
                data_O['intensity'] = data['intensity'];
                data_O['comment_'] = data['comment_'];
                data_O['used_'] = data['used_'];
                spec_O[mass] = data_O;
            # copy spectrum
            peakData_O[frag] = spec_O
            peakData_O_false[frag] = spec_O_false
        #copy out the intensities without the comments
        peakData_intensities_O = {};
        for frag,spec in peakData_O.items():
            spec_tmp = {};
            for mass,v in spec.items():
                spec_tmp[mass]=v['intensity'];
            peakData_intensities_O[frag] = spec_tmp;
        return peakData_O,peakData_O_false,peakData_intensities_O;
    def normalize_peakSpectrum_normMax(self,peakSpectrum_I,scalingFactors_I):
        '''normalize peakSpectrum taken from different m+0, m+1, ... fragments
        using a reference scaling factor'''

        # Input:
        #   peakSpectrum_I = {precursor_fragment:{product_fragment:{product_mass:intensity}}}
        #   scalingFactors_I = {precursor_fragment:intensity}
        # Output: 
        #   peakSpectrum_normalized = {product_fragment:{mass:intensity}}

        '''Algorithm:
        part 1: scale
        for each precursor i:
            for each product j in precursor i:
                for each mass m in product j:
                    peakSpectrum[precursor_i][product_j][m]*scalingFactor[precursor_i]
        part 2: reduce:
        for each product j in all precursors:
            for each mass in product j:
                for each precursor i with product j:
                    peakSpectrum_O[product_j][m] += peakSpectrum[precursor_i][product_j][m]*scalingFactor[precursor_i]
        '''

        precursor_fragments_I = list(peakSpectrum_I.keys());
        precursorSpectrum_dict = {};
        product_fragments_all = [];
        product_mass_all = [];
        # iterate through each precursor fragment
        for precursor in precursor_fragments_I:
            product_fragments_I = list(peakSpectrum_I[precursor].keys());
            productSpectrum_dict = {};
            product_fragments_all.extend(product_fragments_I);
            # iterate through each product fragment
            for product in product_fragments_I:
                spectrum_dict = {};
                product_mass_dict = {};
                product_mass_tmp = [];
                # iterate through each mass
                for k,v in peakSpectrum_I[precursor][product].items():
                    if peakSpectrum_I[precursor][product][k]:
                        spectrum_dict[k] = peakSpectrum_I[precursor][product][k]*scalingFactors_I[precursor];
                    else:
                        spectrum_dict[k] = 0.0;
                    product_mass_tmp.append(k);
                productSpectrum_dict[product] = spectrum_dict;
                product_mass_dict[product] = product_mass_tmp;
                product_mass_all.append(product_mass_dict);
            precursorSpectrum_dict[precursor] = productSpectrum_dict

        # reduce product fragments list
        product_fragments_reduced = list(set(product_fragments_all));
        
        # reduce product masses
        product_mass_combined = {};
        product_mass_reduced = {};
        for product in product_fragments_all:
            product_mass_combined[product] = [];
            for product_mass in product_mass_all:
                if product in product_mass:
                    product_mass_combined[product].extend(product_mass[product]);
            product_mass_reduced[product] = list(set(product_mass_combined[product]));

        peakSpectrum_normalized_O = {};
        # iterate through all common product fragments
        for product in product_fragments_reduced:
            peakSpectrum_normalized_O[product] = None;
            peakSpectrum_normalized_tmp = {};
            # iterate through each mass
            for mass in product_mass_reduced[product]:
                peakSpectrum_normalized_tmp[mass] = 0.0;
                # iterate through each precursor
                for precursor in precursor_fragments_I:
                    if product in precursorSpectrum_dict[precursor]:
                        if mass in precursorSpectrum_dict[precursor][product]:
                            peakSpectrum_normalized_tmp[mass] += precursorSpectrum_dict[precursor][product][mass]
                        else:
                            peakSpectrum_normalized_tmp[mass] += 0.0;
                    else: peakSpectrum_normalized_tmp[mass] += 0.0;
            peakSpectrum_normalized_O[product] = peakSpectrum_normalized_tmp;

        # re-normalize the spectrum to max-normalized spectrum
        intensityListMax = {};
        peakSpectrum_normalized_O_max = {};
        for product,spec in peakSpectrum_normalized_O.items():
             intensityList = [];
             for mass,intensity in spec.items():
                 intensityList.append(intensity);
             intensityListMax = max(intensityList);
             fragmentSpectrum = {};
             for mass,intensity in spec.items():
                 if intensityListMax != 0.0:
                    fragmentSpectrum[mass] = intensity/intensityListMax;
                 else:
                    fragmentSpectrum[mass] = 0.0;
             peakSpectrum_normalized_O_max[product] = fragmentSpectrum;

        return peakSpectrum_normalized_O_max
    def calculate_fragmentSpectrumAccuracy(self, peakSpectrum_normalized_list_I):
        '''calculate the accuracy from the normalized intensity'''
        # Input:
        #   peakSpectrum_normalized_list_I = [{fragment:{mass:intensity}}]
        # Output:
        #   peakSpectrum_accuracy_O = {fragment:float};

        fragments_I = list(peakSpectrum_normalized_list_I[0].keys());
        peakSpectrum_theoretical = self.report_fragmentSpectrum_normMax(fragments_I,True);

        peakSpectrum_accuracy_O = {};
        for frag in fragments_I:
            peakSpectrum_accuracy_O[frag] = None;

            if not peakSpectrum_theoretical[frag]: continue; # no carbons in fragment

            intensityList = [];
            masses = [];
            for peakSpectrum in peakSpectrum_normalized_list_I:
                intensityDict = {};
                peakSpectrumMasses = list(peakSpectrum_theoretical[frag].keys());
                for mass in peakSpectrumMasses:
                    if frag in peakSpectrum and mass in peakSpectrum[frag] and peakSpectrum[frag][mass] > 0.0: 
                        intensityDict[mass] = peakSpectrum[frag][mass];
                    else: 
                        intensityDict[mass] = 0.0;
                    if not mass in masses: masses.append(mass);
                intensityList.append(intensityDict);
                ## uncomment to only compare measured masses
                #intensityDict = {};
                #peakSpectrumMasses = peakSpectrum[frag].keys();
                #for mass in peakSpectrumMasses:
                #    if peakSpectrum[frag][mass] > 0.0: 
                #        intensityDict[mass] = peakSpectrum[frag][mass];
                #        if not mass in masses: masses.append(mass);
                #intensityList.append(intensityDict);
            accuracyLst = [];
            for mass in masses:
                data = [];
                for intensity in intensityList:
                    if intensity[mass]>=0.0:data.append(intensity[mass]);
                if data and peakSpectrum_theoretical[frag][mass]:
                    intensity_array = numpy.array(data);
                    accuracyLst.append(abs(intensity_array.mean() - peakSpectrum_theoretical[frag][mass]))

            accuracyLstMean = None;
            if accuracyLst: 
                accuracyLstMean = numpy.mean(accuracyLst);
                peakSpectrum_accuracy_O[frag] = accuracyLstMean;
            else: peakSpectrum_accuracy_O[frag] = None;

        return peakSpectrum_accuracy_O;
    def calculate_fragmentSpectrumAccuracy_normSum(self, peakSpectrum_normalized_list_I):
        '''calculate the accuracy from the normalized intensity'''
        # Input:
        #   peakSpectrum_normalized_list_I = [{fragment:{mass:intensity}}]
        # Output:
        #   peakSpectrum_accuracy_O = {fragment:float};

        fragments_I = list(peakSpectrum_normalized_list_I[0].keys());
        peakSpectrum_theoretical = self.report_fragmentSpectrum_normSum(fragments_I,True);

        peakSpectrum_accuracy_O = {};
        for frag in fragments_I:
            peakSpectrum_accuracy_O[frag] = None;

            if not peakSpectrum_theoretical[frag]: continue; # no carbons in fragment

            intensityList = [];
            masses = [];
            for peakSpectrum in peakSpectrum_normalized_list_I:
                intensityDict = {};
                peakSpectrumMasses = list(peakSpectrum_theoretical[frag].keys());
                for mass in peakSpectrumMasses:
                    if frag in peakSpectrum and mass in peakSpectrum[frag] and peakSpectrum[frag][mass] > 0.0: 
                        intensityDict[mass] = peakSpectrum[frag][mass];
                    else: 
                        intensityDict[mass] = 0.0;
                    if not mass in masses: masses.append(mass);
                intensityList.append(intensityDict);
                ## uncomment to only compare measured masses
                #intensityDict = {};
                #peakSpectrumMasses = peakSpectrum[frag].keys();
                #for mass in peakSpectrumMasses:
                #    if peakSpectrum[frag][mass] > 0.0: 
                #        intensityDict[mass] = peakSpectrum[frag][mass];
                #        if not mass in masses: masses.append(mass);
                #intensityList.append(intensityDict);
            accuracyLst = [];
            for mass in masses:
                data = [];
                for intensity in intensityList:
                    if intensity[mass]>=0.0:data.append(intensity[mass]);
                if data and peakSpectrum_theoretical[frag][mass]:
                    intensity_array = numpy.array(data);
                    accuracyLst.append(abs(intensity_array.mean() - peakSpectrum_theoretical[frag][mass]))

            accuracyLstMean = None;
            if accuracyLst: 
                accuracyLstMean = numpy.mean(accuracyLst);
                peakSpectrum_accuracy_O[frag] = accuracyLstMean;
            else: peakSpectrum_accuracy_O[frag] = None;

        return peakSpectrum_accuracy_O;
    def make_CSourceMix(self,csources_I, composition_I):
        '''Make a carbon source mix of a specified composition'''
        # Input: (e.g. 80/20 1-13C/U-13C glc)
        #       csources_I = backbone of the csources [['[13C]HO','CH2O','CH2O','CH2O','CH2O','CH3O'],
        #                                              ['[13C]HO','[13C]H2O','[13C]H2O','[13C]H2O','[13C]H2O','[13C]H3O']]
        #       composition_I = composition csources [0.8,0.2]
        # Output: 
        #       emu_O = {strings of emu distribution: spectral list}

        emu_O = {};
        emu_all = [];
        ncsources = len(csources_I)
        for cs in csources_I:
            emu_tmp = {};
            emu_tmp = self.make_EMUDistributionAndCSpectra(cs)
            emu_all.append(emu_tmp);
        for k in list(emu_all[0].keys()):
            spectra_tmp = [];
            spectra_tmp = [0.0]*len(emu_all[0][k])
            for i in range(ncsources):
                for j in range(len(emu_all[i][k])):
                    spectra_tmp[j] += composition_I[i]*emu_all[i][k][j];
            emu_O[k] = spectra_tmp;
        return emu_O; 
    def make_EMUDistributionAndCSpectra(self,csource_I):
        '''Make EMU distribution based on the carbon source'''
        # Input:
        #       csource_I = carbon backbone of the csource
        #                   e.g. 1-13C glc = ['[13C]HO','CH2','CH2','CH2','CH2','CH3O']
        #                        U-13C glc = ['[13C]HO','[13C]H2O','[13C]H2O','[13C]H2O','[13C]H2O','[13C]H3O']
        #                        glc = ['CHO','CH2O','CH2O','CH2O','CH2O','CH3O']
        # Output:
        #       emu_O = {strings of emu distribution: spectral list}

        nC = len(csource_I)
        emu_O = {};
        # iterate through each carbon and change from 0 to 1
        emu_c = nC*'0'; #intialize
        emu_lst = list(emu_c);
        for j in range(nC):
            emu_lst[j] = '1'
            for c in range(j,nC):
                emu_lst_2 = copy.copy(emu_lst)
                emu_lst_2[j] = '0';
                emu_lst_2[c] = '1';
                emu_tmp = copy.copy(emu_lst_2);
                cfrag = [];
                for i in range(c,nC):
                    emu_tmp[c] = '0';
                    emu_tmp[i] = '1';
                    emu_str = 'x' + ''.join(emu_tmp)
                    dfrag = [csource_I[p] for p,n in enumerate(emu_tmp) if n=='1']
                    dfrag_tmp = ''.join(dfrag)
                    #if emu_str.find('0')==-1: #ignore the fully labeled fragment
                    #    continue;
                    spectrum_tmp = self.report_fragmentSpectrum_normSum([dfrag_tmp],round_mass=True)
                    # format from dict into a list:
                    spectrum_tmp_lst = [];
                    spectrum_masses_lst = [];
                    for k,v in spectrum_tmp[dfrag_tmp].items():
                        spectrum_masses_lst.append(k);
                    spectrum_masses_lst.sort();
                    for k in spectrum_masses_lst:
                        spectrum_tmp_lst.append(spectrum_tmp[dfrag_tmp][k]);
                    emu_O[emu_str] = spectrum_tmp_lst;
        
        emu_c = nC*'1'; #intialize
        emu_lst = list(emu_c);
        for j in range(nC-1):
            emu_lst[j] = '0'
            for c in range(j,nC-1):
                emu_lst_2 = copy.copy(emu_lst)
                emu_lst_2[j] = '1';
                emu_lst_2[c] = '0';
                emu_tmp = copy.copy(emu_lst_2);
                cfrag = [];
                for i in range(c,nC-1):
                    emu_tmp[c] = '1';
                    emu_tmp[i] = '0';
                    emu_str = 'x' + ''.join(emu_tmp)
                    dfrag = [csource_I[p] for p,n in enumerate(emu_tmp) if n=='1']
                    dfrag_tmp = ''.join(dfrag)
                    #if emu_str.find('0')==-1: #ignore the fully labeled fragment
                    #    continue;
                    spectrum_tmp = self.report_fragmentSpectrum_normSum([dfrag_tmp],round_mass=True)
                    # format from dict into a list:
                    spectrum_tmp_lst = [];
                    spectrum_masses_lst = [];
                    for k,v in spectrum_tmp[dfrag_tmp].items():
                        spectrum_masses_lst.append(k);
                    spectrum_masses_lst.sort();
                    for k in spectrum_masses_lst:
                        spectrum_tmp_lst.append(spectrum_tmp[dfrag_tmp][k]);
                    emu_O[emu_str] = spectrum_tmp_lst;
        return emu_O;

    def make_fragmentID(self,met_id_I,formula_I,mass_I):
        """Make a unique fragment ID"""
        fragment_id_O = met_id_I + "_" + formula_I + "_" + str(mass_I);
        return fragment_id_O;

    def make_sampleFragmentID(self,sample_name_I,met_id_I,formula_I,mass_I):
        """Make a unique fragment ID"""
        fragment_id_O = sample_name_I + "_" + met_id_I + "_" + formula_I + "_" + str(mass_I);
        return fragment_id_O;
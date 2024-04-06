# Struct to capture/pass settings for run
import os
import ntpath

class SettingsContainer:
    def __init__(self, use_z, anno_only, source_neale, source_bolt, keep_tmp, species, vertebrate, plant, metazoa, fungi): #, bacteria, protist):
        self.use_z = use_z
        self.anno_only = anno_only
        self.source_neale = source_neale
        self.source_bolt = source_bolt
        self.keep_tmp = keep_tmp
        self.species = species
        self.vertebrate = vertebrate
        self.plant = plant
        self.metazoa = metazoa
        self.fungi = fungi
        # self.bacteria = bacteria
        # self.protist = protist
        speciesFlags = [self.vertebrate, self.plant, self.metazoa, self.fungi] #, self.bacteria, self.protist]
        if (self.species == "homo_sapiens") & (sum(speciesFlags) > 0):
            raise KeyError("May not use non-human descriptor flag (e.g. --plant, --bacteria, etc.) without providing --species.")
        if (self.species != "homo_sapiens") & (sum(speciesFlags) == 0):
            raise KeyError("Non-human species provided but non-human descriptor flag (e.g. --plant, --bacteria, etc.) not specified.")
        if sum(speciesFlags) > 1:
            raise KeyError("Only one of the non-human descriptor flags (e.g. --plant, --bacteria, etc.) is allowed for non-human species.")

        return

# Struct to concisely pass all the arguments/filenames/etc to the functions
class ArgumentContainer:
    def __init__(self, gwas, pval, p_col, maf, maf_col, ref_col, alt_col, snp_map, output, cache_dir, cache_version, biotypes, dist, eoi, ref, build, draws, settings):
        self.gwas = gwas
        self.pval = pval
        self.p_col = p_col
        self.maf = maf
        self.maf_col = maf_col
        self.ref_col = ref_col
        self.alt_col = alt_col
        self.snp_map = snp_map

        # Change/assign variables which depend on other settings/options
        if (settings.source_neale == True) & (settings.source_bolt == True):
            raise KeyError("Cannot use both --source_neale and --source_bolt. Please select only one or specify custom GWAS column names with --gwas_{HEADER} option. Exiting.")
        if settings.source_neale == True:
            self.p_col = "pval"
            self.maf_col = "minor_AF"
        if settings.source_bolt == True:
            self.p_col = "P_BOLT_LMM_INF"
            self.maf_col = "MAF"
        if settings.use_z == True:
            self.p_col = "Z"

        if output == None:
            self.output = os.path.splitext(ntpath.basename(gwas))[0]
        else:
            self.output = output

        self.cache_dir = cache_dir
        self.cache_v = cache_version
        self.biotypes = biotypes
        self.dist = dist
        self.eoi = eoi
        self.ref = ref
        self.draws = draws
        self.build = build
        # if settings.species != "homo_sapiens":
        #     if self.build == "38":
        #         print("WARNING: Using non-human --species, but --build specified as 38 (only applicable to human reference genome). Will ignore build specification.")

        return

class PrerankArgumentContainer:
    def __init__(self, input, output, ref_build, buffer, biotypes, topN, pval_col, chr_col, pos_col, pval, dmpN, dmpP, dmpD, chr, excl, incl, score_method, call_peaks):
        self.input = input
        self.output = output

        if output == None:
            self.output = os.path.splitext(ntpath.basename(input))[0]
        else:
            self.output = output

        self.build = ref_build
        self.buffer = buffer
        self.biotypes = biotypes
        self.nTop = topN
        self.vCol = pval_col
        self.cCol = chr_col
        self.pCol = pos_col
        self.pThresh = pval
        self.dmpThresh = dmpN
        self.dmpP = dmpP
        self.dmpD = dmpD

        if self.dmpD == None:
            self.dmpD = buffer

        self.chrReg = chr
        if self.chrReg != None:
            if self.chrReg not in range(1,23):
                raise KeyError('--chr must be selected from chromosome 1 to 22. Specifying multiple chromosomes and/or positions is not supported. Exiting.')

        self.xFile = excl
        self.fFile = incl

        self.scoring = score_method.lower()
        if self.scoring not in ['min', 'mean']:
            raise KeyError('--score_method input invalid. Options are \'min\' or \'mean\'. Exiting.')

        self.call_peaks = call_peaks

        return

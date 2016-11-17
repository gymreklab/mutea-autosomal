import collections
import genotypers
import numpy
import operator
import sys
import vcf

ERRPROB = 0.01 # Probability a read is an error (rather than stutter)

def ERROR(msg):
    sys.stderr.write(msg.strip() + "\n")
    sys.exit(1)

# Returns a dictionary mapping each sample name to a dictionary 
# containing the observed repeat lengths and counts for that sample's reads
def get_str_read_counts(vcf_reader, uselocus=None):
    # Attempt to read record from VCF
    if uselocus is not None:
        loci = list(vcf_reader.fetch(*uselocus))
        loci = [locus for locus in loci if locus.POS == uselocus[1]]
        if len(loci) != 1:
            return False, "N/A", -1, -1, -1, {}, -1, -1, None
        record = loci[0]
    else:
        try:
            record = vcf_reader.next()
        except StopIteration:
            return False, "N/A", -1, -1, -1, {}, -1, -1, None

    # Determine most common frame
    if 'MOTIF' in record.INFO:
        motif_len = len(record.INFO['MOTIF'])
    else:
        motif_len = record.INFO['PERIOD']
    frame_counts = collections.defaultdict(int)
    for sample in record:
        if sample['GT'] is not None:
            key = 'ALLREADS' if hasattr(sample.data, 'ALLREADS') else 'MALLREADS'
            if sample[key] is None:
                continue

            read_count_data  = map(lambda x: (int(x[0]), int(x[1])), map(lambda x: x.split("|"), sample[key].split(";")))
            tot_sample_reads = sum(map(lambda x: x[1], read_count_data))
            for str_length,read_count in read_count_data:
                frame_counts[str_length%motif_len] += read_count*1.0/tot_sample_reads 
    if len(frame_counts) == 0:
        return get_str_read_counts(vcf_reader)
    best_frame,value = max(frame_counts.items(), key=operator.itemgetter(1))

    # Construct a dictionary mapping sample names -> read count dict for each sample without out-of-frame reads
    sample_read_count_dict = {}
    in_frame_count,out_frame_count = 0,0
    for sample in record:
        if sample['GT'] is not None:
            if sample[key] is None:
                continue
            read_count_data = map(lambda x: (int(x[0]), int(x[1])), map(lambda x: x.split("|"), sample[key].split(";")))
            read_sizes      = map(lambda x: x[0], read_count_data)
            all_in_frame    = all(map(lambda x: x%motif_len == best_frame, read_sizes))
            if all_in_frame:
                # Convert bp diffs to repeat diffs
                read_count_data = map(lambda x: ((x[0]-best_frame)/motif_len, x[1]), read_count_data)

                sample_read_count_dict[sample.sample] = dict(read_count_data)
                in_frame_count += 1
            else:
                out_frame_count += 1
    str_key = record.ID if record.ID is not None else record.CHROM+":"+str(record.POS)
    return True,record.CHROM,record.POS,record.INFO['END'],motif_len,sample_read_count_dict, in_frame_count, out_frame_count, str_key

def get_sample_tmrcas(vcf_reader, uselocus=None):
    # Attempt to read record from VCF
    if uselocus is not None:
        loci = list(vcf_reader.fetch(*uselocus))
        loci = [locus for locus in loci if locus.POS == uselocus[1]]
        if len(loci) != 1:
            return {}
        record = loci[0]
    else:
        try:
            record = vcf_reader.next()
        except StopIteration:
            return {}
    tmrcas = {}
    for sample in record:
        try:
            tmrcas[sample.sample] = int(sample["TMRCA"])
        except: ERROR("No TMRCA field found")
    return tmrcas

def get_str_gts_diploid(vcf_reader, uselocus=None):
    # Attempt to read record from VCF
    if uselocus is not None:
        loci = list(vcf_reader.fetch(*uselocus))
        loci = [locus for locus in loci if locus.POS == uselocus[1]]
        if len(loci) != 1:
            return False, {}, -1, -1, None, None
        record = loci[0]
    else:
        try:
            record = vcf_reader.next()
        except StopIteration:
            return False, {}, -1, -1, None, None
    
    # Get all lengths
    all_lens = []
    for sample in record:
        if sample["GT"] is not None:
            for gt in map(int, sample["GT"].split("/")):
                all_lens.append(len(str(record.alleles[gt])))

    # Determine most common frame
    if 'MOTIF' in record.INFO:
        motif_len = len(record.INFO['MOTIF'])
    else:
        motif_len = record.INFO['PERIOD']
    all_lens     = sorted(all_lens)
    frame_counts = motif_len*[0]
    if len(all_lens) == 0:
        exit("No samples have a valid genotype for the STR of interest")
    for length in all_lens:
        frame_counts[length%motif_len] += 1
    best_frame, value = max(enumerate(frame_counts), key=operator.itemgetter(1))

    # Utilize median observed allele as "center"
    all_lens = filter(lambda x: x%motif_len == best_frame, all_lens)
    all_lens = sorted(all_lens)
    center   = all_lens[len(all_lens)/2]
    min_str  = (all_lens[0]  - center)/motif_len
    max_str  = (all_lens[-1] - center)/motif_len

    # Determine genotype posteriors, relative to "center" and in terms of # repeats
    gt_probs = {}
    for sample in record:
        sample_probs = {}
        if sample["GT"] is None: continue
        gt1, gt2 = map(int, sample["GT"].split("/"))
        samplegt = ( (len(record.alleles[gt1])-center)/motif_len, (len(record.alleles[gt2])-center)/motif_len )
        for i in xrange(len(record.alleles)):
            for j in xrange(len(record.alleles)):
                if record.alleles[i] is None or record.alleles[j] is None: continue
                if len(record.alleles[i])%motif_len != best_frame or len(record.alleles[j])%motif_len != best_frame: continue
                repeat_diffs = ( (len(record.alleles[i])-center)/motif_len, (len(record.alleles[j])-center)/motif_len)
                post = int(samplegt==repeat_diffs)
                if post == 0: continue # don't need to record
                if repeat_diffs in sample_probs:
                    sample_probs[repeat_diffs] += post
                else: sample_probs[repeat_diffs] = post
        gt_probs[sample.sample] = sample_probs
        
    if record.ID is not None:
        str_key = record.ID
    else:
        str_key = record.CHROM + ":" + str(record.POS)
    return True, gt_probs, min_str, max_str, str_key, motif_len

def get_str_gts(vcf_reader):
    # Attempt to read record from VCF
    try:
        record = vcf_reader.next()
    except StopIteration:
        return False, {}, 0, 0, ""

    # Ignore missing and heterozygous genotypes
    all_lens  = []
    for sample in record:
        if sample['GT'] is not None:
            all_lens.append(len(str(record.alleles[int(sample['GT'])])))

    # Determine most common frame
    if 'MOTIF' in record.INFO:
        motif_len = len(record.INFO['MOTIF'])
    else:
        motif_len = record.INFO['PERIOD']
    all_lens     = sorted(all_lens)
    frame_counts = motif_len*[0]
    if len(all_lens) == 0:
        exit("No samples have a valid genotype for the STR of interest")
    for length in all_lens:
        frame_counts[length%motif_len] += 1
    best_frame, value = max(enumerate(frame_counts), key=operator.itemgetter(1))
    
    # Utilize median observed allele as 'center'
    all_lens = filter(lambda x: x%motif_len == best_frame, all_lens)
    all_lens = sorted(all_lens)
    center   = all_lens[len(all_lens)/2]
    min_str  = (all_lens[0]  - center)/motif_len
    max_str  = (all_lens[-1] - center)/motif_len

    # Determine posterior genotype likelihoods for samples
    # Samples whose most probable genotype is missing, heterozygous or out-of-frame are not assigned any likelihoods
    # Lengths of genotypes are stored relative to the 'center' and in terms of the number of repeat differences
    gt_probs  = {}
    for sample in record:
        repeat_diffs = []
        diff_LLs     = []
        if sample['GT'] is None:
            continue
        for i in xrange(len(record.alleles)):
            if record.alleles[i] is None:
                continue
            if len(record.alleles[i])%motif_len != best_frame:
                continue
            repeat_diffs.append((len(record.alleles[i])-center)/motif_len)
            diff_LLs.append(sample['GL'][i])

        diff_LLs     = numpy.array(diff_LLs)
        max_LL       = numpy.max(diff_LLs)
        total_LL     = max_LL + numpy.log(numpy.sum(numpy.exp(diff_LLs-max_LL)))
        posteriors   = numpy.exp(diff_LLs - total_LL)
        sample_probs = {}

        # Combine posteriors for repeats with same length
        # Is there are better way than to just add them?
        for i in xrange(len(repeat_diffs)):
            if repeat_diffs[i] in sample_probs:
                sample_probs[repeat_diffs[i]] += posteriors[i]
            else:
                sample_probs[repeat_diffs[i]] = posteriors[i]
        gt_probs[sample.sample] = sample_probs

    if record.ID is not None:
        str_key = record.ID
    else:
        str_key = record.CHROM + ":" + str(record.POS)
    return True,gt_probs,min_str,max_str,str_key

def compute_median(values, counts):
    total_count = sum(counts)
    target_val  = total_count/2.0
    for i in xrange(len(values)):
        if counts[i] > target_val:
            return values[i]
        target_val -= counts[i]
    return values[-1]

def counts_to_centalized_posteriors(sample_read_counts, p_geom, down, up, diploid=False):
    # Compute the minimum and maximum observed allele
    min_allele = min(map(lambda x: min(x.keys()), sample_read_counts.values()))
    max_allele = max(map(lambda x: max(x.keys()), sample_read_counts.values()))

    # Create the stutter genotyper
    genotyper = genotypers.EstStutterGenotyper(_diploid=diploid)
    genotyper.create(down, up, p_geom, min_allele, max_allele)

    # Compute genotype posteriors using genotyper
    gt_posteriors = {}
    for sample,counts in sample_read_counts.items():
        gt_posteriors[sample] = genotyper.get_genotype_posteriors(counts, errprob=ERRPROB)

    # Compute 'central' allele using the allele with the median posterior sum
    gt_counts = collections.defaultdict(int)
    for posteriors in gt_posteriors.values():
        for gt,prob in  posteriors.items():
            if diploid:
                gt_counts[gt[0]] += prob
                gt_counts[gt[1]] += prob
            else: gt_counts[gt] += prob
    count_items = sorted(gt_counts.items())
    center      = compute_median(map(lambda x: x[0], count_items), map(lambda x: x[1], count_items))
    print("CENTRAL ALLELE = %d"%(center))

    # Normalize all genotypes relative to this central allele
    norm_gt_posteriors = {}
    for sample,posteriors in gt_posteriors.items():
        new_posteriors = {}
        for gt,prob in posteriors.items():
            if diploid:
                new_posteriors[(gt[0]-center, gt[1]-center)] = prob
            else: new_posteriors[gt-center] = prob
        norm_gt_posteriors[sample] = new_posteriors
    if diploid:
        return norm_gt_posteriors,min_allele-center,max_allele-center,center
    else: return norm_gt_posteriors,min_allele-center,max_allele-center


def patient_negative_for(patient, variant):
    b = True
    variants = []
    variants.extend(patient['snv'])
    variants.extend(patient['cnv'])
    variants.extend(patient['io'])
    for v in variants:
        # if variant == v['gene'] + ' ' + v['pdot']:
        if variant.startswith(v['gene']):
            b = False
            break
    return b

def get_gene_from_variant(variant):
    return variant.split(' ', 1)[0]


def append_variant(variants, variant):
    gene = variant.split(' ', 1)[0]
    pdot = variant.split(' ', 1)[1]
    if gene=='KRAS' or gene=='NRAS':
        pdot = 'mutation'
    found = False
    for v in variants:
        if v['gene']==gene:
            if pdot not in v['pdots']:
                v['pdots'].append(pdot)
            found = True
            break
    if not found:
        v = {'gene':gene, 'pdots':[pdot]}
        variants.append(v)


def transform_variants(variants):
    l = []
    for v in variants:
        s = ''
        for pdot in v['pdots']:
            if len(s)>0:
                s = s + ', '
            s = s + pdot
        s = v['gene'] + ' ' + s
        l.append(s)
    return l



def get_pertinent_negatives_for_report(patient,evidence_collection):
    variants = []
    myquery = {'approval_index': {'$lte': 3}, 'diseases': patient['ckb_disease']}
    mydocs = evidence_collection.find(myquery)
    for doc in mydocs:
        for variant in doc['variants']:
            if 'CD274' in variant:
                variant = variant.replace('CD274', 'PD-L1')
            if 'amp' in variant:
                variant = variant.replace('amp', 'amplification')
            if 'ins' in variant:
                variant = variant.replace('ins', 'insertion')
            if 'del' in variant:
                variant = variant.replace('del', 'deletion')
            if 'rearrange' in variant:
                variant = variant.replace('rearrange', 'fusion')
            if 'mutant' in variant:
                variant = variant.replace('mutant', 'mutation')
            if 'act mut' in variant:
                variant = variant.replace('act mut', 'mutation')
            if 'inmutation' in variant:
                variant = variant.replace('inmutation', 'mutation')
            if variant == 'ALK positive':
                variant = 'ALK fusion'
            if variant == 'ROS1 positive':
                variant = 'ROS1 fusion'
            if variant == 'RET positive':
                variant = 'RET fusion'
            if variant.startswith('BRAF V600'):
                variant = 'BRAF V600 mutation'
            if variant.startswith('EGFR'):
                variant = 'EGFR activating mutation'
            if 'NTRK' in variant:
                variant = 'NTRK1/2/3 fusions'
            if variant not in variants:
                if variant != 'Unknown unknown' and 'wild-type' not in variant:
                    if patient_negative_for(patient, variant):
                        append_variant(variants,variant)
    variants = transform_variants(variants)
    variants = sorted(variants)

    return variants

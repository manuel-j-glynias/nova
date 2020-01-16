
# genes:{$in: ['KRAS','NRAS']}
# {variants: {'$regex' :'wild-type$'},approval_index: {'$lte': 3}, diseases:'colorectal cancer' }
# {approval_index :1, diseases: {$in: ['stomach cancer', 'hepatocellular carcinoma']} }
import pprint


def get_diseases_query_terms(diseases):
    terms = []
    for d in diseases:
         # terms.append("'" + d + "'")
         terms.append(d)
    return terms


def get_wild_types_for_disease(patient,evidence_collection):
    myquery = {'approval_index': {'$lte': 3}, 'diseases': {'$in': get_diseases_query_terms(patient['ckb_diseases'])}, 'variants': {'$regex' :'wild-type$'}}
    mydocs = evidence_collection.find(myquery)
    variants = []
    for doc in mydocs:
        for variant in doc['variants']:
            if 'wild-type' in variant and variant not in variants:
                variants.append(variant)
    return variants


def patient_negative_for(patient, variant):
    b = True
    patient_variants = []
    patient_variants.extend(patient['snv'])
    patient_variants.extend(patient['cnv'])
    for v in patient_variants:
        # if variant == v['gene'] + ' ' + v['pdot']:
        if variant.startswith(v['gene']):
            b = False
            break
    return b


def patient_positive_for(patient, variant):
    for v in patient['reportable_markers']:
        if variant == v['gene'] + ' ' + v['pdot']:
            return True, v
    return False,None



def add_WT_to_patient(patient,wt):
    for variant in wt:
        if patient_negative_for(patient,variant):
            pos = variant.find('wild')
            gene = variant[:pos-1]
            wt_variant = {'gene':gene, 'pdot':'wild-type', 'type':'snv'}
            patient['reportable_markers'].append(wt_variant)



def get_variants_for_evidence_query(patient):
    variants = []
    variant_dict = {}
    for v in patient['reportable_markers']:
        variant = v['gene'] + ' ' + v['pdot']
        variants.append(variant)
        variant_dict[variant] = v
        if v['type'] == 'io':
            if v['gene']== 'PD-L1':
                variant = 'CD274 ' +  v['pdot']
                variants.append(variant)
                variant_dict[variant] = v
        else:
            if 'is_oncogene' in v:
                if v['is_oncogene'] and v['is_activating']:
                    variant = v['gene'] + ' act mut'
                    variants.append(variant)
                    variant_dict[variant] = v

                    # hack for stupid jax labeling issue
                    variant = v['gene'] + ' inact mut'
                    variants.append(variant)
                    variant_dict[variant] = v

                    variant = v['gene'] + ' mutant'
                    variants.append(variant)
                    variant_dict[variant] = v
                else:
                    variant = v['gene'] + ' inact mut'
                    variants.append(variant)
                    variant_dict[variant] = v

                    # hack for stupid jax labeling issue
                    variant = v['gene'] + ' act mut'
                    variants.append(variant)
                    variant_dict[variant] = v

                    variant = v['gene'] + ' mutant'
                    variants.append(variant)
                    variant_dict[variant] = v
    patient['ckb_variants'] = variants
    return variants, variant_dict


def has_variant(v, indication):
    b = False
    for variant in indication['variants']:
        if v['gene']== variant['gene'] and v['pdot']==variant['pdot']:
            b = True
            break
    return b

def maybe_add_variant(evidence, indication):
    if evidence['indication']==indication['indication']:
        for v in evidence['variants']:
            if not has_variant(v,indication):
                indication['variants'].append(v)


def maybe_append_evidence(indications, evidence):
    therapy = evidence['therapy']
    b = True
    for indication in indications:
        if indication['therapy']==therapy:
            b = False
            maybe_add_variant(evidence,indication)
            break
    if b:
        indications.append(evidence)


def get_evdience_for_disease_and_markers(patient,evidence_collection,approval_index,diseases):
    variants,variant_dict  = get_variants_for_evidence_query(patient)

    myquery = {'approval_index': {'$lte': approval_index}, 'diseases': {'$in': get_diseases_query_terms(diseases)},
               'variants': {'$in' : variants} }
    mydocs = evidence_collection.find(myquery).sort("approval_index", 1)
    for doc in mydocs:
        evidence = get_evidence_from_doc(doc, variant_dict)
        evidence['disease'] = evidence['indication']
        evidence['off_label'] = False
        if doc['contraindicated']:
            evidence['contraindicated'] = True
            patient['contraindications'].append(evidence)
        else:
            evidence['contraindicated'] = False
            maybe_append_evidence(patient['indications'],evidence)


def therapy_not_on_label(therapy, indications):
    b = True
    for indication in indications:
        if indication['therapy']==therapy:
            b = False
            break
    return b


# def get_evdience_for_markers(patient,evidence_collection,approval_index):
#     patient['off_label'] = []
#     variants,variant_dict  = get_variants_for_evidence_query(patient)
#
#     myquery = {'approval_index': {'$lte': approval_index},
#                'variants': {'$in' : variants} }
#     mydocs = evidence_collection.find(myquery).sort("approval_index", 1)
#     for doc in mydocs:
#         evidence = get_evidence_from_doc(doc, variant_dict)
#     if therapy_not_on_label(evidence['therapy'],patient['indications']):
#         evidence['off_label'] = True
#         evidence['contraindicated'] = False
#         patient['off_label'].append(evidence)


def get_evidence_from_doc(doc, variant_dict):
    evidence = {}
    evidence['therapy'] = doc['therapy']
    evidence['indication'] = doc['indication']
    evidence['drugs'] = doc['drugs']
    evidence['drug_tree'] = doc['drug_tree']
    evidence['evidence_type'] = doc['evidence_type']
    evidence['efficacy_evidence'] = doc['efficacy_evidence']
    evidence['approval_status'] = doc['approval_status']
    evidence['approval_index'] = doc['approval_index']
    evidence['ampCapAscoEvidenceLevel'] = doc['ampCapAscoEvidenceLevel']
    evidence['ampCapAscoInferredTier'] = doc['ampCapAscoInferredTier']
    evidence['cap'] = evidence['ampCapAscoInferredTier'] + evidence['ampCapAscoEvidenceLevel']
    evidence['response_type'] = doc['response_type']
    evidence['availability'] = '-'
    evidence['io'] = False
    evidence['variants'] = []
    variants = doc['variants']
    for variant in variants:
        # b,v = patient_positive_for(patient, variant)
        if variant in variant_dict:
            v = variant_dict[variant]
            evidence['variants'].append(v)
            if v['type'] == 'io':
                evidence['io'] = True
    return evidence


# approval_index:5 = Phase II
def get_evdience_for_breast_cancer(patient, evidence_collection):
    diseases = ['Her2-receptor positive breast cancer',
                'Her2-receptor negative breast cancer',
                'triple-receptor negative breast cancer',
                'estrogen-receptor positive breast cancer',
                'progesterone-receptor positive breast cancer']

    get_evdience_for_disease_and_markers(patient, evidence_collection, 5, patient['ckb_diseases'])
    get_evdience_for_disease_and_markers(patient, evidence_collection, 5, diseases)



def find_therapies(patient, evidence_collection):
    wt = get_wild_types_for_disease(patient, evidence_collection)
    add_WT_to_patient(patient, wt)
    patient['indications'] = []
    patient['contraindications'] = []
    patient['cdx_io'] = []


    if patient['OmniDisease']=='Breast Cancer':
        get_evdience_for_breast_cancer(patient,evidence_collection)
    else:
        get_evdience_for_disease_and_markers(patient, evidence_collection,5,patient['ckb_diseases'])

    # get_evdience_for_markers(patient, evidence_collection, 5)   no off-label

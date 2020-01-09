
# genes:{$in: ['KRAS','NRAS']}
# {variants: {'$regex' :'wild-type$'},approval_index: {'$lte': 3}, diseases:'colorectal cancer' }
import pprint


def get_wild_types_for_disease(patient,evidence_collection):
    myquery = {'approval_index': {'$lte': 3}, 'diseases': patient['ckb_disease'], 'variants': {'$regex' :'wild-type$'}}
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


def get_evdience_for_disease_and_markers(patient,evidence_collection,approval_index):
    variants,variant_dict  = get_variants_for_evidence_query(patient)

    myquery = {'approval_index': {'$lte': approval_index}, 'diseases': patient['ckb_disease'],
               'variants': {'$in' : variants} }
    mydocs = evidence_collection.find(myquery).sort("approval_index", 1)
    for doc in mydocs:
        evidence = get_evidence_from_doc(doc, variant_dict)
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


def get_evdience_for_markers(patient,evidence_collection,approval_index):
    patient['off_label'] = []
    variants,variant_dict  = get_variants_for_evidence_query(patient)

    myquery = {'approval_index': {'$lte': approval_index},
               'variants': {'$in' : variants} }
    mydocs = evidence_collection.find(myquery).sort("approval_index", 1)
    for doc in mydocs:
        evidence = get_evidence_from_doc(doc, variant_dict)
    if therapy_not_on_label(evidence['therapy'],patient['indications']):
        evidence['off_label'] = True
        evidence['contraindicated'] = False
        patient['off_label'].append(evidence)


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
def find_therapies(patient, evidence_collection):
    wt = get_wild_types_for_disease(patient, evidence_collection)
    add_WT_to_patient(patient, wt)
    patient['indications'] = []
    patient['contraindications'] = []
    patient['cdx_io'] = []

    get_evdience_for_disease_and_markers(patient, evidence_collection,5)

    get_evdience_for_markers(patient, evidence_collection, 5)
    # print()
    # print(patient['ckb_disease'])
    # for v in patient['reportable_markers']:
    #     print(v['gene'] + ' ' + v['pdot'])
    #
    # for indication in patient['indications']:
    #     print(indication['therapy'],indication['ampCapAscoEvidenceLevel'],indication['ampCapAscoInferredTier'],
    #           indication['variants'])
    # print('off_label')
    # for indication in patient['off_label']:
    #     print(indication['therapy'],indication['ampCapAscoEvidenceLevel'],indication['ampCapAscoInferredTier'],
    #           indication['indication'],indication['variants'])



# def add_pathologists_comments(patient):
#     pathologist_comments = 'This 44-year old woman with lung non-small cell carcinoma has one genomic variant (EGFR ' \
#                            'L858R) and one positive immunotherapy marker (PD-L1) with companion diagnostic indications. Osimertinib, a third-generation EGFR inhibitor, is an approved frontline targeted therapy for this patient. Patients with EGFR mutations may also respond to single agent anti-PD-1 therapy following progression on EGFR therapy. The high tumor mutational burden in this patient may indicate response to combination ipilimumab + nivolumab or single agent nivolumab immunotherapy. This patient does have a mutation in STK11, and emerging evidence indicates they are less likely to respond to anti-PD-1 therapy. There are two additional genomic variants detected with therapeutic options in clinical trials (TP53 E62X, BRCA2 K585R), as well as highly expressed immuno-oncology markers with targets in clinical trials.'
#     patient['pathologist_comments'] = pathologist_comments